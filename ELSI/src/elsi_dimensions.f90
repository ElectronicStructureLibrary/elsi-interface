!Copyright (c) 2016, ELSI consortium
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
! * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer in the
!   documentation and/or other materials provided with the distribution.
! * Neither the name of the ELSI project nor the
!   names of its contributors may be used to endorse or promote products
!   derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!> 
!! This module contains variables accessible in ELSI and related modules.
!!

module ELSI_DIMENSIONS

  use iso_c_binding
  use f_ppexsi_interface

  implicit none
  
  !> Solver (ELPA=0,LIBOMM=1,PEXSI=2,CHESS=3)
  integer :: method = -1

  !> Matrix storage format
  character(5) :: matrix_format

  !> Global matrix dimension
  integer :: n_g_rank
  
  !> Block size
  integer :: n_b_rows
  integer :: n_b_cols

  !> Processor grid
  integer :: n_p_rows
  integer :: n_p_cols
  
  !> Local matrix dimension
  integer :: n_l_rows
  integer :: n_l_cols

  !> MPI variables
  integer :: myid                   !< Local process ID
  integer :: n_procs                !< Number of MPI tasks
  integer :: mpi_comm_global        !< Global MPI communicator
  integer :: mpierr                 !< Error handler
  logical :: mpi_is_setup = .false. !< Is MPI initialized? 

  !> BLACS variables
  integer :: blacs_ctxt               !< Blacs context
  integer :: sc_desc(9)               !< Blacs descriptor
  integer :: mpi_comm_row             !< Row communicatior for ELPA
  integer :: mpi_comm_col             !< Column communicatior for ELPA
  integer :: my_p_row                 !< Process row position
  integer :: my_p_col                 !< Process column position
  logical :: external_blacs = .true.  !< Is BLACS initialized outside ELSI?
  logical :: blacs_is_setup = .false. !< Is BLACS initialized?

  !> Sparse matrix information
  integer :: n_g_nonzero !< Number of nonzeros in the global matrix
  integer :: n_l_nonzero !< Number of nonzeros in the local matrix

  !> Overlap
  logical :: overlap_is_unity = .false. !< Is overlap matrix unity?
  integer :: n_states                   !< Number of states to be computed

  !> libOMM variables
  logical :: new_overlap          !< Is a new overlap matrix provided?
  real*8  :: total_energy         !< Total energy value
  logical :: C_matrix_initialized !< Is there an initial guess of coefficients?
  integer :: omm_flavour = -1     !< "Flavours": 0 = Basic, 1 = Cholesky S, Cholesky U in S, 3 = preconditioning
  real*8  :: scale_kinetic        !< Scaling of the kinetic energy matrix
  logical :: calc_ed              !< Calculate the energy weigthed density matrix?
  real*8  :: eta                  !< Eigenspectrum shift parameter
  real*8  :: min_tol              !< Tolerance for minimization
  integer :: nk_times_nspin = -1  !< n_k_points * n_spin
  integer :: i_k_spin = -1        !< Combined k_point spin index
  logical :: omm_verbose          !< OMM output level
  logical :: do_dealloc           !< Deallocate internal storage?

  !> PEXSI Variables
  integer(c_intptr_t)    :: pexsi_plan !< PEXSI plan
  type(f_ppexsi_options) :: pexsi_options
  integer(c_int)         :: pexsi_info
  integer(c_int)         :: pexsi_output_file_index
  real(c_double)         :: mu_pexsi   !< Chemical potential
  real(c_double)         :: n_electrons_pexsi
  real(c_double)         :: mu_min_inertia
  real(c_double)         :: mu_max_inertia
  integer(c_int)         :: n_total_inertia_iter
  integer(c_int)         :: n_total_pexsi_iter
  real(c_double)         :: e_tot_h
  real(c_double)         :: e_tot_s
  real(c_double)         :: f_tot

  !> Method names
  enum, bind( C )
    enumerator :: ELPA, LIBOMM, PEXSI, CHESS
  end enum

contains

!>
!! Print message
!!
subroutine elsi_print(message)
   
   implicit none

   character(len=*), intent(in) :: message

   character(LEN=4096) :: string_message
   integer :: i_task

   do i_task = 0, n_procs - 1
      if(myid == i_task) then
         write(string_message, "(1X,'*** Proc',I5,': ',A)") &
               & myid, trim(message)

         write(*,'(A)') trim(string_message)
      endif
      call MPI_BARRIER(mpi_comm_global,mpierr)
   enddo

end subroutine 

!>
!! Clean shutdown.
!!
subroutine elsi_stop(message, caller)

   implicit none
   include "mpif.h"

   character(len=*), intent(in) :: message
   character(len=*), intent(in) :: caller

   character(LEN=4096) :: string_message
   integer :: i_task

   do i_task = 0, n_procs-1
      if(myid == i_task) then
         write(string_message, "(1X,'*** Proc',I5,' in ',A,': ',A)") &
               & myid, trim(caller), trim(message)
         write(*,'(A)') trim(string_message)
      endif
      call MPI_BARRIER(mpi_comm_global,mpierr)
   enddo

   if(n_procs > 1) then
      call MPI_Abort(mpi_comm_global, 0, mpierr)
   endif
     
   stop

end subroutine

!>
!! Print status of variables.
!!
subroutine elsi_variable_status()

   implicit none
   include "mpif.h"

   character(LEN=4096) :: string_message
   integer :: i_task

   if(myid == 0) then
      write(string_message, "(1X,'*** Proc',I5,&
            &' : Matrixsize global ',I5,' x ',I5)")&
            & myid, n_g_rank, n_g_rank
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'*** Proc',I5,&
            &' : Non_zero elements global ',I16)")&
            & myid, n_g_nonzero
      write(*,'(A)') trim(string_message)
        
      write(string_message, "(1X,'*** Proc',I5,&
            &' : Processgrid global ',I5,' x ',I5)")&
            & myid, n_p_rows, n_p_cols
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'*** Proc',I5,&
            &' : Number of States ',I5)")&
            & myid, n_states
      write(*,'(A)') trim(string_message)
        
      write(string_message, "(1X,'*** Proc',I5,&
            &' : Overlap is Unity ? ',L2)") &
            & myid, overlap_is_unity
      write(*,'(A)') trim(string_message)
   endif

   call MPI_BARRIER(mpi_comm_global, mpierr)

   do i_task = 0, n_procs
      if(i_task == myid) then
         write(string_message, "(1X,'*** Proc',I5,&
               &' : Matrixsize local ',I5,' x ',I5)") &
               & myid, n_l_rows, n_l_cols
         write(*,'(A)') trim(string_message)
         
         write(string_message, "(1X,'*** Proc',I5,&
               &' : Blocksize global ',I5,' x ',I5)") &
               & myid, n_b_rows, n_b_cols
         write(*,'(A)') trim(string_message)

         write(string_message, "(1X,'*** Proc',I5,&
               &' : Non_zero elements local ',I16)")&
               & myid, n_l_nonzero
         write(*,'(A)') trim(string_message)
      endif

      call MPI_BARRIER(mpi_comm_global, mpierr)
   enddo
     
end subroutine

!>
!! Set PEXSI variables to default.
!!
subroutine elsi_set_pexsi_default_options()

   implicit none
   include "mpif.h"

   ! Use the PEXSI Default options
   call f_ppexsi_set_default_options(pexsi_options)

   ! But do some modification to be more secure in convergence:

   ! Initial guess of lower bound for mu
   pexsi_options%muMin0   = -0.5d0
   ! Initial guess of upper bound for mu
   pexsi_options%muMax0   = 0.5d0
   ! Initial guess for mu (for the solver) (AG)
   pexsi_options%mu0      = 0.0d0
   ! An upper bound for the spectral radius of S−1H
   pexsi_options%deltaE   = 20d0
   ! Number of terms in the pole expansion
   pexsi_options%numPole  = 40
   ! Temperature, in the same unit as H
   pexsi_options%temperature = 0.0095 ! 3000K
   ! Safe guard criterion in terms of the chemical potential to
   ! reinvoke the inertia counting procedure.
   pexsi_options%muPEXSISafeGuard = 0.2d0
   ! Stopping criterion of the PEXSI iteration in terms of the number of
   ! electrons compared to numElectronExact
   pexsi_options%numElectronPEXSITolerance = 1d-3
   ! Output
   pexsi_options%verbosity = 1
   ! Maximum number of Iterations
   pexsi_options%maxPEXSIIter = 10

   call elsi_print_pexsi_options()

end subroutine

!>
!! Print status of PEXSI variables. 
!!
subroutine elsi_print_pexsi_options()

   implicit none
   include "mpif.h"

   character(LEN=4096) :: string_message

   if(myid == 0) then
      write(string_message, "(1X,'Temperature (H) ',F10.4)") &
            & pexsi_options%temperature
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'Spectral Gap (H) ',F10.4)") &
            & pexsi_options%gap
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'Delta E (H) ',F10.4)") &
            & pexsi_options%deltaE
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'Number of Poles ',I5)") &
            & pexsi_options%numPole
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'Use inertial count ',I2)") &
            & pexsi_options%isInertiaCount
      write(*,'(A)') trim(string_message)
      
      write(string_message, "(1X,'Max Pexsi Iterations ',I5)") &
            & pexsi_options%maxPEXSIIter
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'Minimal Mu (H) ',F10.4)") &
            & pexsi_options%muMin0
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'Maximal Mu (H) ',F10.4)") &
            & pexsi_options%muMax0
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'Start Mu (H) ',F10.4)") &
            & pexsi_options%mu0
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'Tolerance Mu (H) ',E10.1)") &
            & pexsi_options%muInertiaTolerance
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'Expansion Mu (H) ',F10.4)") &
            & pexsi_options%muInertiaExpansion
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'Safe guard Mu (H) ',F10.4)") &
            & pexsi_options%muPexsiSafeGuard
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'Tolerance Electrons (H) ',E10.1)") &
            & pexsi_options%numElectronPEXSITolerance
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'Sparcity ',F7.3)") &
            & (1d0 * n_g_nonzero / n_g_rank) / n_g_rank
      write(*,'(A)') trim(string_message)

   endif

end subroutine

!>
!! Set libOMM variables to default.
!!
subroutine elsi_set_omm_default_options()

   implicit none
   include "mpif.h"

   !< How to perform the calculation
   !! 0 = Basic
   !! 1 = Cholesky factorisation of S requested
   !! 2 = Cholesky already performed, U is provided in S
   !! 3 = Use preconditioning based on the energy density
   omm_flavour = 2
   !< Scaling of the kinetic energy matrix
   scale_kinetic = 5d0
   !< Calculate the energy weigthed density matrix?
   calc_ed = .false.
   !< Eigenspectrum shift parameter
   eta = 0.0d0
   !< Tolerance for minimization
   min_tol = 1.0d-9
   !< n_k_points * n_spin
   nk_times_nspin = 1
   !< Combined k_point spin index
   i_k_spin = 1
   !< OMM output level
   omm_verbose = .true.
   !< Deallocate internal storage?
   do_dealloc = .false.

   call elsi_print_omm_options()

end subroutine

!>
!! Print status of libOMM variables.
!!
subroutine elsi_print_omm_options()

   implicit none
   include "mpif.h"

   character(LEN=4096) :: string_message

   if(myid == 0) then
      write(string_message, "(1X,'  | Eta (H) ',F10.4)") eta
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'  | Kinetic scaling (H) ',F10.4)") scale_kinetic
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'  | Tolerance (H) ',E10.1)") min_tol
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'  | libOMM Flavour ',I1)") omm_flavour
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'  | Verbose ',L1)") omm_verbose
      write(*,'(A)') trim(string_message)

      write(string_message, "(1X,'  | Dealloc  ',L1)") do_dealloc
      write(*,'(A)') trim(string_message)
   endif

end subroutine

end module ELSI_DIMENSIONS
