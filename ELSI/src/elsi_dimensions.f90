!Copyright (c) 2015, ELSI consortium
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
!! This module contains variables accessible in ELSI and related modules
!!

module ELSI_DIMENSIONS

  use iso_c_binding
  use f_ppexsi_interface

  implicit none
  
  !>Method for EV Solver (ELPA=1,OMM=2,PEXSI=3)
  integer :: method = -1

  !>Mode for EV Solver (REAL_VALUES=1,COMPLEX_VALUES=2)
  integer :: mode = -1

  !> Global Matrix Dimensions
  integer :: n_g_rank       !< Rank of eigenvalue problem 
  
  !> Global Matrix Blocks
  integer :: n_b_rows       !< Number of rows 
  integer :: n_b_cols       !< Number of columns 

  !> Processor grid
  integer :: n_p_rows       !< Number of processor rows 
  integer :: n_p_cols       !< Number of processor columns
  
  !> Local Matrix Dimensions
  integer :: n_l_rows        !< Number of rows for this process
  integer :: n_l_cols        !< Number of columns for this process

  !> MPI variables
  integer :: myid            !< local process id
  integer :: n_procs         !< number of mpi processes
  integer :: mpi_comm_global !< global mpi communicator
  integer :: mpierr          !< mpi error handler
  logical :: external_mpi = .False. !< Has somebody else initialzed mpi
  logical :: mpi_is_setup = .False. !< Has somebody initialzed mpi

  !> BLACS variables
  integer :: blacs_ctxt     !< local blacs context
  integer :: sc_desc(9)     !< blacs descriptor
  integer :: mpi_comm_row   !< row    communicatior
  integer :: mpi_comm_col   !< column communicatior
  integer :: my_p_row       !< process row    position
  integer :: my_p_col       !< process column position
  integer :: blacs_info     !< Info code for blacs related calls 
  logical :: external_blacs = .False. !< Has somebody else initialized blacs? 
  logical :: blacs_is_setup = .False. !< Has somebody initialized blacs? 

  !> HDF5 variables
  integer :: h5err          !< HDF5 error code

  !> Sparse Matrix information
  integer :: n_g_nonzero    !< Number of nonzero elements of the global matrix
  integer :: n_l_nonzero    !< Number of nonzero elements of the local matrix

  !> Overlap
  logical :: overlap_is_unity = .True. !< Is the overlap unity
  integer :: n_eigenvectors            !< Number of eigenvectors to be calculated

  !> ELPA variables
  !< This parameter sets the threshold when to switch from ELPA2 to ELPA1
  real*8, parameter :: elpa_step_switch = 1.d0 
  
  !> OMM variables
  !< Is a new overlap used for OMM?
  logical :: new_overlap = .False.
  !< Total energy value
  real*8  :: total_energy
  !< Do we provide an initial guess for the Wannier functions?
  logical :: C_matrix_initialized = .False.
  !< How do we perform the calculation
  !! 0 = Basic
  !! 1 = Cholesky factorisation of S requested
  !! 2 = Cholesky already performed, U is provided in S
  !! 3 = Use preconditioning based on the energy density
  integer :: omm_flavour = -1
  !< scaling of the kinetic energy matrix, recommended round 5 Ha
  real*8 :: scale_kinetic = 5d0
  !< Calculate the energy weigthed density matrix
  logical :: calc_ed = .False.
  !< Eigenspectrum shift parameter
  real*8 :: eta = 0d0
  !< Tolerance for minimization
  real*8 :: min_tol = 1.0d-9
  !< n_k_points * n_spin
  integer :: nk_times_nspin = 1
  !< combined k_point spin index
  integer :: i_k_spin = 1
  !< Shall omm talk to us?
  logical :: omm_verbose = .True.
  !< Shall omm deallocate all internal?
  logical :: do_dealloc = .True.

  !> PEXSI Variables
  !< Pexsi Plan
  integer(c_intptr_t)   :: pexsi_plan
  !< Pexsi Options
  type(f_ppexsi_options):: pexsi_options
  !< Pexsi Info
  integer(c_int)        :: pexsi_info
  !< Pexsi output file
  integer(c_int)        :: pexsi_output_file_index
  !< Pexsi chemical potential
  real(c_double)        :: mu_pexsi
  !< Pexsi number of electrons
  real(c_double)        :: n_electrons_pexsi
  real(c_double)        :: mu_min_inertia
  real(c_double)        :: mu_max_inertia
  integer(c_int)        :: n_total_inertia_iter
  integer(c_int)        :: n_total_pexsi_iter
  real(c_double)        :: e_tot_h
  real(c_double)        :: e_tot_s
  real(c_double)        :: f_tot


  !> Physics Variables
  !< The number of electrons
  real*8  :: n_electrons

  !> Method names
  enum, bind( C )
    enumerator :: ELPA, OMM_DENSE, PEXSI
  end enum

  enum, bind( C )
    enumerator :: REAL_VALUES, COMPLEX_VALUES
  end enum   


  contains

!>
!! Elsi print of process
!!
  subroutine elsi_print(message)
   
      implicit none

      character(len=*), intent(in) :: message

      character(LEN=4096) :: string_message

      integer :: i_task

       do i_task = 0, n_procs - 1
        if (myid == i_task) then
          write(string_message, "(1X,'*** Proc',I5,': ',A)") &
            & myid, trim(message)

          write(*,'(A)') trim(string_message)
        end if
        call MPI_BARRIER(mpi_comm_global,mpierr)
      end do

  end subroutine 

!>
!! Clean shutdown in case of error
!!
subroutine elsi_stop(message, caller)

      implicit none
      include "mpif.h"

      character(len=*), intent(in) :: message
      character(len=*), intent(in) :: caller

      character(LEN=4096) :: string_message
      integer :: i_task

      do i_task = 0, n_procs - 1
        if (myid == i_task) then
          write(string_message, "(1X,'*** Proc',I5,' in ',A,': ',A)") &
             & myid, trim(caller), trim(message)
          write(*,'(A)') trim(string_message)
        end if
        call MPI_BARRIER(mpi_comm_global,mpierr)
      end do

      if (n_procs > 1) then
         call MPI_Abort(mpi_comm_global, 0, mpierr)
      end if
        
      stop

end subroutine elsi_stop

!>
!! Give Status of ELSI variables 
!!
subroutine elsi_variable_status()

      implicit none
      include "mpif.h"

      character(LEN=4096) :: string_message
      integer :: i_task

      do i_task = 0, n_procs
         if (i_task == myid) then
            write(string_message, "(1X,'*** Proc',I5,&
            ' : Matrixsize global ',I5,' x ',I5)") &
              & myid, n_g_rank, n_g_rank
            write(*,'(A)') trim(string_message)
            write(string_message, "(1X,'*** Proc',I5,&
            ' : Matrixsize local ',I5,' x ',I5)") &
              & myid, n_l_rows, n_l_cols
            write(*,'(A)') trim(string_message)
            write(string_message, "(1X,'*** Proc',I5,&
            ' : Blocksize global ',I5,' x ',I5)") &
              & myid, n_b_rows, n_b_cols
            write(*,'(A)') trim(string_message)
            write(string_message, "(1X,'*** Proc',I5,&
            ' : Processgrid global ',I5,' x ',I5)") &
              & myid, n_p_rows, n_p_cols
            write(*,'(A)') trim(string_message)
         end if

         call MPI_BARRIER(mpi_comm_global,mpierr)
      end do
        
end subroutine elsi_variable_status


end module ELSI_DIMENSIONS
