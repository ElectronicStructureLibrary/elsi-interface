! Copyright (c) 2015-2017, the ELSI team. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!  * Neither the name of the "ELectronic Structure Infrastructure" project nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
! OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!>
!! This module provides interfaces to PEXSI.
!!
module ELSI_PEXSI

   use ISO_C_BINDING
   use ELSI_CONSTANTS, only: BLACS_DENSE,UNSET
   use ELSI_DIMENSIONS, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4
   use ELSI_TIMERS
   use ELSI_UTILS
   use F_PPEXSI_INTERFACE

   implicit none
   private

   public :: elsi_init_pexsi
   public :: elsi_solve_evp_pexsi
   public :: elsi_set_pexsi_default
   public :: elsi_print_pexsi_options

contains

!>
!! This routine initializes PEXSI and its processor grid.
!!
subroutine elsi_init_pexsi(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4) :: n_rows_tmp
   integer(kind=i4) :: n_groups
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_init_pexsi"

   ! Safety
   if(elsi_h%pexsi_driver == 1) then
      elsi_h%n_mu_points = 1
   endif

   if(elsi_h%n_elsi_calls == 1) then
      n_groups = elsi_h%pexsi_options%numPole*elsi_h%n_mu_points

      if(elsi_h%n_p_per_pole_pexsi == UNSET) then
         if(mod(elsi_h%n_procs,n_groups) == 0) then
            elsi_h%n_p_per_pole_pexsi = elsi_h%n_procs/n_groups

            call elsi_statement_print("  PEXSI parallel over poles.",elsi_h)
            write(info_str,"(A,I13)") "  | Number of MPI tasks per pole: ",&
               elsi_h%n_p_per_pole_pexsi
            call elsi_statement_print(info_str,elsi_h)
         else
            call elsi_stop("  PEXSI not parallel over poles. High"//&
                    " performance of PEXSI is expected if the number"//&
                    " of MPI tasks is a multiple of the number of"//&
                    " PEXSI poles. Please adjust either the number"//&
                    " of MPI tasks, or the number of poles."//&
                    " Exiting...",elsi_h,caller)
         endif
      endif

      ! Set square-like process grid for selected inversion of each pole
      do n_rows_tmp = nint(sqrt(real(elsi_h%n_p_per_pole_pexsi))),2,-1
         if(mod(elsi_h%n_p_per_pole_pexsi,n_rows_tmp) == 0) exit
      enddo

      elsi_h%n_p_rows_pexsi = n_rows_tmp
      elsi_h%n_p_cols_pexsi = elsi_h%n_p_per_pole_pexsi/elsi_h%n_p_rows_pexsi

      ! PEXSI process grid
      elsi_h%my_p_col_pexsi = mod(elsi_h%myid,elsi_h%n_p_per_pole_pexsi)
      elsi_h%my_p_row_pexsi = elsi_h%myid/elsi_h%n_p_per_pole_pexsi

      ! PEXSI uses a pure block distribution in the first process row
      elsi_h%n_b_rows_pexsi = elsi_h%n_g_size

      ! The last process holds all remaining columns
      elsi_h%n_b_cols_pexsi = elsi_h%n_g_size/elsi_h%n_p_per_pole_pexsi
      if(elsi_h%my_p_col_pexsi == elsi_h%n_p_per_pole_pexsi-1) then
         elsi_h%n_b_cols_pexsi = elsi_h%n_g_size-&
                                    (elsi_h%n_p_per_pole_pexsi-1)*elsi_h%n_b_cols_pexsi
      endif

      elsi_h%n_l_rows_pexsi = elsi_h%n_b_rows_pexsi
      elsi_h%n_l_cols_pexsi = elsi_h%n_b_cols_pexsi

      ! Only master process outputs
      if(elsi_h%myid == 0) then
         elsi_h%pexsi_output_file_index = 0
      else
         elsi_h%pexsi_output_file_index = -1
      endif

      elsi_h%pexsi_plan = f_ppexsi_plan_initialize(elsi_h%mpi_comm,&
                             elsi_h%n_p_rows_pexsi,elsi_h%n_p_cols_pexsi,&
                             elsi_h%pexsi_output_file_index,elsi_h%pexsi_info)

      if(elsi_h%pexsi_info /= 0) then
         call elsi_stop(" PEXSI plan initialization failed. Exiting...",elsi_h,caller)
      endif

      elsi_h%pexsi_started = .true.
   endif

end subroutine

!>
!! This routine interfaces to PEXSI.
!!
subroutine elsi_solve_evp_pexsi(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   real(kind=r8), save :: this_pexsi_tol = 1.0e-2_r8
   integer(kind=i4)    :: mpierr
   character*200       :: info_str

   character*40, parameter :: caller = "elsi_solve_evp_pexsi"

   call elsi_start_density_matrix_time(elsi_h)

   if(elsi_h%pexsi_driver == 1) then
      if(elsi_h%small_pexsi_tol) then
         elsi_h%pexsi_options%numElectronPEXSITolerance = this_pexsi_tol

         write(info_str,"(A,E10.1)") "  | Current tolerance of number of electrons: ",&
            this_pexsi_tol
         call elsi_statement_print(info_str,elsi_h)
      endif
   endif

   if(elsi_h%n_elsi_calls == 1) then
      elsi_h%pexsi_options%isSymbolicFactorize = 1
   else
      elsi_h%pexsi_options%isSymbolicFactorize = 0
   endif

   if(.not. allocated(elsi_h%e_den_mat_pexsi)) then
      call elsi_allocate(elsi_h,elsi_h%e_den_mat_pexsi,elsi_h%nnz_l_pexsi,"e_den_mat_pexsi",caller)
   endif
   elsi_h%e_den_mat_pexsi = 0.0_r8

   if(.not. allocated(elsi_h%f_den_mat_pexsi)) then
      call elsi_allocate(elsi_h,elsi_h%f_den_mat_pexsi,elsi_h%nnz_l_pexsi,"f_den_mat_pexsi",caller)
   endif
   elsi_h%f_den_mat_pexsi = 0.0_r8

   ! Load sparse matrices for PEXSI
   if(elsi_h%overlap_is_unit) then
      call f_ppexsi_load_real_hs_matrix(elsi_h%pexsi_plan,elsi_h%pexsi_options,elsi_h%n_g_size,&
              elsi_h%nnz_g,elsi_h%nnz_l_pexsi,elsi_h%n_l_cols_pexsi,elsi_h%col_ptr_ccs,&
              elsi_h%row_ind_ccs,elsi_h%ham_real_ccs,1,elsi_h%ovlp_real_ccs,elsi_h%pexsi_info)
   else
      call f_ppexsi_load_real_hs_matrix(elsi_h%pexsi_plan,elsi_h%pexsi_options,elsi_h%n_g_size,&
              elsi_h%nnz_g,elsi_h%nnz_l_pexsi,elsi_h%n_l_cols_pexsi,elsi_h%col_ptr_ccs,&
              elsi_h%row_ind_ccs,elsi_h%ham_real_ccs,0,elsi_h%ovlp_real_ccs,elsi_h%pexsi_info)
   endif

   if(elsi_h%pexsi_info /= 0) then
      call elsi_stop(" PEXSI not able to load H/S matrix. Exiting...",elsi_h,caller)
   endif

   if(elsi_h%pexsi_options%isInertiaCount == 0) then
      call elsi_statement_print("  PEXSI inertia counting skipped",elsi_h)
   endif

   ! Solve the eigenvalue problem
   call elsi_statement_print("  Starting PEXSI density matrix solver",elsi_h)

   if(elsi_h%pexsi_driver == 1) then
      call f_ppexsi_dft_driver(elsi_h%pexsi_plan,elsi_h%pexsi_options,elsi_h%n_electrons,&
              elsi_h%mu,elsi_h%n_electrons_pexsi,elsi_h%mu_min_inertia,elsi_h%mu_max_inertia,&
              elsi_h%n_total_inertia_iter,elsi_h%n_total_pexsi_iter,elsi_h%pexsi_info)
   else
      call f_ppexsi_dft_driver3(elsi_h%pexsi_plan,elsi_h%pexsi_options,elsi_h%n_electrons,&
              2,elsi_h%n_mu_points,elsi_h%mu,elsi_h%n_electrons_pexsi,&
              elsi_h%n_total_inertia_iter,elsi_h%pexsi_info)
   endif

   if(elsi_h%pexsi_info /= 0) then
      call elsi_stop(" PEXSI DFT driver not able to solve problem. Exiting...",elsi_h,caller)
   endif

   ! Reuse chemical potential
   if(elsi_h%pexsi_driver == 1) then
      if(abs(elsi_h%mu-elsi_h%pexsi_options%mu0) > 1.0e-3_r8) then
         elsi_h%pexsi_options%isInertiaCount = 1
      else
         elsi_h%pexsi_options%isInertiaCount = 0
      endif

      elsi_h%pexsi_options%mu0 = elsi_h%mu

      if(elsi_h%small_pexsi_tol) then
         if(abs(elsi_h%n_electrons-elsi_h%n_electrons_pexsi) < this_pexsi_tol) then
            if(1.0e-1_r8*this_pexsi_tol > elsi_h%final_pexsi_tol) then
               this_pexsi_tol = 1.0e-1_r8*this_pexsi_tol
            else
               this_pexsi_tol = elsi_h%final_pexsi_tol
            endif
         endif
      endif
   else
      elsi_h%pexsi_options%muMin0 = elsi_h%pexsi_options%muMin0-1.0_r8
      elsi_h%pexsi_options%muMax0 = elsi_h%pexsi_options%muMax0+1.0_r8
   endif

   ! Get the results
   if((elsi_h%my_p_row_pexsi == 0) .or. (elsi_h%matrix_storage_format == BLACS_DENSE)) then
      if(elsi_h%pexsi_driver == 1) then
         call f_ppexsi_retrieve_real_dft_matrix(elsi_h%pexsi_plan,elsi_h%den_mat_ccs,&
                 elsi_h%e_den_mat_pexsi,elsi_h%f_den_mat_pexsi,elsi_h%energy_hdm,&
                 elsi_h%energy_sedm,elsi_h%free_energy,elsi_h%pexsi_info)
      else
         call f_ppexsi_retrieve_real_dft_matrix2(elsi_h%pexsi_plan,elsi_h%den_mat_ccs,&
                 elsi_h%e_den_mat_pexsi,elsi_h%f_den_mat_pexsi,elsi_h%energy_hdm,&
                 elsi_h%energy_sedm,elsi_h%free_energy,elsi_h%pexsi_info)
      endif
   endif

   if(elsi_h%pexsi_info /= 0) then
      call elsi_stop(" PEXSI not able to retrieve solution. Exiting...",elsi_h,caller)
   endif

   call MPI_Barrier(elsi_h%mpi_comm,mpierr)
   call elsi_stop_density_matrix_time(elsi_h)

end subroutine

!>
!! This routine sets default PEXSI parameters.
!!
subroutine elsi_set_pexsi_default(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   character*40, parameter :: caller = "elsi_set_pexsi_default"

   ! Use the PEXSI Default options
   call f_ppexsi_set_default_options(elsi_h%pexsi_options)

   ! Use 1 process in symbolic factorization
   elsi_h%pexsi_options%npSymbFact = 1

   ! PEXSI DFT driver
   elsi_h%pexsi_driver = 2

   ! Number of poles
   elsi_h%pexsi_options%numPole = 20

   ! Number of mu points if using Moussa's pole expansion
   elsi_h%n_mu_points = 2

   ! Output level
   elsi_h%pexsi_options%verbosity = 0

end subroutine

!>
!! This routine prints PEXSI settings.
!!
subroutine elsi_print_pexsi_options(elsi_h)

   implicit none

   type(elsi_handle), intent(in) :: elsi_h !< Handle

   character*200 :: info_str

   character*40, parameter :: caller = "elsi_print_pexsi_options"

   write(info_str,"(A)") "  PEXSI settings (in the same unit of Hamiltonian):"
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | PEXSI DFT driver ',I5)") &
      elsi_h%pexsi_driver
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Temperature ',F10.4)") &
      elsi_h%pexsi_options%temperature
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Spectral gap ',F10.4)") &
      elsi_h%pexsi_options%gap
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Number of poles ',I5)") &
      elsi_h%pexsi_options%numPole
   call elsi_statement_print(info_str,elsi_h)

   if(elsi_h%pexsi_driver == 2) then
      write(info_str,"(1X,' | Number of mu points ',I5)") &
         elsi_h%n_mu_points
      call elsi_statement_print(info_str,elsi_h)
   endif

   if(elsi_h%pexsi_driver == 1) then
      write(info_str,"(1X,' | Max PEXSI iterations ',I5)") &
         elsi_h%pexsi_options%maxPEXSIIter
      call elsi_statement_print(info_str,elsi_h)

      write(info_str,"(1X,' | Initial guess of chemical potential ',F10.4)") &
         elsi_h%pexsi_options%mu0
      call elsi_statement_print(info_str,elsi_h)

      write(info_str,"(1X,' | Safeguard of chemical potential ',F10.4)") &
         elsi_h%pexsi_options%muPexsiSafeGuard
      call elsi_statement_print(info_str,elsi_h)

      write(info_str,"(1X,' | Tolerance of number of electrons ',E10.1)") &
         elsi_h%pexsi_options%numElectronPEXSITolerance
      call elsi_statement_print(info_str,elsi_h)
   endif

   write(info_str,"(1X,' | Lower bound of chemical potential ',F10.4)") &
      elsi_h%pexsi_options%muMin0
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Upper bound of chemical potential ',F10.4)") &
      elsi_h%pexsi_options%muMax0
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Tolerance of chemical potential ',E10.1)") &
      elsi_h%pexsi_options%muInertiaTolerance
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Number of processes for symbolic factorization ',I5)") &
      elsi_h%pexsi_options%npSymbFact
   call elsi_statement_print(info_str,elsi_h)

end subroutine

end module ELSI_PEXSI
