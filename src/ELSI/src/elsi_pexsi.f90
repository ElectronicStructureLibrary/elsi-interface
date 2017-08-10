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

   use ELSI_CONSTANTS, only: BLACS_DENSE,REAL_VALUES,COMPLEX_VALUES,UNSET
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS
   use F_PPEXSI_INTERFACE

   implicit none
   private

   public :: elsi_init_pexsi
   public :: elsi_solve_evp_pexsi
   public :: elsi_compute_edm_pexsi
   public :: elsi_set_pexsi_default
   public :: elsi_print_pexsi_options

contains

!>
!! This routine initializes PEXSI and its processor grid.
!!
subroutine elsi_init_pexsi(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4) :: n_rows_tmp
   integer(kind=i4) :: n_groups
   integer(kind=i4) :: output_id
   integer(kind=i4) :: ierr
   integer(kind=i4) :: mpierr
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_init_pexsi"

   if(elsi_h%n_elsi_calls == 1) then
      if(elsi_h%n_p_per_pole == UNSET) then
         elsi_h%n_p_per_pole = elsi_h%n_procs/(elsi_h%pexsi_options%numPole*&
                                  elsi_h%pexsi_options%nPoints)
      endif

      write(info_str,"(1X,' | Number of MPI tasks per pole                   ',I10)") &
         elsi_h%n_p_per_pole
      call elsi_statement_print(info_str,elsi_h)

      ! Set square-like process grid for selected inversion of each pole
      do n_rows_tmp = nint(sqrt(real(elsi_h%n_p_per_pole))),2,-1
         if(mod(elsi_h%n_p_per_pole,n_rows_tmp) == 0) exit
      enddo

      elsi_h%n_p_rows_pexsi = n_rows_tmp
      elsi_h%n_p_cols_pexsi = elsi_h%n_p_per_pole/n_rows_tmp

      ! PEXSI process grid
      elsi_h%my_p_col_pexsi = mod(elsi_h%myid,elsi_h%n_p_per_pole)
      elsi_h%my_p_row_pexsi = elsi_h%myid/elsi_h%n_p_per_pole

      ! Point parallelization
      elsi_h%n_p_per_point = elsi_h%n_procs/elsi_h%pexsi_options%nPoints
      elsi_h%my_point      = elsi_h%myid/elsi_h%n_p_per_point
      elsi_h%myid_point    = mod(elsi_h%myid,elsi_h%n_p_per_point)

      ! PEXSI MPI communicators
      call MPI_Comm_split(elsi_h%mpi_comm,elsi_h%my_p_col_pexsi,&
              elsi_h%my_p_row_pexsi,elsi_h%comm_among_pole,mpierr)

      call MPI_Comm_split(elsi_h%mpi_comm,elsi_h%my_p_row_pexsi,&
              elsi_h%my_p_col_pexsi,elsi_h%comm_in_pole,mpierr)

      call MPI_Comm_split(elsi_h%mpi_comm,elsi_h%myid_point,&
              elsi_h%my_point,elsi_h%comm_among_point,mpierr)

      call MPI_Comm_split(elsi_h%mpi_comm,elsi_h%my_point,&
              elsi_h%myid_point,elsi_h%comm_in_point,mpierr)

      ! 1D block distribution
      elsi_h%n_l_cols_sp = elsi_h%n_basis/elsi_h%n_p_per_pole

      ! The last process holds all remaining columns
      if(elsi_h%my_p_col_pexsi == elsi_h%n_p_per_pole-1) then
         elsi_h%n_l_cols_sp = elsi_h%n_basis-&
                                 (elsi_h%n_p_per_pole-1)*elsi_h%n_l_cols_sp
      endif

      ! Only one process outputs
      if(elsi_h%myid_all == 0) then
         output_id = 0
      else
         output_id = -1
      endif

      elsi_h%pexsi_plan = f_ppexsi_plan_initialize(elsi_h%mpi_comm,&
                             elsi_h%n_p_rows_pexsi,elsi_h%n_p_cols_pexsi,&
                             output_id,ierr)

      if(ierr /= 0) then
         call elsi_stop(" PEXSI initialization failed. Exiting...",&
                  elsi_h,caller)
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

   real(kind=r8)    :: ne_drv
   real(kind=r8)    :: mu_range
   real(kind=r8)    :: shift_width
   real(kind=r8)    :: local_energy
   real(kind=r8)    :: factor_min
   real(kind=r8)    :: factor_max
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   integer(kind=i4) :: n_iner_steps
   integer(kind=i4) :: n_shift
   integer(kind=i4) :: aux_min
   integer(kind=i4) :: aux_max
   integer(kind=i4) :: i
   integer(kind=i4) :: idx
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: ierr
   logical          :: converged
   character*200    :: info_str

   real(kind=r8),    allocatable :: shifts(:)
   real(kind=r8),    allocatable :: inertias(:)
   real(kind=r8),    allocatable :: ne_lower(:)
   real(kind=r8),    allocatable :: ne_upper(:)
   real(kind=r8),    allocatable :: tmp_real(:)
   complex(kind=r8), allocatable :: tmp_complex(:)
   real(kind=r8),    allocatable :: send_buffer(:)
   complex(kind=r8), allocatable :: send_buffer_complex(:)

   real(kind=r8), external :: ddot

   character*40,  parameter :: caller = "elsi_solve_evp_pexsi"

   ! Load sparse matrices for PEXSI
   select case(elsi_h%matrix_data_type)
   case(REAL_VALUES)
      if(elsi_h%ovlp_is_unit) then
         call f_ppexsi_load_real_hs_matrix(elsi_h%pexsi_plan,&
                 elsi_h%pexsi_options,elsi_h%n_basis,elsi_h%nnz_g,&
                 elsi_h%nnz_l_sp,elsi_h%n_l_cols_sp,elsi_h%col_ptr_ccs,&
                 elsi_h%row_ind_ccs,elsi_h%ham_real_ccs,1,&
                 elsi_h%ovlp_real_ccs,ierr)
      else
         call f_ppexsi_load_real_hs_matrix(elsi_h%pexsi_plan,&
                 elsi_h%pexsi_options,elsi_h%n_basis,elsi_h%nnz_g,&
                 elsi_h%nnz_l_sp,elsi_h%n_l_cols_sp,elsi_h%col_ptr_ccs,&
                 elsi_h%row_ind_ccs,elsi_h%ham_real_ccs,0,&
                 elsi_h%ovlp_real_ccs,ierr)
      endif
   case(COMPLEX_VALUES)
      if(elsi_h%ovlp_is_unit) then
         call f_ppexsi_load_complex_hs_matrix(elsi_h%pexsi_plan,&
                 elsi_h%pexsi_options,elsi_h%n_basis,elsi_h%nnz_g,&
                 elsi_h%nnz_l_sp,elsi_h%n_l_cols_sp,elsi_h%col_ptr_ccs,&
                 elsi_h%row_ind_ccs,elsi_h%ham_complex_ccs,1,&
                 elsi_h%ovlp_complex_ccs,ierr)
      else
         call f_ppexsi_load_complex_hs_matrix(elsi_h%pexsi_plan,&
                 elsi_h%pexsi_options,elsi_h%n_basis,elsi_h%nnz_g,&
                 elsi_h%nnz_l_sp,elsi_h%n_l_cols_sp,elsi_h%col_ptr_ccs,&
                 elsi_h%row_ind_ccs,elsi_h%ham_complex_ccs,0,&
                 elsi_h%ovlp_complex_ccs,ierr)
      endif
   end select

   if(ierr /= 0) then
      call elsi_stop(" PEXSI load matrices failed. Exiting...",elsi_h,caller)
   endif

   call elsi_statement_print("  Starting PEXSI density matrix solver",elsi_h)

   ! Symbolic factorization
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_get_time(elsi_h,t0)

      select case(elsi_h%matrix_data_type)
      case(REAL_VALUES)
         call f_ppexsi_symbolic_factorize_real_symmetric_matrix(&
                 elsi_h%pexsi_plan,elsi_h%pexsi_options,ierr)

         call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(&
                 elsi_h%pexsi_plan,elsi_h%pexsi_options,ierr)
      case(COMPLEX_VALUES)
         call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(&
                 elsi_h%pexsi_plan,elsi_h%pexsi_options,ierr)

         call f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix(&
                 elsi_h%pexsi_plan,elsi_h%pexsi_options,ierr)
      end select

      call elsi_get_time(elsi_h,t1)

      write(info_str,"('  Finished symbolic factorization')")
      call elsi_statement_print(info_str,elsi_h)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_statement_print(info_str,elsi_h)
   endif

   if(ierr /= 0) then
      call elsi_stop(" Symbolic factorization failed. Exiting...",elsi_h,caller)
   endif

   ! Inertia counting
   call elsi_get_time(elsi_h,t0)

   n_iner_steps = 0
   mu_range = elsi_h%pexsi_options%muMax0-elsi_h%pexsi_options%muMin0
   n_shift = max(10,elsi_h%n_procs/elsi_h%n_p_per_pole)

   call elsi_allocate(elsi_h,shifts,n_shift,"shifts",caller)
   call elsi_allocate(elsi_h,inertias,n_shift,"inertias",caller)
   call elsi_allocate(elsi_h,ne_lower,n_shift,"ne_lower",caller)
   call elsi_allocate(elsi_h,ne_upper,n_shift,"ne_upper",caller)

   if(.not. elsi_h%spin_is_set) then
      if(elsi_h%n_spins == 2) then
         elsi_h%spin_degen = 1.0_r8
      else
         elsi_h%spin_degen = 2.0_r8
      endif
   endif

   do while((n_iner_steps < 10) .and. &
            (mu_range > elsi_h%pexsi_options%muInertiaTolerance))
      n_iner_steps = n_iner_steps+1

      shift_width = mu_range/(n_shift-1)

      ne_lower = 0.0_r8
      ne_upper = elsi_h%n_basis*elsi_h%spin_degen

      do i = 1,n_shift
         shifts(i)   = elsi_h%pexsi_options%muMin0+(i-1)*shift_width
      enddo

      select case(elsi_h%matrix_data_type)
      case(REAL_VALUES)
         call f_ppexsi_inertia_count_real_matrix(elsi_h%pexsi_plan,&
                 elsi_h%pexsi_options,n_shift,shifts,inertias,ierr)
      case(COMPLEX_VALUES)
         call f_ppexsi_inertia_count_complex_matrix(elsi_h%pexsi_plan,&
                 elsi_h%pexsi_options,n_shift,shifts,inertias,ierr)
      end select

      inertias = inertias*elsi_h%spin_degen*elsi_h%i_weight

      ! Get global inertias
      if(elsi_h%n_spins*elsi_h%n_kpts > 1) then
         call elsi_allocate(elsi_h,send_buffer,n_shift,"send_buffer",caller)

         if(elsi_h%myid == 0) then
            send_buffer = inertias
         else
            send_buffer = 0.0_r8
         endif

         call MPI_Allreduce(send_buffer,inertias,n_shift,mpi_real8,&
                 mpi_sum,elsi_h%mpi_comm_all,mpierr)

         call elsi_deallocate(elsi_h,send_buffer,"send_buffer")
      endif

      idx = ceiling(3*elsi_h%pexsi_options%temperature/shift_width)

      do i = idx+1,n_shift
         ne_lower(i) = 0.5_r8*(inertias(i-idx)+inertias(i))
         ne_upper(i-idx) = ne_lower(i)
      enddo

      aux_min = 1
      aux_max = n_shift

      do i = 2,n_shift-1
         if((ne_upper(i) < elsi_h%n_electrons) .and. &
            (ne_upper(i+1) .ge. elsi_h%n_electrons))  then
            aux_min = i
         endif

         if((ne_lower(i) > elsi_h%n_electrons) .and. &
            (ne_lower(i-1) .le. elsi_h%n_electrons)) then
            aux_max = i
         endif
      enddo

      if((aux_min == 1) .and. (aux_max == n_shift)) then
         exit
      else
         elsi_h%pexsi_options%muMin0 = shifts(aux_min)
         elsi_h%pexsi_options%muMax0 = shifts(aux_max)
         mu_range = elsi_h%pexsi_options%muMax0-elsi_h%pexsi_options%muMin0
      endif
   enddo

   call elsi_deallocate(elsi_h,shifts,"shifts")
   call elsi_deallocate(elsi_h,inertias,"inertias")
   call elsi_deallocate(elsi_h,ne_lower,"ne_lower")
   call elsi_deallocate(elsi_h,ne_upper,"ne_upper")

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished inertia counting')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

   if(ierr /= 0) then
      call elsi_stop(" Inertia counting failed. Exiting...",elsi_h,caller)
   endif

   ! Fermi operator expansion
   call elsi_get_time(elsi_h,t0)

   shift_width = mu_range/(elsi_h%pexsi_options%nPoints+1)

   call elsi_allocate(elsi_h,shifts,elsi_h%pexsi_options%nPoints,&
           "shifts",caller)

   do i = 1,elsi_h%pexsi_options%nPoints
      shifts(i) = elsi_h%pexsi_options%muMin0+i*shift_width
   enddo

   do i = 1,elsi_h%pexsi_options%nPoints
      elsi_h%mu = shifts(i)

      if(elsi_h%my_point == i-1) then
         select case(elsi_h%matrix_data_type)
         case(REAL_VALUES)
            call f_ppexsi_calculate_fermi_operator_real3(elsi_h%pexsi_plan,&
                    elsi_h%pexsi_options,elsi_h%mu,elsi_h%n_electrons,&
                    elsi_h%ne_pexsi,ne_drv,ierr)
         case(COMPLEX_VALUES)
            call f_ppexsi_calculate_fermi_operator_complex(elsi_h%pexsi_plan,&
                    elsi_h%pexsi_options,elsi_h%mu,elsi_h%n_electrons,&
                    elsi_h%ne_pexsi,ne_drv,ierr)
         end select
      endif
   enddo

   call elsi_allocate(elsi_h,send_buffer,elsi_h%pexsi_options%nPoints,&
           "send_buffer",caller)
   if(elsi_h%n_elsi_calls == 1) then
      call elsi_allocate(elsi_h,elsi_h%ne_vec,elsi_h%pexsi_options%nPoints,&
              "ne_vec",caller)
   endif

   send_buffer(elsi_h%my_point+1) = elsi_h%ne_pexsi*elsi_h%i_weight

   call MPI_Allreduce(send_buffer,elsi_h%ne_vec,elsi_h%pexsi_options%nPoints,&
           mpi_real8,mpi_sum,elsi_h%comm_among_point,mpierr)

   ! Get global number of electrons
   if(elsi_h%n_spins*elsi_h%n_kpts > 1) then
      if(elsi_h%myid == 0) then
         send_buffer = elsi_h%ne_vec
      else
         send_buffer = 0.0_r8
      endif

      call MPI_Allreduce(send_buffer,elsi_h%ne_vec,&
              elsi_h%pexsi_options%nPoints,mpi_real8,mpi_sum,&
              elsi_h%mpi_comm_all,mpierr)
   endif

   call elsi_deallocate(elsi_h,send_buffer,"send_buffer")

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished Fermi operator calculation')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

   if(ierr /= 0) then
      call elsi_stop(" Fermi operator calculation failed. Exiting...",&
              elsi_h,caller)
   endif

   ! Get density matrix
   call elsi_get_time(elsi_h,t0)

   select case(elsi_h%matrix_data_type)
   case(REAL_VALUES)
      call elsi_allocate(elsi_h,tmp_real,elsi_h%nnz_l_sp,&
              "tmp_real",caller)

      call f_ppexsi_retrieve_real_dm(elsi_h%pexsi_plan,tmp_real,&
              local_energy,ierr)
   case(COMPLEX_VALUES)
      call elsi_allocate(elsi_h,tmp_complex,elsi_h%nnz_l_sp,&
              "tmp_complex",caller)

      call f_ppexsi_retrieve_complex_dm(elsi_h%pexsi_plan,tmp_complex,&
              local_energy,ierr)
   end select

   if(ierr /= 0) then
      call elsi_stop(" Retrieving density matirx failed. Exiting...",&
              elsi_h,caller)
   endif

   ! Check convergence
   converged = .false.
   aux_min = 0
   aux_max = n_shift+1

   do i = 1,elsi_h%pexsi_options%nPoints
      if(elsi_h%ne_vec(i) < elsi_h%n_electrons-&
         elsi_h%pexsi_options%numElectronPEXSITolerance) then
         elsi_h%pexsi_options%muMin0 = shifts(i)
         aux_min = i
      endif
   enddo

   do i = elsi_h%pexsi_options%nPoints,1,-1
      if(elsi_h%ne_vec(i) > elsi_h%n_electrons+&
         elsi_h%pexsi_options%numElectronPEXSITolerance) then
         elsi_h%pexsi_options%muMax0 = shifts(i)
         aux_max = i
      endif
   enddo

   if(elsi_h%pexsi_options%nPoints == 1) then
      ! Scale density matrix
      select case(elsi_h%matrix_data_type)
      case(REAL_VALUES)
         tmp_real = (elsi_h%n_electrons/elsi_h%ne_pexsi)*tmp_real
      case(COMPLEX_VALUES)
         tmp_complex = (elsi_h%n_electrons/elsi_h%ne_pexsi)*tmp_complex
      end select

      converged = .true.
      elsi_h%mu = shifts(1)
   else
      ! Safety check
      if(aux_min == 0) then
         aux_min = 1

         if(aux_max .le. aux_min) then
            aux_max = 2
         endif
      endif

      if(aux_max == n_shift+1) then
         aux_max = n_shift

         if(aux_min .ge. aux_max) then
            aux_min = n_shift-1
         endif
      endif

      do i = aux_min,aux_max
         if(abs(elsi_h%ne_vec(i)-elsi_h%n_electrons) < &
            elsi_h%pexsi_options%numElectronPEXSITolerance) then
            elsi_h%mu = shifts(i)
            converged = .true.

            select case(elsi_h%matrix_data_type)
            case(REAL_VALUES)
               call MPI_Bcast(tmp_real,elsi_h%nnz_l_sp,mpi_real8,i,&
                       elsi_h%comm_among_point,mpierr)
            case(COMPLEX_VALUES)
               call MPI_Bcast(tmp_complex,elsi_h%nnz_l_sp,mpi_complex16,&
                       i,elsi_h%comm_among_point,mpierr)
            end select

            exit
         endif
      enddo
   endif

   ! Interpolation
   if(.not. converged) then
      ! Chemical potential
      elsi_h%mu = shifts(aux_min)+(elsi_h%n_electrons-elsi_h%ne_vec(aux_min))/&
                     (elsi_h%ne_vec(aux_max)-elsi_h%ne_vec(aux_min))*&
                     (shifts(aux_max)-shifts(aux_min))

      ! Density matrix
      factor_min = (elsi_h%ne_vec(aux_max)-elsi_h%n_electrons)/&
                      (elsi_h%ne_vec(aux_max)-elsi_h%ne_vec(aux_min))
      factor_max = (elsi_h%n_electrons-elsi_h%ne_vec(aux_min))/&
                      (elsi_h%ne_vec(aux_max)-elsi_h%ne_vec(aux_min))

      select case(elsi_h%matrix_data_type)
      case(REAL_VALUES)
         call elsi_allocate(elsi_h,send_buffer,elsi_h%nnz_l_sp,&
                 "send_buffer",caller)

         if(elsi_h%my_point == aux_min-1) then
            send_buffer = factor_min*tmp_real
         elseif(elsi_h%my_point == aux_max-1) then
            send_buffer = factor_max*tmp_real
         endif

         call MPI_Allreduce(send_buffer,tmp_real,elsi_h%nnz_l_sp,mpi_real8,&
                 mpi_sum,elsi_h%comm_among_point,mpierr)

         if(elsi_h%my_p_row_pexsi == 0) then
            elsi_h%dm_real_ccs = tmp_real
         endif

         call elsi_deallocate(elsi_h,send_buffer,"send_buffer")
         call elsi_deallocate(elsi_h,tmp_real,"tmp_real")
      case(COMPLEX_VALUES)
         call elsi_allocate(elsi_h,send_buffer_complex,elsi_h%nnz_l_sp,&
                 "send_buffer_complex",caller)

         if(elsi_h%my_point == aux_min-1) then
            send_buffer_complex = factor_min*tmp_complex
         elseif(elsi_h%my_point == aux_max-1) then
            send_buffer_complex = factor_max*tmp_complex
         endif

         call MPI_Allreduce(send_buffer_complex,tmp_complex,elsi_h%nnz_l_sp,&
                 mpi_complex16,mpi_sum,elsi_h%comm_among_point,mpierr)

         if(elsi_h%my_p_row_pexsi == 0) then
            elsi_h%dm_complex_ccs = tmp_complex
         endif

         call elsi_deallocate(elsi_h,send_buffer_complex,"send_buffer_complex")
         call elsi_deallocate(elsi_h,tmp_complex,"tmp_complex")
      end select
   endif

   call elsi_deallocate(elsi_h,shifts,"shifts")

   ! Compute energy = Tr(H * DM)
   if(elsi_h%my_p_row_pexsi == 0) then
      local_energy = ddot(elsi_h%nnz_l_sp,elsi_h%ham_real_ccs,1,&
                        elsi_h%dm_real_ccs,1)

      call MPI_Reduce(local_energy,elsi_h%energy_hdm,1,mpi_real8,&
              mpi_sum,0,elsi_h%comm_in_pole,mpierr)
   endif

   call MPI_Bcast(elsi_h%energy_hdm,1,mpi_real8,0,elsi_h%mpi_comm,mpierr)

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished density matrix interpolation')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

   call MPI_Barrier(elsi_h%mpi_comm,mpierr)

end subroutine

!> 
!! This routine computes the energy-weighted density matrix.
!! 
subroutine elsi_compute_edm_pexsi(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   real(kind=r8)    :: mu_range
   real(kind=r8)    :: shift_width
   real(kind=r8)    :: local_energy
   real(kind=r8)    :: factor_min
   real(kind=r8)    :: factor_max
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   integer(kind=i4) :: n_shift
   integer(kind=i4) :: aux_min
   integer(kind=i4) :: aux_max
   integer(kind=i4) :: i
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: ierr
   logical          :: converged
   character*200    :: info_str

   real(kind=r8),    allocatable :: shifts(:)
   real(kind=r8),    allocatable :: tmp_real(:)
   complex(kind=r8), allocatable :: tmp_complex(:)
   real(kind=r8),    allocatable :: send_buffer(:)
   complex(kind=r8), allocatable :: send_buffer_complex(:)

   character*40, parameter :: caller = "elsi_compute_edm_pexsi"

   call elsi_get_time(elsi_h,t0)

   select case(elsi_h%matrix_data_type)
   case(REAL_VALUES)
      call f_ppexsi_calculate_edm_correction_real(elsi_h%pexsi_plan,&
              elsi_h%pexsi_options,ierr)
   case(COMPLEX_VALUES)
      call f_ppexsi_calculate_edm_correction_complex(elsi_h%pexsi_plan,&
              elsi_h%pexsi_options,ierr)
   end select

   if(ierr /= 0) then
      call elsi_stop(" Energy density matrix correction failed. Exiting...",&
              elsi_h,caller)
   endif

   ! Get energy density matrix
   select case(elsi_h%matrix_data_type)
   case(REAL_VALUES)
      call elsi_allocate(elsi_h,tmp_real,elsi_h%nnz_l_sp,&
              "tmp_real",caller)

      call f_ppexsi_retrieve_real_edm(elsi_h%pexsi_plan,tmp_real,&
              local_energy,ierr)
   case(COMPLEX_VALUES)
      call elsi_allocate(elsi_h,tmp_complex,elsi_h%nnz_l_sp,&
              "tmp_complex",caller)

      call f_ppexsi_retrieve_complex_edm(elsi_h%pexsi_plan,tmp_complex,&
              local_energy,ierr)
   end select

   if(ierr /= 0) then
      call elsi_stop(" Retrieving energy density matirx failed. Exiting...",&
              elsi_h,caller)
   endif

   ! Check convergence
   mu_range    = elsi_h%pexsi_options%muMax0-elsi_h%pexsi_options%muMin0
   n_shift     = max(10,elsi_h%n_procs/elsi_h%n_p_per_pole)
   shift_width = mu_range/(n_shift-1)
   converged   = .false.
   aux_min     = 0
   aux_max     = n_shift+1

   call elsi_allocate(elsi_h,shifts,n_shift,"shifts",caller)

   do i = 1,n_shift
      shifts(i)   = elsi_h%pexsi_options%muMin0+(i-1)*shift_width
   enddo

   do i = 1,elsi_h%pexsi_options%nPoints
      if(elsi_h%ne_vec(i) < elsi_h%n_electrons-&
         elsi_h%pexsi_options%numElectronPEXSITolerance) then
         elsi_h%pexsi_options%muMin0 = shifts(i)
         aux_min = i
      endif
   enddo

   do i = elsi_h%pexsi_options%nPoints,1,-1
      if(elsi_h%ne_vec(i) > elsi_h%n_electrons+&
         elsi_h%pexsi_options%numElectronPEXSITolerance) then
         elsi_h%pexsi_options%muMax0 = shifts(i)
         aux_max = i
      endif
   enddo

   if(elsi_h%pexsi_options%nPoints == 1) then
      ! Scale energy density matrix
      select case(elsi_h%matrix_data_type)
      case(REAL_VALUES)
         tmp_real = (elsi_h%n_electrons/elsi_h%ne_pexsi)*tmp_real
      case(COMPLEX_VALUES)
         tmp_complex = (elsi_h%n_electrons/elsi_h%ne_pexsi)*tmp_complex
      end select

      converged = .true.
      elsi_h%mu = shifts(1)
   else
      ! Safety check
      if(aux_min == 0) then
         aux_min = 1

         if(aux_max .le. aux_min) then
            aux_max = 2
         endif
      endif

      if(aux_max == n_shift+1) then
         aux_max = n_shift

         if(aux_min .ge. aux_max) then
            aux_min = n_shift-1
         endif
      endif

      do i = aux_min,aux_max
         if(abs(elsi_h%ne_vec(i)-elsi_h%n_electrons) < &
            elsi_h%pexsi_options%numElectronPEXSITolerance) then
            elsi_h%mu = shifts(i)
            converged = .true.

            select case(elsi_h%matrix_data_type)
            case(REAL_VALUES)
               call MPI_Bcast(tmp_real,elsi_h%nnz_l_sp,mpi_real8,i,&
                       elsi_h%comm_among_point,mpierr)
            case(COMPLEX_VALUES)
               call MPI_Bcast(tmp_complex,elsi_h%nnz_l_sp,mpi_complex16,&
                       i,elsi_h%comm_among_point,mpierr)
            end select

            exit
         endif
      enddo
   endif

   ! Interpolation
   if(.not. converged) then
      ! Energy density matrix
      factor_min = (elsi_h%ne_vec(aux_max)-elsi_h%n_electrons)/&
                      (elsi_h%ne_vec(aux_max)-elsi_h%ne_vec(aux_min))
      factor_max = (elsi_h%n_electrons-elsi_h%ne_vec(aux_min))/&
                      (elsi_h%ne_vec(aux_max)-elsi_h%ne_vec(aux_min))

      select case(elsi_h%matrix_data_type)
      case(REAL_VALUES)
         call elsi_allocate(elsi_h,send_buffer,elsi_h%nnz_l_sp,&
                 "send_buffer",caller)

         if(elsi_h%my_point == aux_min-1) then
            send_buffer = factor_min*tmp_real
         elseif(elsi_h%my_point == aux_max-1) then
            send_buffer = factor_max*tmp_real
         endif

         call MPI_Allreduce(send_buffer,tmp_real,elsi_h%nnz_l_sp,mpi_real8,&
                 mpi_sum,elsi_h%comm_among_point,mpierr)

         if(elsi_h%my_p_row_pexsi == 0) then
            elsi_h%dm_real_ccs = tmp_real
         endif

         call elsi_deallocate(elsi_h,send_buffer,"send_buffer")
         call elsi_deallocate(elsi_h,tmp_real,"tmp_real")

      case(COMPLEX_VALUES)
         call elsi_allocate(elsi_h,send_buffer_complex,elsi_h%nnz_l_sp,&
                 "send_buffer_complex",caller)

         if(elsi_h%my_point == aux_min-1) then
            send_buffer_complex = factor_min*tmp_complex
         elseif(elsi_h%my_point == aux_max-1) then
            send_buffer_complex = factor_max*tmp_complex
         endif

         call MPI_Allreduce(send_buffer_complex,tmp_complex,elsi_h%nnz_l_sp,&
                 mpi_complex16,mpi_sum,elsi_h%comm_among_point,mpierr)

         if(elsi_h%my_p_row_pexsi == 0) then
            elsi_h%dm_complex_ccs = tmp_complex
         endif

         call elsi_deallocate(elsi_h,send_buffer_complex,"send_buffer_complex")
         call elsi_deallocate(elsi_h,tmp_complex,"tmp_complex")
      end select
   endif

   call elsi_deallocate(elsi_h,shifts,"shifts")

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished energy density matrix calculation')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

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

   write(info_str,"(1X,' | Temperature                                    ',E10.2)") &
      elsi_h%pexsi_options%temperature
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Spectral gap                                   ',F10.3)") &
      elsi_h%pexsi_options%gap
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Number of poles                                ',I10)") &
      elsi_h%pexsi_options%numPole
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Number of mu points                            ',I10)") &
      elsi_h%pexsi_options%nPoints
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Lower bound of chemical potential              ',E10.2)") &
      elsi_h%pexsi_options%muMin0
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Upper bound of chemical potential              ',E10.2)") &
      elsi_h%pexsi_options%muMax0
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Stopping criterion of inertia counting         ',E10.2)") &
      elsi_h%pexsi_options%muInertiaTolerance
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Number of processes for symbolic factorization ',I10)") &
      elsi_h%pexsi_options%npSymbFact
   call elsi_statement_print(info_str,elsi_h)

end subroutine

end module ELSI_PEXSI
