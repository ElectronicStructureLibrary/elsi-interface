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

   use ELSI_CONSTANTS, only: REAL_VALUES,COMPLEX_VALUES,UNSET
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS
   use F_PPEXSI_INTERFACE

   implicit none

   private

   public :: elsi_set_pexsi_default
   public :: elsi_init_pexsi
   public :: elsi_solve_evp_pexsi
   public :: elsi_compute_edm_pexsi

contains

!>
!! This routine initializes PEXSI and its processor grid.
!!
subroutine elsi_init_pexsi(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: n_rows_tmp
   integer(kind=i4) :: output_id
   integer(kind=i4) :: ierr
   integer(kind=i4) :: mpierr
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_init_pexsi"

   if(e_h%n_elsi_calls == 1) then
      if(e_h%n_p_per_pole == UNSET) then
         e_h%n_p_per_pole = e_h%n_procs/(e_h%pexsi_options%numPole*&
                               e_h%pexsi_options%nPoints)
      endif

      write(info_str,"('  | MPI tasks per pole         ',I10)")&
         e_h%n_p_per_pole
      call elsi_statement_print(info_str,e_h)

      ! Set square-like process grid for selected inversion of each pole
      do n_rows_tmp = nint(sqrt(real(e_h%n_p_per_pole))),2,-1
         if(mod(e_h%n_p_per_pole,n_rows_tmp) == 0) exit
      enddo

      e_h%n_p_rows_pexsi = n_rows_tmp
      e_h%n_p_cols_pexsi = e_h%n_p_per_pole/n_rows_tmp

      ! PEXSI process grid
      e_h%my_p_col_pexsi = mod(e_h%myid,e_h%n_p_per_pole)
      e_h%my_p_row_pexsi = e_h%myid/e_h%n_p_per_pole

      ! Point parallelization
      e_h%n_p_per_point = e_h%n_procs/e_h%pexsi_options%nPoints
      e_h%my_point      = e_h%myid/e_h%n_p_per_point
      e_h%myid_point    = mod(e_h%myid,e_h%n_p_per_point)

      ! PEXSI MPI communicators
      call MPI_Comm_split(e_h%mpi_comm,e_h%my_p_col_pexsi,e_h%my_p_row_pexsi,&
              e_h%comm_among_pole,mpierr)

      call MPI_Comm_split(e_h%mpi_comm,e_h%my_p_row_pexsi,e_h%my_p_col_pexsi,&
              e_h%comm_in_pole,mpierr)

      call MPI_Comm_split(e_h%mpi_comm,e_h%myid_point,e_h%my_point,&
              e_h%comm_among_point,mpierr)

      call MPI_Comm_split(e_h%mpi_comm,e_h%my_point,e_h%myid_point,&
              e_h%comm_in_point,mpierr)

      if(.not. e_h%sparsity_ready) then
         ! Set up 1D block distribution
         e_h%n_l_cols_sp = e_h%n_basis/e_h%n_p_per_pole

         ! The last process holds all remaining columns
         if(e_h%my_p_col_pexsi == e_h%n_p_per_pole-1) then
            e_h%n_l_cols_sp = e_h%n_basis-(e_h%n_p_per_pole-1)*e_h%n_l_cols_sp
         endif
      endif

      ! Only one process outputs
      if(e_h%myid_all == 0) then
         output_id = 0
      else
         output_id = -1
      endif

      e_h%pexsi_plan = f_ppexsi_plan_initialize(e_h%mpi_comm,&
                          e_h%n_p_rows_pexsi,e_h%n_p_cols_pexsi,output_id,ierr)

      if(ierr /= 0) then
         call elsi_stop(" Initialization failed. Exiting...",e_h,caller)
      endif

      e_h%pexsi_started = .true.
   endif

end subroutine

!>
!! This routine interfaces to PEXSI.
!!
subroutine elsi_solve_evp_pexsi(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   real(kind=r8)    :: ne_drv
   real(kind=r8)    :: mu_range
   real(kind=r8)    :: shift_width
   real(kind=r8)    :: local_energy
   real(kind=r8)    :: factor_min
   real(kind=r8)    :: factor_max
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   complex(kind=r8) :: local_cmplx
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
   complex(kind=r8), allocatable :: tmp_cmplx(:)
   real(kind=r8),    allocatable :: send_buf(:)
   complex(kind=r8), allocatable :: send_buf_cmplx(:)

   real(kind=r8),    external :: ddot
!   complex(kind=r8), external :: zdotu

   character*40,  parameter :: caller = "elsi_solve_evp_pexsi"

   ! Load sparse matrices for PEXSI
   select case(e_h%matrix_data_type)
   case(REAL_VALUES)
      if(e_h%ovlp_is_unit) then
         call f_ppexsi_load_real_hs_matrix(e_h%pexsi_plan,e_h%pexsi_options,&
                 e_h%n_basis,e_h%nnz_g,e_h%nnz_l_sp,e_h%n_l_cols_sp,&
                 e_h%col_ptr_ccs,e_h%row_ind_ccs,e_h%ham_real_ccs,1,&
                 e_h%ovlp_real_ccs,ierr)
      else
         call f_ppexsi_load_real_hs_matrix(e_h%pexsi_plan,e_h%pexsi_options,&
                 e_h%n_basis,e_h%nnz_g,e_h%nnz_l_sp,e_h%n_l_cols_sp,&
                 e_h%col_ptr_ccs,e_h%row_ind_ccs,e_h%ham_real_ccs,0,&
                 e_h%ovlp_real_ccs,ierr)
      endif
   case(COMPLEX_VALUES)
      if(e_h%ovlp_is_unit) then
         call f_ppexsi_load_complex_hs_matrix(e_h%pexsi_plan,e_h%pexsi_options,&
                 e_h%n_basis,e_h%nnz_g,e_h%nnz_l_sp,e_h%n_l_cols_sp,&
                 e_h%col_ptr_ccs,e_h%row_ind_ccs,e_h%ham_cmplx_ccs,1,&
                 e_h%ovlp_cmplx_ccs,ierr)
      else
         call f_ppexsi_load_complex_hs_matrix(e_h%pexsi_plan,e_h%pexsi_options,&
                 e_h%n_basis,e_h%nnz_g,e_h%nnz_l_sp,e_h%n_l_cols_sp,&
                 e_h%col_ptr_ccs,e_h%row_ind_ccs,e_h%ham_cmplx_ccs,0,&
                 e_h%ovlp_cmplx_ccs,ierr)
      endif
   end select

   if(ierr /= 0) then
      call elsi_stop(" Failed to load matrices. Exiting...",e_h,caller)
   endif

   call elsi_statement_print("  Starting PEXSI density matrix solver",e_h)

   ! Symbolic factorization
   if(e_h%n_elsi_calls == 1) then
      call elsi_get_time(e_h,t0)

      select case(e_h%matrix_data_type)
      case(REAL_VALUES)
         call f_ppexsi_symbolic_factorize_real_symmetric_matrix(e_h%pexsi_plan,&
                 e_h%pexsi_options,ierr)

         call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(&
                 e_h%pexsi_plan,e_h%pexsi_options,ierr)
      case(COMPLEX_VALUES)
         call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(&
                 e_h%pexsi_plan,e_h%pexsi_options,ierr)

         call f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix(&
                 e_h%pexsi_plan,e_h%pexsi_options,e_h%ovlp_cmplx_ccs,ierr)
      end select

      call elsi_get_time(e_h,t1)

      write(info_str,"('  Finished symbolic factorization')")
      call elsi_statement_print(info_str,e_h)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_statement_print(info_str,e_h)
   endif

   if(ierr /= 0) then
      call elsi_stop(" Symbolic factorization failed. Exiting...",e_h,caller)
   endif

   ! Inertia counting
   call elsi_get_time(e_h,t0)

   n_iner_steps = 0
   mu_range = e_h%pexsi_options%muMax0-e_h%pexsi_options%muMin0
   n_shift = max(10,e_h%n_procs/e_h%n_p_per_pole)

   call elsi_allocate(e_h,shifts,n_shift,"shifts",caller)
   call elsi_allocate(e_h,inertias,n_shift,"inertias",caller)
   call elsi_allocate(e_h,ne_lower,n_shift,"ne_lower",caller)
   call elsi_allocate(e_h,ne_upper,n_shift,"ne_upper",caller)

   if(.not. e_h%spin_is_set) then
      if(e_h%n_spins == 2) then
         e_h%spin_degen = 1.0_r8
      else
         e_h%spin_degen = 2.0_r8
      endif
   endif

   do while(n_iner_steps < 10 .and. &
            mu_range > e_h%pexsi_options%muInertiaTolerance)
      n_iner_steps = n_iner_steps+1

      shift_width = mu_range/(n_shift-1)

      ne_lower = 0.0_r8
      ne_upper = e_h%n_basis*e_h%spin_degen

      do i = 1,n_shift
         shifts(i)   = e_h%pexsi_options%muMin0+(i-1)*shift_width
      enddo

      select case(e_h%matrix_data_type)
      case(REAL_VALUES)
         call f_ppexsi_inertia_count_real_matrix(e_h%pexsi_plan,&
                 e_h%pexsi_options,n_shift,shifts,inertias,ierr)
      case(COMPLEX_VALUES)
         call f_ppexsi_inertia_count_complex_matrix(e_h%pexsi_plan,&
                 e_h%pexsi_options,n_shift,shifts,inertias,ierr)
      end select

      inertias = inertias*e_h%spin_degen*e_h%i_weight

      ! Get global inertias
      if(e_h%n_spins*e_h%n_kpts > 1) then
         call elsi_allocate(e_h,send_buf,n_shift,"send_buf",caller)

         if(e_h%myid == 0) then
            send_buf = inertias
         else
            send_buf = 0.0_r8
         endif

         call MPI_Allreduce(send_buf,inertias,n_shift,mpi_real8,mpi_sum,&
                 e_h%mpi_comm_all,mpierr)

         call elsi_deallocate(e_h,send_buf,"send_buf")
      endif

      idx = ceiling(3*e_h%pexsi_options%temperature/shift_width)

      do i = idx+1,n_shift
         ne_lower(i) = 0.5_r8*(inertias(i-idx)+inertias(i))
         ne_upper(i-idx) = ne_lower(i)
      enddo

      aux_min = 1
      aux_max = n_shift

      do i = 2,n_shift-1
         if(ne_upper(i) < e_h%n_electrons .and. &
            ne_upper(i+1) >= e_h%n_electrons)  then
            aux_min = i
         endif

         if(ne_lower(i) > e_h%n_electrons .and. &
            ne_lower(i-1) <= e_h%n_electrons) then
            aux_max = i
         endif
      enddo

      if(aux_min == 1 .and. aux_max == n_shift) then
         exit
      else
         e_h%pexsi_options%muMin0 = shifts(aux_min)
         e_h%pexsi_options%muMax0 = shifts(aux_max)
         mu_range = e_h%pexsi_options%muMax0-e_h%pexsi_options%muMin0
      endif
   enddo

   call elsi_deallocate(e_h,shifts,"shifts")
   call elsi_deallocate(e_h,inertias,"inertias")
   call elsi_deallocate(e_h,ne_lower,"ne_lower")
   call elsi_deallocate(e_h,ne_upper,"ne_upper")

   call elsi_get_time(e_h,t1)

   if(n_iner_steps > 0) then
      write(info_str,"('  Finished inertia counting')")
      call elsi_statement_print(info_str,e_h)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_statement_print(info_str,e_h)
   endif

   if(ierr /= 0) then
      call elsi_stop(" Inertia counting failed. Exiting...",e_h,caller)
   endif

   ! Fermi operator expansion
   call elsi_get_time(e_h,t0)

   shift_width = mu_range/(e_h%pexsi_options%nPoints+1)

   call elsi_allocate(e_h,shifts,e_h%pexsi_options%nPoints,"shifts",caller)

   do i = 1,e_h%pexsi_options%nPoints
      shifts(i) = e_h%pexsi_options%muMin0+i*shift_width
   enddo

   do i = 1,e_h%pexsi_options%nPoints
      e_h%mu = shifts(i)

      if(e_h%my_point == i-1) then
         select case(e_h%matrix_data_type)
         case(REAL_VALUES)
            call f_ppexsi_calculate_fermi_operator_real3(e_h%pexsi_plan,&
                    e_h%pexsi_options,e_h%mu,e_h%n_electrons,e_h%ne_pexsi,&
                    ne_drv,ierr)
         case(COMPLEX_VALUES)
            call f_ppexsi_calculate_fermi_operator_complex(e_h%pexsi_plan,&
                    e_h%pexsi_options,e_h%mu,e_h%n_electrons,e_h%ne_pexsi,&
                    ne_drv,ierr)
         end select
      endif
   enddo

   call elsi_allocate(e_h,send_buf,e_h%pexsi_options%nPoints,"send_buf",caller)

   if(e_h%n_elsi_calls == 1) then
      call elsi_allocate(e_h,e_h%ne_vec,e_h%pexsi_options%nPoints,"ne_vec",&
              caller)
   endif

   send_buf(e_h%my_point+1) = e_h%ne_pexsi*e_h%i_weight

   call MPI_Allreduce(send_buf,e_h%ne_vec,e_h%pexsi_options%nPoints,mpi_real8,&
           mpi_sum,e_h%comm_among_point,mpierr)

   ! Get global number of electrons
   if(e_h%n_spins*e_h%n_kpts > 1) then
      if(e_h%myid == 0) then
         send_buf = e_h%ne_vec
      else
         send_buf = 0.0_r8
      endif

      call MPI_Allreduce(send_buf,e_h%ne_vec,e_h%pexsi_options%nPoints,&
              mpi_real8,mpi_sum,e_h%mpi_comm_all,mpierr)
   endif

   call elsi_deallocate(e_h,send_buf,"send_buf")

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished Fermi operator calculation')")
   call elsi_statement_print(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,e_h)

   if(ierr /= 0) then
      call elsi_stop(" Fermi operator calculation failed. Exiting...",e_h,&
              caller)
   endif

   ! Get density matrix
   call elsi_get_time(e_h,t0)

   select case(e_h%matrix_data_type)
   case(REAL_VALUES)
      call elsi_allocate(e_h,tmp_real,e_h%nnz_l_sp,"tmp_real",caller)

      call f_ppexsi_retrieve_real_dm(e_h%pexsi_plan,tmp_real,local_energy,ierr)
   case(COMPLEX_VALUES)
      call elsi_allocate(e_h,tmp_cmplx,e_h%nnz_l_sp,"tmp_cmplx",caller)

      call f_ppexsi_retrieve_complex_dm(e_h%pexsi_plan,tmp_cmplx,local_energy,&
              ierr)
   end select

   if(ierr /= 0) then
      call elsi_stop(" Failed to get density matirx. Exiting...",e_h,caller)
   endif

   ! Check convergence
   converged = .false.
   aux_min = 0
   aux_max = e_h%pexsi_options%nPoints+1

   do i = 1,e_h%pexsi_options%nPoints
      if(e_h%ne_vec(i) < e_h%n_electrons-&
         e_h%pexsi_options%numElectronPEXSITolerance) then
         e_h%pexsi_options%muMin0 = shifts(i)
         aux_min = i
      endif
   enddo

   do i = e_h%pexsi_options%nPoints,1,-1
      if(e_h%ne_vec(i) > e_h%n_electrons+&
         e_h%pexsi_options%numElectronPEXSITolerance) then
         e_h%pexsi_options%muMax0 = shifts(i)
         aux_max = i
      endif
   enddo

   if(e_h%pexsi_options%nPoints == 1) then
      ! Scale density matrix
      select case(e_h%matrix_data_type)
      case(REAL_VALUES)
         tmp_real = (e_h%n_electrons/e_h%ne_pexsi)*tmp_real
      case(COMPLEX_VALUES)
         tmp_cmplx = (e_h%n_electrons/e_h%ne_pexsi)*tmp_cmplx
      end select

      converged = .true.
      e_h%mu = shifts(1)
   else
      ! Safety check
      if(aux_min == 0) then
         aux_min = 1

         if(aux_max <= aux_min) then
            aux_max = 2
         endif
      endif

      if(aux_max == e_h%pexsi_options%nPoints+1) then
         aux_max = e_h%pexsi_options%nPoints

         if(aux_min >= aux_max) then
            aux_min = e_h%pexsi_options%nPoints-1
         endif
      endif

      do i = aux_min,aux_max
         if(abs(e_h%ne_vec(i)-e_h%n_electrons) < &
            e_h%pexsi_options%numElectronPEXSITolerance) then
            e_h%mu = shifts(i)
            converged = .true.

            select case(e_h%matrix_data_type)
            case(REAL_VALUES)
               call MPI_Bcast(tmp_real,e_h%nnz_l_sp,mpi_real8,i,&
                       e_h%comm_among_point,mpierr)
            case(COMPLEX_VALUES)
               call MPI_Bcast(tmp_cmplx,e_h%nnz_l_sp,mpi_complex16,i,&
                       e_h%comm_among_point,mpierr)
            end select

            exit
         endif
      enddo
   endif

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Chemical potential
      e_h%mu = shifts(aux_min)+(e_h%n_electrons-e_h%ne_vec(aux_min))/&
                  (e_h%ne_vec(aux_max)-e_h%ne_vec(aux_min))*&
                  (shifts(aux_max)-shifts(aux_min))

      ! Density matrix
      factor_min = (e_h%ne_vec(aux_max)-e_h%n_electrons)/&
                      (e_h%ne_vec(aux_max)-e_h%ne_vec(aux_min))
      factor_max = (e_h%n_electrons-e_h%ne_vec(aux_min))/&
                      (e_h%ne_vec(aux_max)-e_h%ne_vec(aux_min))

      select case(e_h%matrix_data_type)
      case(REAL_VALUES)
         call elsi_allocate(e_h,send_buf,e_h%nnz_l_sp,"send_buf",caller)

         if(e_h%my_point == aux_min-1) then
            send_buf = factor_min*tmp_real
         elseif(e_h%my_point == aux_max-1) then
            send_buf = factor_max*tmp_real
         endif

         call MPI_Allreduce(send_buf,tmp_real,e_h%nnz_l_sp,mpi_real8,mpi_sum,&
                 e_h%comm_among_point,mpierr)

         if(e_h%my_p_row_pexsi == 0) then
            e_h%dm_real_ccs = tmp_real
         endif

         call elsi_deallocate(e_h,send_buf,"send_buf")
         call elsi_deallocate(e_h,tmp_real,"tmp_real")
      case(COMPLEX_VALUES)
         call elsi_allocate(e_h,send_buf_cmplx,e_h%nnz_l_sp,"send_buf_cmplx",&
                 caller)

         if(e_h%my_point == aux_min-1) then
            send_buf_cmplx = factor_min*tmp_cmplx
         elseif(e_h%my_point == aux_max-1) then
            send_buf_cmplx = factor_max*tmp_cmplx
         endif

         call MPI_Allreduce(send_buf_cmplx,tmp_cmplx,e_h%nnz_l_sp,&
                 mpi_complex16,mpi_sum,e_h%comm_among_point,mpierr)

         if(e_h%my_p_row_pexsi == 0) then
            e_h%dm_cmplx_ccs = tmp_cmplx
         endif

         call elsi_deallocate(e_h,send_buf_cmplx,"send_buf_cmplx")
         call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")
      end select
   endif

   call elsi_deallocate(e_h,shifts,"shifts")

   ! Compute energy = Tr(H*DM)
   if(e_h%my_p_row_pexsi == 0) then
      select case(e_h%matrix_data_type)
      case(REAL_VALUES)
         local_energy = ddot(e_h%nnz_l_sp,e_h%ham_real_ccs,1,e_h%dm_real_ccs,1)
      case(COMPLEX_VALUES)
         ! The following lines are equivalent to "zdotu" in LAPACK,
         ! which is simple but problematic in some version of LAPACK
         local_cmplx = (0.0_r8,0.0_r8)

         do i = 1,e_h%nnz_l_sp
            local_cmplx = local_cmplx+e_h%ham_cmplx_ccs(i)*e_h%dm_cmplx_ccs(i)
         enddo

         local_energy = dble(local_cmplx)
      end select

      call MPI_Reduce(local_energy,e_h%energy_hdm,1,mpi_real8,mpi_sum,0,&
              e_h%comm_in_pole,mpierr)
   endif

   call MPI_Bcast(e_h%energy_hdm,1,mpi_real8,0,e_h%mpi_comm,mpierr)

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished density matrix correction')")
   call elsi_statement_print(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,e_h)

   call MPI_Barrier(e_h%mpi_comm,mpierr)

end subroutine

!>
!! This routine computes the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_pexsi(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   real(kind=r8)    :: mu_range
   real(kind=r8)    :: shift_width
   real(kind=r8)    :: local_energy
   real(kind=r8)    :: factor_min
   real(kind=r8)    :: factor_max
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   integer(kind=i4) :: aux_min
   integer(kind=i4) :: aux_max
   integer(kind=i4) :: i
   integer(kind=i4) :: mpierr
   integer(kind=i4) :: ierr
   logical          :: converged
   character*200    :: info_str

   real(kind=r8),    allocatable :: shifts(:)
   real(kind=r8),    allocatable :: tmp_real(:)
   complex(kind=r8), allocatable :: tmp_cmplx(:)
   real(kind=r8),    allocatable :: send_buf(:)
   complex(kind=r8), allocatable :: send_buf_cmplx(:)

   character*40, parameter :: caller = "elsi_compute_edm_pexsi"

   call elsi_get_time(e_h,t0)

   select case(e_h%matrix_data_type)
   case(REAL_VALUES)
      call f_ppexsi_calculate_edm_correction_real(e_h%pexsi_plan,&
              e_h%pexsi_options,ierr)
   case(COMPLEX_VALUES)
      call f_ppexsi_calculate_edm_correction_complex(e_h%pexsi_plan,&
              e_h%pexsi_options,ierr)
   end select

   if(ierr /= 0) then
      call elsi_stop(" Energy density matrix correction failed. Exiting...",&
              e_h,caller)
   endif

   ! Get energy density matrix
   select case(e_h%matrix_data_type)
   case(REAL_VALUES)
      call elsi_allocate(e_h,tmp_real,e_h%nnz_l_sp,"tmp_real",caller)

      call f_ppexsi_retrieve_real_edm(e_h%pexsi_plan,tmp_real,local_energy,ierr)
   case(COMPLEX_VALUES)
      call elsi_allocate(e_h,tmp_cmplx,e_h%nnz_l_sp,"tmp_cmplx",caller)

      call f_ppexsi_retrieve_complex_edm(e_h%pexsi_plan,tmp_cmplx,local_energy,&
              ierr)
   end select

   if(ierr /= 0) then
      call elsi_stop(" Failed to get energy density matirx failed. Exiting...",&
              e_h,caller)
   endif

   ! Check convergence
   mu_range    = e_h%pexsi_options%muMax0-e_h%pexsi_options%muMin0
   shift_width = mu_range/(e_h%pexsi_options%nPoints+1)
   converged   = .false.
   aux_min     = 0
   aux_max     = e_h%pexsi_options%nPoints+1

   call elsi_allocate(e_h,shifts,e_h%pexsi_options%nPoints,"shifts",caller)

   do i = 1,e_h%pexsi_options%nPoints
      shifts(i) = e_h%pexsi_options%muMin0+i*shift_width
   enddo

   do i = 1,e_h%pexsi_options%nPoints
      if(e_h%ne_vec(i) < e_h%n_electrons-&
         e_h%pexsi_options%numElectronPEXSITolerance) then
         e_h%pexsi_options%muMin0 = shifts(i)
         aux_min = i
      endif
   enddo

   do i = e_h%pexsi_options%nPoints,1,-1
      if(e_h%ne_vec(i) > e_h%n_electrons+&
         e_h%pexsi_options%numElectronPEXSITolerance) then
         e_h%pexsi_options%muMax0 = shifts(i)
         aux_max = i
      endif
   enddo

   if(e_h%pexsi_options%nPoints == 1) then
      ! Scale energy density matrix
      select case(e_h%matrix_data_type)
      case(REAL_VALUES)
         tmp_real = (e_h%n_electrons/e_h%ne_pexsi)*tmp_real
      case(COMPLEX_VALUES)
         tmp_cmplx = (e_h%n_electrons/e_h%ne_pexsi)*tmp_cmplx
      end select

      converged = .true.
      e_h%mu = shifts(1)
   else
      ! Safety check
      if(aux_min == 0) then
         aux_min = 1

         if(aux_max <= aux_min) then
            aux_max = 2
         endif
      endif

      if(aux_max == e_h%pexsi_options%nPoints+1) then
         aux_max = e_h%pexsi_options%nPoints

         if(aux_min >= aux_max) then
            aux_min = e_h%pexsi_options%nPoints-1
         endif
      endif

      do i = aux_min,aux_max
         if(abs(e_h%ne_vec(i)-e_h%n_electrons) < &
            e_h%pexsi_options%numElectronPEXSITolerance) then
            e_h%mu = shifts(i)
            converged = .true.

            select case(e_h%matrix_data_type)
            case(REAL_VALUES)
               call MPI_Bcast(tmp_real,e_h%nnz_l_sp,mpi_real8,i,&
                       e_h%comm_among_point,mpierr)
            case(COMPLEX_VALUES)
               call MPI_Bcast(tmp_cmplx,e_h%nnz_l_sp,mpi_complex16,i,&
                       e_h%comm_among_point,mpierr)
            end select

            exit
         endif
      enddo
   endif

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Energy density matrix
      factor_min = (e_h%ne_vec(aux_max)-e_h%n_electrons)/&
                      (e_h%ne_vec(aux_max)-e_h%ne_vec(aux_min))
      factor_max = (e_h%n_electrons-e_h%ne_vec(aux_min))/&
                      (e_h%ne_vec(aux_max)-e_h%ne_vec(aux_min))

      select case(e_h%matrix_data_type)
      case(REAL_VALUES)
         call elsi_allocate(e_h,send_buf,e_h%nnz_l_sp,"send_buf",caller)

         if(e_h%my_point == aux_min-1) then
            send_buf = factor_min*tmp_real
         elseif(e_h%my_point == aux_max-1) then
            send_buf = factor_max*tmp_real
         endif

         call MPI_Allreduce(send_buf,tmp_real,e_h%nnz_l_sp,mpi_real8,mpi_sum,&
                 e_h%comm_among_point,mpierr)

         if(e_h%my_p_row_pexsi == 0) then
            e_h%dm_real_ccs = tmp_real
         endif

         call elsi_deallocate(e_h,send_buf,"send_buf")
         call elsi_deallocate(e_h,tmp_real,"tmp_real")
      case(COMPLEX_VALUES)
         call elsi_allocate(e_h,send_buf_cmplx,e_h%nnz_l_sp,"send_buf_cmplx",&
                 caller)

         if(e_h%my_point == aux_min-1) then
            send_buf_cmplx = factor_min*tmp_cmplx
         elseif(e_h%my_point == aux_max-1) then
            send_buf_cmplx = factor_max*tmp_cmplx
         endif

         call MPI_Allreduce(send_buf_cmplx,tmp_cmplx,e_h%nnz_l_sp,&
                 mpi_complex16,mpi_sum,e_h%comm_among_point,mpierr)

         if(e_h%my_p_row_pexsi == 0) then
            e_h%dm_cmplx_ccs = tmp_cmplx
         endif

         call elsi_deallocate(e_h,send_buf_cmplx,"send_buf_cmplx")
         call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")
      end select
   endif

   call elsi_deallocate(e_h,shifts,"shifts")

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished energy density matrix calculation')")
   call elsi_statement_print(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,e_h)

end subroutine

!>
!! This routine sets default PEXSI parameters.
!!
subroutine elsi_set_pexsi_default(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character*40, parameter :: caller = "elsi_set_pexsi_default"

   ! Use the PEXSI Default options
   call f_ppexsi_set_default_options(e_h%pexsi_options)

   ! Use 1 process in symbolic factorization
   e_h%pexsi_options%npSymbFact = 1

end subroutine

end module ELSI_PEXSI
