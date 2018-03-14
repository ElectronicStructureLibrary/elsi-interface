! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides interfaces to PEXSI.
!!
module ELSI_PEXSI

   use ELSI_CONSTANTS,     only: UNSET
   use ELSI_DATATYPE,      only: elsi_handle
   use ELSI_IO,            only: elsi_say
   use ELSI_MALLOC,        only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,           only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8,&
                                 mpi_complex16
   use ELSI_PRECISION,     only: r8,i4
   use ELSI_TIMINGS,       only: elsi_get_time
   use F_PPEXSI_INTERFACE, only: f_ppexsi_plan_initialize,&
                                 f_ppexsi_set_default_options,&
                                 f_ppexsi_load_real_hs_matrix,&
                                 f_ppexsi_load_complex_hs_matrix,&
                                 f_ppexsi_symbolic_factorize_real_symmetric_matrix,&
                                 f_ppexsi_symbolic_factorize_complex_symmetric_matrix,&
                                 f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix,&
                                 f_ppexsi_inertia_count_real_matrix,&
                                 f_ppexsi_inertia_count_complex_matrix,&
                                 f_ppexsi_calculate_fermi_operator_real3,&
                                 f_ppexsi_calculate_fermi_operator_complex,&
                                 f_ppexsi_calculate_edm_correction_real,&
                                 f_ppexsi_calculate_edm_correction_complex,&
                                 f_ppexsi_retrieve_real_dm,&
                                 f_ppexsi_retrieve_complex_dm,&
                                 f_ppexsi_retrieve_real_edm,&
                                 f_ppexsi_retrieve_complex_edm,&
                                 f_ppexsi_plan_finalize

   implicit none

   private

   public :: elsi_set_pexsi_default
   public :: elsi_init_pexsi
   public :: elsi_solve_evp_pexsi_real
   public :: elsi_compute_edm_pexsi_real
   public :: elsi_solve_evp_pexsi_cmplx
   public :: elsi_compute_edm_pexsi_cmplx

contains

!>
!! This routine initializes PEXSI and its processor grid.
!!
subroutine elsi_init_pexsi(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   integer(kind=i4)   :: n_rows_tmp
   integer(kind=i4)   :: output_id
   integer(kind=i4)   :: ierr
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_init_pexsi"

   if(e_h%handle_ready) then
      e_h%handle_changed = .true.
   endif

   if(e_h%n_elsi_calls == 1) then
      e_h%pexsi_options%spin = e_h%spin_degen

      if(e_h%pexsi_np_per_pole == UNSET) then
         e_h%pexsi_np_per_pole = e_h%n_procs/(e_h%pexsi_options%numPole*&
                                    e_h%pexsi_options%nPoints)
      endif

      write(info_str,"('  | MPI tasks per pole         ',I10)")&
         e_h%pexsi_np_per_pole
      call elsi_say(e_h,info_str)

      ! Set square-like process grid for selected inversion of each pole
      do n_rows_tmp = nint(sqrt(real(e_h%pexsi_np_per_pole))),2,-1
         if(mod(e_h%pexsi_np_per_pole,n_rows_tmp) == 0) exit
      enddo

      e_h%pexsi_n_prow = n_rows_tmp
      e_h%pexsi_n_pcol = e_h%pexsi_np_per_pole/n_rows_tmp

      ! PEXSI process grid
      e_h%pexsi_my_pcol = mod(e_h%myid,e_h%pexsi_np_per_pole)
      e_h%pexsi_my_prow = e_h%myid/e_h%pexsi_np_per_pole

      ! Point parallelization
      e_h%pexsi_np_per_point = e_h%n_procs/e_h%pexsi_options%nPoints
      e_h%pexsi_my_point     = e_h%myid/e_h%pexsi_np_per_point
      e_h%pexsi_myid_point   = mod(e_h%myid,e_h%pexsi_np_per_point)

      ! PEXSI MPI communicators
      call MPI_Comm_split(e_h%mpi_comm,e_h%pexsi_my_pcol,e_h%pexsi_my_prow,&
              e_h%pexsi_comm_among_pole,ierr)

      call elsi_check_mpi(e_h,"MPI_Comm_split",ierr,caller)

      call MPI_Comm_split(e_h%mpi_comm,e_h%pexsi_my_prow,e_h%pexsi_my_pcol,&
              e_h%pexsi_comm_in_pole,ierr)

      call elsi_check_mpi(e_h,"MPI_Comm_split",ierr,caller)

      call MPI_Comm_split(e_h%mpi_comm,e_h%pexsi_myid_point,e_h%pexsi_my_point,&
              e_h%pexsi_comm_among_point,ierr)

      call elsi_check_mpi(e_h,"MPI_Comm_split",ierr,caller)

      call MPI_Comm_split(e_h%mpi_comm,e_h%pexsi_my_point,e_h%pexsi_myid_point,&
              e_h%pexsi_comm_in_point,ierr)

      call elsi_check_mpi(e_h,"MPI_Comm_split",ierr,caller)

      if(.not. e_h%pexsi_csc_ready) then
         ! Set up 1D block distribution
         e_h%n_lcol_sp1 = e_h%n_basis/e_h%pexsi_np_per_pole

         ! The last process holds all remaining columns
         if(e_h%pexsi_my_pcol == e_h%pexsi_np_per_pole-1) then
            e_h%n_lcol_sp1 = e_h%n_basis-&
                                (e_h%pexsi_np_per_pole-1)*e_h%n_lcol_sp1
         endif
      endif

      ! Only one process outputs
      if(e_h%myid_all == 0) then
         output_id = 0
      else
         output_id = -1
      endif

      e_h%pexsi_plan = f_ppexsi_plan_initialize(e_h%mpi_comm,e_h%pexsi_n_prow,&
                          e_h%pexsi_n_pcol,output_id,ierr)

      if(ierr /= 0) then
         call elsi_stop(" Initialization failed.",e_h,caller)
      endif

      if(e_h%n_lcol_sp == UNSET) then
         e_h%n_lcol_sp = e_h%n_lcol_sp1
      endif

      e_h%pexsi_started = .true.
   endif

end subroutine

!>
!! This routine interfaces to PEXSI.
!!
subroutine elsi_solve_evp_pexsi_real(e_h,ham,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: ham(e_h%nnz_l_sp1)
   real(kind=r8),     intent(inout) :: ovlp(e_h%nnz_l_sp1)
   real(kind=r8),     intent(inout) :: dm(e_h%nnz_l_sp1)

   real(kind=r8)      :: ne_drv
   real(kind=r8)      :: mu_range
   real(kind=r8)      :: shift_width
   real(kind=r8)      :: local_energy
   real(kind=r8)      :: factor_min
   real(kind=r8)      :: factor_max
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   integer(kind=i4)   :: n_iner_steps
   integer(kind=i4)   :: n_shift
   integer(kind=i4)   :: aux_min
   integer(kind=i4)   :: aux_max
   integer(kind=i4)   :: i
   integer(kind=i4)   :: idx
   integer(kind=i4)   :: ierr
   logical            :: converged
   character(len=200) :: info_str

   real(kind=r8), allocatable :: shifts(:)
   real(kind=r8), allocatable :: inertias(:)
   real(kind=r8), allocatable :: ne_lower(:)
   real(kind=r8), allocatable :: ne_upper(:)
   real(kind=r8), allocatable :: tmp_real(:)
   real(kind=r8), allocatable :: send_buf(:)

   real(kind=r8), external :: ddot

   character(len=40), parameter :: caller = "elsi_solve_evp_pexsi_real"

   ! Load sparse matrices for PEXSI
   if(e_h%ovlp_is_unit) then
      call f_ppexsi_load_real_hs_matrix(e_h%pexsi_plan,e_h%pexsi_options,&
              e_h%n_basis,e_h%nnz_g,e_h%nnz_l_sp1,e_h%n_lcol_sp1,&
              e_h%col_ptr_pexsi,e_h%row_ind_pexsi,ham,1,ovlp,ierr)
   else
      call f_ppexsi_load_real_hs_matrix(e_h%pexsi_plan,e_h%pexsi_options,&
              e_h%n_basis,e_h%nnz_g,e_h%nnz_l_sp1,e_h%n_lcol_sp1,&
              e_h%col_ptr_pexsi,e_h%row_ind_pexsi,ham,0,ovlp,ierr)
   endif

   if(ierr /= 0) then
      call elsi_stop(" Failed to load matrices.",e_h,caller)
   endif

   call elsi_say(e_h,"  Starting PEXSI density matrix solver")

   ! Symbolic factorization
   if(e_h%n_elsi_calls == 1) then
      call elsi_get_time(t0)

      call f_ppexsi_symbolic_factorize_real_symmetric_matrix(e_h%pexsi_plan,&
              e_h%pexsi_options,ierr)

      call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(e_h%pexsi_plan,&
              e_h%pexsi_options,ierr)

      call elsi_get_time(t1)

      write(info_str,"('  Finished symbolic factorization')")
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_say(e_h,info_str)
   endif

   if(ierr /= 0) then
      call elsi_stop(" Symbolic factorization failed.",e_h,caller)
   endif

   ! Inertia counting
   call elsi_get_time(t0)

   n_iner_steps = 0
   mu_range     = e_h%pexsi_options%muMax0-e_h%pexsi_options%muMin0
   n_shift      = max(10,e_h%n_procs/e_h%pexsi_np_per_pole)

   call elsi_allocate(e_h,shifts,n_shift,"shifts",caller)
   call elsi_allocate(e_h,inertias,n_shift,"inertias",caller)
   call elsi_allocate(e_h,ne_lower,n_shift,"ne_lower",caller)
   call elsi_allocate(e_h,ne_upper,n_shift,"ne_upper",caller)

   do while(n_iner_steps < 10 .and. &
            mu_range > e_h%pexsi_options%muInertiaTolerance)
      n_iner_steps = n_iner_steps+1

      shift_width = mu_range/(n_shift-1)

      ne_lower = 0.0_r8
      ne_upper = e_h%n_basis*e_h%spin_degen

      do i = 1,n_shift
         shifts(i) = e_h%pexsi_options%muMin0+(i-1)*shift_width
      enddo

      call f_ppexsi_inertia_count_real_matrix(e_h%pexsi_plan,e_h%pexsi_options,&
              n_shift,shifts,inertias,ierr)

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
                 e_h%mpi_comm_all,ierr)

         call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

         call elsi_deallocate(e_h,send_buf,"send_buf")
      endif

      idx = ceiling(3*e_h%pexsi_options%temperature/shift_width)

      do i = idx+1,n_shift
         ne_lower(i)     = 0.5_r8*(inertias(i-idx)+inertias(i))
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

   call elsi_get_time(t1)

   if(n_iner_steps > 0) then
      write(info_str,"('  Finished inertia counting')")
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_say(e_h,info_str)
   endif

   if(ierr /= 0) then
      call elsi_stop(" Inertia counting failed.",e_h,caller)
   endif

   ! Fermi operator expansion
   call elsi_get_time(t0)

   shift_width = mu_range/(e_h%pexsi_options%nPoints+1)

   call elsi_allocate(e_h,shifts,e_h%pexsi_options%nPoints,"shifts",caller)

   do i = 1,e_h%pexsi_options%nPoints
      shifts(i) = e_h%pexsi_options%muMin0+i*shift_width
   enddo

   do i = 1,e_h%pexsi_options%nPoints
      e_h%mu = shifts(i)

      if(e_h%pexsi_my_point == i-1) then
         call f_ppexsi_calculate_fermi_operator_real3(e_h%pexsi_plan,&
                 e_h%pexsi_options,e_h%mu,e_h%n_electrons,e_h%pexsi_ne,ne_drv,&
                 ierr)
      endif
   enddo

   call elsi_allocate(e_h,send_buf,e_h%pexsi_options%nPoints,"send_buf",caller)

   if(e_h%n_elsi_calls == 1) then
      call elsi_allocate(e_h,e_h%ne_vec_pexsi,e_h%pexsi_options%nPoints,&
              "ne_vec_pexsi",caller)
   endif

   send_buf(e_h%pexsi_my_point+1) = e_h%pexsi_ne*e_h%i_weight

   call MPI_Allreduce(send_buf,e_h%ne_vec_pexsi,e_h%pexsi_options%nPoints,&
           mpi_real8,mpi_sum,e_h%pexsi_comm_among_point,ierr)

   call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

   ! Get global number of electrons
   if(e_h%n_spins*e_h%n_kpts > 1) then
      if(e_h%myid == 0) then
         send_buf = e_h%ne_vec_pexsi
      else
         send_buf = 0.0_r8
      endif

      call MPI_Allreduce(send_buf,e_h%ne_vec_pexsi,e_h%pexsi_options%nPoints,&
              mpi_real8,mpi_sum,e_h%mpi_comm_all,ierr)

      call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)
   endif

   call elsi_deallocate(e_h,send_buf,"send_buf")

   call elsi_get_time(t1)

   write(info_str,"('  Finished Fermi operator calculation')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

   if(ierr /= 0) then
      call elsi_stop(" Fermi operator calculation failed.",e_h,caller)
   endif

   ! Get density matrix
   call elsi_get_time(t0)

   call elsi_allocate(e_h,tmp_real,e_h%nnz_l_sp1,"tmp_real",caller)

   call f_ppexsi_retrieve_real_dm(e_h%pexsi_plan,tmp_real,local_energy,ierr)

   if(ierr /= 0) then
      call elsi_stop(" Failed to get density matirx.",e_h,caller)
   endif

   ! Check convergence
   converged = .false.
   aux_min   = 0
   aux_max   = e_h%pexsi_options%nPoints+1

   do i = 1,e_h%pexsi_options%nPoints
      if(e_h%ne_vec_pexsi(i) < e_h%n_electrons-&
         e_h%pexsi_options%numElectronPEXSITolerance) then
         e_h%pexsi_options%muMin0 = shifts(i)
         aux_min = i
      endif
   enddo

   do i = e_h%pexsi_options%nPoints,1,-1
      if(e_h%ne_vec_pexsi(i) > e_h%n_electrons+&
         e_h%pexsi_options%numElectronPEXSITolerance) then
         e_h%pexsi_options%muMax0 = shifts(i)
         aux_max = i
      endif
   enddo

   if(e_h%pexsi_options%nPoints == 1) then
      ! Scale density matrix
      tmp_real  = (e_h%n_electrons/e_h%pexsi_ne)*tmp_real
      converged = .true.
      e_h%mu    = shifts(1)
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
         if(abs(e_h%ne_vec_pexsi(i)-e_h%n_electrons) < &
            e_h%pexsi_options%numElectronPEXSITolerance) then
            e_h%mu    = shifts(i)
            converged = .true.

            call MPI_Bcast(tmp_real,e_h%nnz_l_sp1,mpi_real8,i-1,&
                    e_h%pexsi_comm_among_point,ierr)

            call elsi_check_mpi(e_h,"MPI_Bcast",ierr,caller)

            exit
         endif
      enddo
   endif

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Chemical potential
      e_h%mu = shifts(aux_min)+(e_h%n_electrons-e_h%ne_vec_pexsi(aux_min))/&
                  (e_h%ne_vec_pexsi(aux_max)-e_h%ne_vec_pexsi(aux_min))*&
                  (shifts(aux_max)-shifts(aux_min))

      ! Density matrix
      factor_min = (e_h%ne_vec_pexsi(aux_max)-e_h%n_electrons)/&
                      (e_h%ne_vec_pexsi(aux_max)-e_h%ne_vec_pexsi(aux_min))
      factor_max = (e_h%n_electrons-e_h%ne_vec_pexsi(aux_min))/&
                      (e_h%ne_vec_pexsi(aux_max)-e_h%ne_vec_pexsi(aux_min))

      call elsi_allocate(e_h,send_buf,e_h%nnz_l_sp1,"send_buf",caller)

      if(e_h%pexsi_my_point == aux_min-1) then
         send_buf = factor_min*tmp_real
      elseif(e_h%pexsi_my_point == aux_max-1) then
         send_buf = factor_max*tmp_real
      endif

      call MPI_Allreduce(send_buf,tmp_real,e_h%nnz_l_sp1,mpi_real8,mpi_sum,&
              e_h%pexsi_comm_among_point,ierr)

      call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(e_h,send_buf,"send_buf")
   endif

   if(e_h%pexsi_my_prow == 0) then
      dm = tmp_real
   endif

   call elsi_deallocate(e_h,tmp_real,"tmp_real")
   call elsi_deallocate(e_h,shifts,"shifts")

   ! Compute energy = Tr(H*DM)
   if(e_h%pexsi_my_prow == 0) then
      local_energy = ddot(e_h%nnz_l_sp1,ham,1,dm,1)

      call MPI_Reduce(local_energy,e_h%energy_hdm,1,mpi_real8,mpi_sum,0,&
              e_h%pexsi_comm_in_pole,ierr)

      call elsi_check_mpi(e_h,"MPI_Reduce",ierr,caller)
   endif

   call MPI_Bcast(e_h%energy_hdm,1,mpi_real8,0,e_h%mpi_comm,ierr)

   call elsi_get_time(t1)

   write(info_str,"('  Finished density matrix correction')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

   call MPI_Barrier(e_h%mpi_comm,ierr)

   call elsi_check_mpi(e_h,"MPI_Barrier",ierr,caller)

end subroutine

!>
!! This routine computes the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_pexsi_real(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: edm(e_h%nnz_l_sp1)

   real(kind=r8)      :: mu_range
   real(kind=r8)      :: shift_width
   real(kind=r8)      :: local_energy
   real(kind=r8)      :: factor_min
   real(kind=r8)      :: factor_max
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   integer(kind=i4)   :: aux_min
   integer(kind=i4)   :: aux_max
   integer(kind=i4)   :: i
   integer(kind=i4)   :: ierr
   logical            :: converged
   character(len=200) :: info_str

   real(kind=r8), allocatable :: shifts(:)
   real(kind=r8), allocatable :: tmp_real(:)
   real(kind=r8), allocatable :: send_buf(:)

   character(len=40), parameter :: caller = "elsi_compute_edm_pexsi_real"

   call elsi_get_time(t0)

   call f_ppexsi_calculate_edm_correction_real(e_h%pexsi_plan,&
           e_h%pexsi_options,ierr)

   if(ierr /= 0) then
      call elsi_stop(" Energy density matrix correction failed.",e_h,caller)
   endif

   ! Get energy density matrix
   call elsi_allocate(e_h,tmp_real,e_h%nnz_l_sp1,"tmp_real",caller)

   call f_ppexsi_retrieve_real_edm(e_h%pexsi_plan,tmp_real,local_energy,ierr)

   if(ierr /= 0) then
      call elsi_stop(" Failed to get energy density matirx.",e_h,caller)
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
      if(e_h%ne_vec_pexsi(i) < e_h%n_electrons-&
         e_h%pexsi_options%numElectronPEXSITolerance) then
         e_h%pexsi_options%muMin0 = shifts(i)
         aux_min = i
      endif
   enddo

   do i = e_h%pexsi_options%nPoints,1,-1
      if(e_h%ne_vec_pexsi(i) > e_h%n_electrons+&
         e_h%pexsi_options%numElectronPEXSITolerance) then
         e_h%pexsi_options%muMax0 = shifts(i)
         aux_max = i
      endif
   enddo

   if(e_h%pexsi_options%nPoints == 1) then
      ! Scale energy density matrix
      tmp_real  = (e_h%n_electrons/e_h%pexsi_ne)*tmp_real
      converged = .true.
      e_h%mu    = shifts(1)
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
         if(abs(e_h%ne_vec_pexsi(i)-e_h%n_electrons) < &
            e_h%pexsi_options%numElectronPEXSITolerance) then
            e_h%mu    = shifts(i)
            converged = .true.

            call MPI_Bcast(tmp_real,e_h%nnz_l_sp1,mpi_real8,i-1,&
                    e_h%pexsi_comm_among_point,ierr)

            call elsi_check_mpi(e_h,"MPI_Bcast",ierr,caller)

            exit
         endif
      enddo
   endif

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Energy density matrix
      factor_min = (e_h%ne_vec_pexsi(aux_max)-e_h%n_electrons)/&
                      (e_h%ne_vec_pexsi(aux_max)-e_h%ne_vec_pexsi(aux_min))
      factor_max = (e_h%n_electrons-e_h%ne_vec_pexsi(aux_min))/&
                      (e_h%ne_vec_pexsi(aux_max)-e_h%ne_vec_pexsi(aux_min))

      call elsi_allocate(e_h,send_buf,e_h%nnz_l_sp1,"send_buf",caller)

      if(e_h%pexsi_my_point == aux_min-1) then
         send_buf = factor_min*tmp_real
      elseif(e_h%pexsi_my_point == aux_max-1) then
         send_buf = factor_max*tmp_real
      endif

      call MPI_Allreduce(send_buf,tmp_real,e_h%nnz_l_sp1,mpi_real8,mpi_sum,&
              e_h%pexsi_comm_among_point,ierr)

      call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(e_h,send_buf,"send_buf")
   endif

   if(e_h%pexsi_my_prow == 0) then
      edm = tmp_real
   endif

   call elsi_deallocate(e_h,tmp_real,"tmp_real")
   call elsi_deallocate(e_h,shifts,"shifts")

   call elsi_get_time(t1)

   write(info_str,"('  Finished energy density matrix calculation')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine interfaces to PEXSI.
!!
subroutine elsi_solve_evp_pexsi_cmplx(e_h,ham,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(inout) :: ham(e_h%nnz_l_sp1)
   complex(kind=r8),  intent(inout) :: ovlp(e_h%nnz_l_sp1)
   complex(kind=r8),  intent(inout) :: dm(e_h%nnz_l_sp1)

   real(kind=r8)      :: ne_drv
   real(kind=r8)      :: mu_range
   real(kind=r8)      :: shift_width
   real(kind=r8)      :: local_energy
   real(kind=r8)      :: factor_min
   real(kind=r8)      :: factor_max
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   complex(kind=r8)   :: local_cmplx
   integer(kind=i4)   :: n_iner_steps
   integer(kind=i4)   :: n_shift
   integer(kind=i4)   :: aux_min
   integer(kind=i4)   :: aux_max
   integer(kind=i4)   :: i
   integer(kind=i4)   :: idx
   integer(kind=i4)   :: ierr
   logical            :: converged
   character(len=200) :: info_str

   real(kind=r8),    allocatable :: shifts(:)
   real(kind=r8),    allocatable :: inertias(:)
   real(kind=r8),    allocatable :: ne_lower(:)
   real(kind=r8),    allocatable :: ne_upper(:)
   complex(kind=r8), allocatable :: tmp_cmplx(:)
   real(kind=r8),    allocatable :: send_buf(:)
   complex(kind=r8), allocatable :: send_buf_cmplx(:)

   complex(kind=r8), external :: zdotu

   character(len=40), parameter :: caller = "elsi_solve_evp_pexsi_cmplx"

   ! Load sparse matrices for PEXSI
   if(e_h%ovlp_is_unit) then
      call f_ppexsi_load_complex_hs_matrix(e_h%pexsi_plan,e_h%pexsi_options,&
              e_h%n_basis,e_h%nnz_g,e_h%nnz_l_sp1,e_h%n_lcol_sp1,&
              e_h%col_ptr_pexsi,e_h%row_ind_pexsi,ham,1,ovlp,ierr)
   else
      call f_ppexsi_load_complex_hs_matrix(e_h%pexsi_plan,e_h%pexsi_options,&
              e_h%n_basis,e_h%nnz_g,e_h%nnz_l_sp1,e_h%n_lcol_sp1,&
              e_h%col_ptr_pexsi,e_h%row_ind_pexsi,ham,0,ovlp,ierr)
   endif

   if(ierr /= 0) then
      call elsi_stop(" Failed to load matrices.",e_h,caller)
   endif

   call elsi_say(e_h,"  Starting PEXSI density matrix solver")

   ! Symbolic factorization
   if(e_h%n_elsi_calls == 1) then
      call elsi_get_time(t0)

      call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(&
              e_h%pexsi_plan,e_h%pexsi_options,ierr)

      call f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix(&
              e_h%pexsi_plan,e_h%pexsi_options,ovlp,ierr)

      call elsi_get_time(t1)

      write(info_str,"('  Finished symbolic factorization')")
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_say(e_h,info_str)
   endif

   if(ierr /= 0) then
      call elsi_stop(" Symbolic factorization failed.",e_h,caller)
   endif

   ! Inertia counting
   call elsi_get_time(t0)

   n_iner_steps = 0
   mu_range     = e_h%pexsi_options%muMax0-e_h%pexsi_options%muMin0
   n_shift      = max(10,e_h%n_procs/e_h%pexsi_np_per_pole)

   call elsi_allocate(e_h,shifts,n_shift,"shifts",caller)
   call elsi_allocate(e_h,inertias,n_shift,"inertias",caller)
   call elsi_allocate(e_h,ne_lower,n_shift,"ne_lower",caller)
   call elsi_allocate(e_h,ne_upper,n_shift,"ne_upper",caller)

   do while(n_iner_steps < 10 .and. &
            mu_range > e_h%pexsi_options%muInertiaTolerance)
      n_iner_steps = n_iner_steps+1

      shift_width = mu_range/(n_shift-1)

      ne_lower = 0.0_r8
      ne_upper = e_h%n_basis*e_h%spin_degen

      do i = 1,n_shift
         shifts(i) = e_h%pexsi_options%muMin0+(i-1)*shift_width
      enddo

      call f_ppexsi_inertia_count_complex_matrix(e_h%pexsi_plan,&
              e_h%pexsi_options,n_shift,shifts,inertias,ierr)

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
                 e_h%mpi_comm_all,ierr)

         call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

         call elsi_deallocate(e_h,send_buf,"send_buf")
      endif

      idx = ceiling(3*e_h%pexsi_options%temperature/shift_width)

      do i = idx+1,n_shift
         ne_lower(i)     = 0.5_r8*(inertias(i-idx)+inertias(i))
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

   call elsi_get_time(t1)

   if(n_iner_steps > 0) then
      write(info_str,"('  Finished inertia counting')")
      call elsi_say(e_h,info_str)
      write(info_str,"('  | Time :',F10.3,' s')") t1-t0
      call elsi_say(e_h,info_str)
   endif

   if(ierr /= 0) then
      call elsi_stop(" Inertia counting failed.",e_h,caller)
   endif

   ! Fermi operator expansion
   call elsi_get_time(t0)

   shift_width = mu_range/(e_h%pexsi_options%nPoints+1)

   call elsi_allocate(e_h,shifts,e_h%pexsi_options%nPoints,"shifts",caller)

   do i = 1,e_h%pexsi_options%nPoints
      shifts(i) = e_h%pexsi_options%muMin0+i*shift_width
   enddo

   do i = 1,e_h%pexsi_options%nPoints
      e_h%mu = shifts(i)

      if(e_h%pexsi_my_point == i-1) then
         call f_ppexsi_calculate_fermi_operator_complex(e_h%pexsi_plan,&
                 e_h%pexsi_options,e_h%mu,e_h%n_electrons,e_h%pexsi_ne,&
                 ne_drv,ierr)
      endif
   enddo

   call elsi_allocate(e_h,send_buf,e_h%pexsi_options%nPoints,"send_buf",caller)

   if(e_h%n_elsi_calls == 1) then
      call elsi_allocate(e_h,e_h%ne_vec_pexsi,e_h%pexsi_options%nPoints,&
              "ne_vec_pexsi",caller)
   endif

   send_buf(e_h%pexsi_my_point+1) = e_h%pexsi_ne*e_h%i_weight

   call MPI_Allreduce(send_buf,e_h%ne_vec_pexsi,e_h%pexsi_options%nPoints,&
           mpi_real8,mpi_sum,e_h%pexsi_comm_among_point,ierr)

   call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

   ! Get global number of electrons
   if(e_h%n_spins*e_h%n_kpts > 1) then
      if(e_h%myid == 0) then
         send_buf = e_h%ne_vec_pexsi
      else
         send_buf = 0.0_r8
      endif

      call MPI_Allreduce(send_buf,e_h%ne_vec_pexsi,e_h%pexsi_options%nPoints,&
              mpi_real8,mpi_sum,e_h%mpi_comm_all,ierr)

      call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)
   endif

   call elsi_deallocate(e_h,send_buf,"send_buf")

   call elsi_get_time(t1)

   write(info_str,"('  Finished Fermi operator calculation')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

   if(ierr /= 0) then
      call elsi_stop(" Fermi operator calculation failed.",e_h,caller)
   endif

   ! Get density matrix
   call elsi_get_time(t0)

   call elsi_allocate(e_h,tmp_cmplx,e_h%nnz_l_sp1,"tmp_cmplx",caller)

   call f_ppexsi_retrieve_complex_dm(e_h%pexsi_plan,tmp_cmplx,local_energy,ierr)

   if(ierr /= 0) then
      call elsi_stop(" Failed to get density matirx.",e_h,caller)
   endif

   ! Check convergence
   converged = .false.
   aux_min   = 0
   aux_max   = e_h%pexsi_options%nPoints+1

   do i = 1,e_h%pexsi_options%nPoints
      if(e_h%ne_vec_pexsi(i) < e_h%n_electrons-&
         e_h%pexsi_options%numElectronPEXSITolerance) then
         e_h%pexsi_options%muMin0 = shifts(i)
         aux_min = i
      endif
   enddo

   do i = e_h%pexsi_options%nPoints,1,-1
      if(e_h%ne_vec_pexsi(i) > e_h%n_electrons+&
         e_h%pexsi_options%numElectronPEXSITolerance) then
         e_h%pexsi_options%muMax0 = shifts(i)
         aux_max = i
      endif
   enddo

   if(e_h%pexsi_options%nPoints == 1) then
      ! Scale density matrix
      tmp_cmplx = (e_h%n_electrons/e_h%pexsi_ne)*tmp_cmplx
      converged = .true.
      e_h%mu    = shifts(1)
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
         if(abs(e_h%ne_vec_pexsi(i)-e_h%n_electrons) < &
            e_h%pexsi_options%numElectronPEXSITolerance) then
            e_h%mu    = shifts(i)
            converged = .true.

            call MPI_Bcast(tmp_cmplx,e_h%nnz_l_sp1,mpi_complex16,i-1,&
                    e_h%pexsi_comm_among_point,ierr)

            call elsi_check_mpi(e_h,"MPI_Bcast",ierr,caller)

            exit
         endif
      enddo
   endif

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Chemical potential
      e_h%mu = shifts(aux_min)+(e_h%n_electrons-e_h%ne_vec_pexsi(aux_min))/&
                  (e_h%ne_vec_pexsi(aux_max)-e_h%ne_vec_pexsi(aux_min))*&
                  (shifts(aux_max)-shifts(aux_min))

      ! Density matrix
      factor_min = (e_h%ne_vec_pexsi(aux_max)-e_h%n_electrons)/&
                      (e_h%ne_vec_pexsi(aux_max)-e_h%ne_vec_pexsi(aux_min))
      factor_max = (e_h%n_electrons-e_h%ne_vec_pexsi(aux_min))/&
                      (e_h%ne_vec_pexsi(aux_max)-e_h%ne_vec_pexsi(aux_min))

      call elsi_allocate(e_h,send_buf_cmplx,e_h%nnz_l_sp1,"send_buf_cmplx",&
              caller)

      if(e_h%pexsi_my_point == aux_min-1) then
         send_buf_cmplx = factor_min*tmp_cmplx
      elseif(e_h%pexsi_my_point == aux_max-1) then
         send_buf_cmplx = factor_max*tmp_cmplx
      endif

      call MPI_Allreduce(send_buf_cmplx,tmp_cmplx,e_h%nnz_l_sp1,mpi_complex16,&
              mpi_sum,e_h%pexsi_comm_among_point,ierr)

      call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(e_h,send_buf_cmplx,"send_buf_cmplx")
   endif

   if(e_h%pexsi_my_prow == 0) then
      dm = tmp_cmplx
   endif

   call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")
   call elsi_deallocate(e_h,shifts,"shifts")

   ! Compute energy = Tr(H*DM)
   if(e_h%pexsi_my_prow == 0) then
      local_cmplx  = zdotu(e_h%nnz_l_sp1,ham,1,dm,1)
      local_energy = real(local_cmplx,kind=r8)

      call MPI_Reduce(local_energy,e_h%energy_hdm,1,mpi_real8,mpi_sum,0,&
              e_h%pexsi_comm_in_pole,ierr)

      call elsi_check_mpi(e_h,"MPI_Reduce",ierr,caller)
   endif

   call MPI_Bcast(e_h%energy_hdm,1,mpi_real8,0,e_h%mpi_comm,ierr)

   call elsi_check_mpi(e_h,"MPI_Bcast",ierr,caller)

   call elsi_get_time(t1)

   write(info_str,"('  Finished density matrix correction')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

   call MPI_Barrier(e_h%mpi_comm,ierr)

   call elsi_check_mpi(e_h,"MPI_Barrier",ierr,caller)

end subroutine

!>
!! This routine computes the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_pexsi_cmplx(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(inout) :: edm(e_h%nnz_l_sp1)

   real(kind=r8)      :: mu_range
   real(kind=r8)      :: shift_width
   real(kind=r8)      :: local_energy
   real(kind=r8)      :: factor_min
   real(kind=r8)      :: factor_max
   real(kind=r8)      :: t0
   real(kind=r8)      :: t1
   integer(kind=i4)   :: aux_min
   integer(kind=i4)   :: aux_max
   integer(kind=i4)   :: i
   integer(kind=i4)   :: ierr
   logical            :: converged
   character(len=200) :: info_str

   real(kind=r8),    allocatable :: shifts(:)
   complex(kind=r8), allocatable :: tmp_cmplx(:)
   complex(kind=r8), allocatable :: send_buf_cmplx(:)

   character(len=40), parameter :: caller = "elsi_compute_edm_pexsi_cmplx"

   call elsi_get_time(t0)

   call f_ppexsi_calculate_edm_correction_complex(e_h%pexsi_plan,&
           e_h%pexsi_options,ierr)

   if(ierr /= 0) then
      call elsi_stop(" Energy density matrix correction failed.",e_h,caller)
   endif

   ! Get energy density matrix
   call elsi_allocate(e_h,tmp_cmplx,e_h%nnz_l_sp1,"tmp_cmplx",caller)

   call f_ppexsi_retrieve_complex_edm(e_h%pexsi_plan,tmp_cmplx,local_energy,&
           ierr)

   if(ierr /= 0) then
      call elsi_stop(" Failed to get energy density matirx.",e_h,caller)
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
      if(e_h%ne_vec_pexsi(i) < e_h%n_electrons-&
         e_h%pexsi_options%numElectronPEXSITolerance) then
         e_h%pexsi_options%muMin0 = shifts(i)
         aux_min = i
      endif
   enddo

   do i = e_h%pexsi_options%nPoints,1,-1
      if(e_h%ne_vec_pexsi(i) > e_h%n_electrons+&
         e_h%pexsi_options%numElectronPEXSITolerance) then
         e_h%pexsi_options%muMax0 = shifts(i)
         aux_max = i
      endif
   enddo

   if(e_h%pexsi_options%nPoints == 1) then
      ! Scale energy density matrix
      tmp_cmplx = (e_h%n_electrons/e_h%pexsi_ne)*tmp_cmplx
      converged = .true.
      e_h%mu    = shifts(1)
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
         if(abs(e_h%ne_vec_pexsi(i)-e_h%n_electrons) < &
            e_h%pexsi_options%numElectronPEXSITolerance) then
            e_h%mu    = shifts(i)
            converged = .true.

            call MPI_Bcast(tmp_cmplx,e_h%nnz_l_sp1,mpi_complex16,i-1,&
                    e_h%pexsi_comm_among_point,ierr)

            call elsi_check_mpi(e_h,"MPI_Bcast",ierr,caller)

            exit
         endif
      enddo
   endif

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Energy density matrix
      factor_min = (e_h%ne_vec_pexsi(aux_max)-e_h%n_electrons)/&
                      (e_h%ne_vec_pexsi(aux_max)-e_h%ne_vec_pexsi(aux_min))
      factor_max = (e_h%n_electrons-e_h%ne_vec_pexsi(aux_min))/&
                      (e_h%ne_vec_pexsi(aux_max)-e_h%ne_vec_pexsi(aux_min))

      call elsi_allocate(e_h,send_buf_cmplx,e_h%nnz_l_sp1,"send_buf_cmplx",&
              caller)

      if(e_h%pexsi_my_point == aux_min-1) then
         send_buf_cmplx = factor_min*tmp_cmplx
      elseif(e_h%pexsi_my_point == aux_max-1) then
         send_buf_cmplx = factor_max*tmp_cmplx
      endif

      call MPI_Allreduce(send_buf_cmplx,tmp_cmplx,e_h%nnz_l_sp1,mpi_complex16,&
              mpi_sum,e_h%pexsi_comm_among_point,ierr)

      call elsi_check_mpi(e_h,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(e_h,send_buf_cmplx,"send_buf_cmplx")
   endif

   if(e_h%pexsi_my_prow == 0) then
      edm = tmp_cmplx
   endif

   call elsi_deallocate(e_h,tmp_cmplx,"tmp_cmplx")
   call elsi_deallocate(e_h,shifts,"shifts")

   call elsi_get_time(t1)

   write(info_str,"('  Finished energy density matrix calculation')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine sets default PEXSI parameters.
!!
subroutine elsi_set_pexsi_default(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h

   character(len=40), parameter :: caller = "elsi_set_pexsi_default"

   if(e_h%handle_ready) then
      e_h%handle_changed = .true.
   endif

   ! Use the PEXSI Default options
   call f_ppexsi_set_default_options(e_h%pexsi_options)

   ! Use 1 process in symbolic factorization
   e_h%pexsi_options%npSymbFact = 1

   ! Use parallel matrix ordering by ParMETIS/PtScotch
   ! Note: must use serial ordering on some platform (segfault otherwise)
   e_h%pexsi_options%ordering = 0

end subroutine

end module ELSI_PEXSI
