! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides interfaces to PEXSI.
!!
module ELSI_PEXSI

   use ELSI_CONSTANTS,     only: UNSET,PEXSI_CSC
   use ELSI_DATATYPE,      only: elsi_param_t,elsi_basic_t
   use ELSI_IO,            only: elsi_say,elsi_get_time
   use ELSI_MALLOC,        only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,           only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8,&
                                 mpi_complex16
   use ELSI_PRECISION,     only: r8,i4
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

   public :: elsi_init_pexsi
   public :: elsi_set_pexsi_default
   public :: elsi_cleanup_pexsi
   public :: elsi_solve_pexsi
   public :: elsi_compute_edm_pexsi

   interface elsi_solve_pexsi
      module procedure elsi_solve_pexsi_real
      module procedure elsi_solve_pexsi_cmplx
   end interface

   interface elsi_compute_edm_pexsi
      module procedure elsi_compute_edm_pexsi_real
      module procedure elsi_compute_edm_pexsi_cmplx
   end interface

contains

!>
!! This routine initializes PEXSI and its processor grid.
!!
subroutine elsi_init_pexsi(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh

   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: log_id
   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_init_pexsi"

   if(.not. ph%pexsi_started) then
      ph%pexsi_options%spin = ph%spin_degen

      ! Point parallelization
      ph%pexsi_np_per_point = bh%n_procs/ph%pexsi_options%nPoints
      ph%pexsi_my_point     = bh%myid/ph%pexsi_np_per_point
      ph%pexsi_myid_point   = mod(bh%myid,ph%pexsi_np_per_point)

      if(ph%pexsi_np_per_pole == UNSET) then
         j = nint(sqrt(real(ph%pexsi_np_per_point,kind=r8)))

         do i = ph%pexsi_np_per_point,j,-1
            if(mod(ph%pexsi_np_per_point,i) == 0) then
               if((ph%pexsi_np_per_point/i) > ph%pexsi_options%numPole) then
                  exit
               endif

               ph%pexsi_np_per_pole = i
            endif
         enddo
      endif

      ! Set square-like process grid for selected inversion of each pole
      do i = nint(sqrt(real(ph%pexsi_np_per_pole,kind=r8))),2,-1
         if(mod(ph%pexsi_np_per_pole,i) == 0) then
            exit
         endif
      enddo

      ! PEXSI process grid
      ph%pexsi_n_prow  = i
      ph%pexsi_n_pcol  = ph%pexsi_np_per_pole/i
      ph%pexsi_my_prow = bh%myid/ph%pexsi_np_per_pole
      ph%pexsi_my_pcol = mod(bh%myid,ph%pexsi_np_per_pole)

      ! PEXSI MPI communicators
      call MPI_Comm_split(bh%comm,ph%pexsi_my_prow,ph%pexsi_my_pcol,&
              ph%pexsi_comm_intra_pole,ierr)

      call elsi_check_mpi(bh,"MPI_Comm_split",ierr,caller)

      call MPI_Comm_split(bh%comm,ph%pexsi_my_pcol,ph%pexsi_my_prow,&
              ph%pexsi_comm_inter_pole,ierr)

      call elsi_check_mpi(bh,"MPI_Comm_split",ierr,caller)

      call MPI_Comm_split(bh%comm,ph%pexsi_myid_point,ph%pexsi_my_point,&
              ph%pexsi_comm_inter_point,ierr)

      call elsi_check_mpi(bh,"MPI_Comm_split",ierr,caller)

      if(.not. bh%pexsi_csc_ready) then
         ! Set up 1D block distribution
         bh%n_lcol_sp1 = ph%n_basis/ph%pexsi_np_per_pole

         ! The last process holds all remaining columns
         if(ph%pexsi_my_pcol == ph%pexsi_np_per_pole-1) then
            bh%n_lcol_sp1 = ph%n_basis-(ph%pexsi_np_per_pole-1)*bh%n_lcol_sp1
         endif
      endif

      ! Only one process outputs
      if(bh%myid_all == 0) then
         log_id = 0
      else
         log_id = -1
      endif

      ph%pexsi_plan = f_ppexsi_plan_initialize(bh%comm,ph%pexsi_n_prow,&
                         ph%pexsi_n_pcol,log_id,ierr)

      if(ierr /= 0) then
         call elsi_stop(bh,"Initialization failed.",caller)
      endif

      if(bh%n_lcol_sp == UNSET) then
         bh%n_lcol_sp = bh%n_lcol_sp1
      endif

      ph%pexsi_started = .true.
   endif

end subroutine

!>
!! This routine interfaces to PEXSI.
!!
subroutine elsi_solve_pexsi_real(ph,bh,row_ind,col_ptr,ne_vec,ham,ovlp,dm)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8),      intent(out)   :: ne_vec(ph%pexsi_options%nPoints)
   real(kind=r8),      intent(in)    :: ham(bh%nnz_l_sp1)
   real(kind=r8),      intent(in)    :: ovlp(bh%nnz_l_sp1)
   real(kind=r8),      intent(out)   :: dm(bh%nnz_l_sp1)

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

   character(len=*), parameter :: caller = "elsi_solve_pexsi_real"

   ! Load sparse matrices for PEXSI
   if(ph%ovlp_is_unit) then
      call f_ppexsi_load_real_hs_matrix(ph%pexsi_plan,ph%pexsi_options,&
              ph%n_basis,bh%nnz_g,bh%nnz_l_sp1,bh%n_lcol_sp1,col_ptr,row_ind,&
              ham,1,ovlp,ierr)
   else
      call f_ppexsi_load_real_hs_matrix(ph%pexsi_plan,ph%pexsi_options,&
              ph%n_basis,bh%nnz_g,bh%nnz_l_sp1,bh%n_lcol_sp1,col_ptr,row_ind,&
              ham,0,ovlp,ierr)
   endif

   if(ierr /= 0) then
      call elsi_stop(bh,"Failed to load matrices.",caller)
   endif

   write(info_str,"(2X,A)") "Starting PEXSI density matrix solver"
   call elsi_say(bh,info_str)

   ! Symbolic factorization
   if(ph%n_calls == 1) then
      call elsi_get_time(t0)

      call f_ppexsi_symbolic_factorize_real_symmetric_matrix(ph%pexsi_plan,&
              ph%pexsi_options,ierr)

      call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(ph%pexsi_plan,&
              ph%pexsi_options,ierr)

      call elsi_get_time(t1)

      write(info_str,"(2X,A)") "Finished symbolic factorization"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,info_str)
   endif

   if(ierr /= 0) then
      call elsi_stop(bh,"Symbolic factorization failed.",caller)
   endif

   ! Inertia counting
   mu_range = ph%pexsi_options%muMax0-ph%pexsi_options%muMin0

   if(mu_range > ph%pexsi_options%muInertiaTolerance) then
      call elsi_get_time(t0)

      n_iner_steps = 0
      n_shift      = max(10,bh%n_procs/ph%pexsi_np_per_pole)

      call elsi_allocate(bh,shifts,n_shift,"shifts",caller)
      call elsi_allocate(bh,inertias,n_shift,"inertias",caller)
      call elsi_allocate(bh,ne_lower,n_shift,"ne_lower",caller)
      call elsi_allocate(bh,ne_upper,n_shift,"ne_upper",caller)

      do while(n_iner_steps < 10 .and.&
               mu_range > ph%pexsi_options%muInertiaTolerance)
         n_iner_steps = n_iner_steps+1
         shift_width  = mu_range/(n_shift-1)
         ne_lower     = 0.0_r8
         ne_upper     = ph%n_basis*ph%spin_degen

         do i = 1,n_shift
            shifts(i) = ph%pexsi_options%muMin0+(i-1)*shift_width
         enddo

         call f_ppexsi_inertia_count_real_matrix(ph%pexsi_plan,&
                 ph%pexsi_options,n_shift,shifts,inertias,ierr)

         inertias = inertias*ph%spin_degen*ph%i_weight

         ! Get global inertias
         if(ph%n_spins*ph%n_kpts > 1) then
            call elsi_allocate(bh,send_buf,n_shift,"send_buf",caller)

            if(bh%myid == 0) then
               send_buf = inertias
            else
               send_buf = 0.0_r8
            endif

            call MPI_Allreduce(send_buf,inertias,n_shift,mpi_real8,mpi_sum,&
                    bh%comm_all,ierr)

            call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

            call elsi_deallocate(bh,send_buf,"send_buf")
         endif

         idx = ceiling(3*ph%pexsi_options%temperature/shift_width)

         do i = idx+1,n_shift
            ne_lower(i)     = 0.5_r8*(inertias(i-idx)+inertias(i))
            ne_upper(i-idx) = ne_lower(i)
         enddo

         aux_min = 1
         aux_max = n_shift

         do i = 2,n_shift-1
            if(ne_upper(i) < ph%n_electrons .and.&
               ne_upper(i+1) >= ph%n_electrons) then
               aux_min = i
            endif

            if(ne_lower(i) > ph%n_electrons .and.&
               ne_lower(i-1) <= ph%n_electrons) then
               aux_max = i
            endif
         enddo

         if(aux_min == 1 .and. aux_max == n_shift) then
            exit
         else
            ph%pexsi_options%muMin0 = shifts(aux_min)
            ph%pexsi_options%muMax0 = shifts(aux_max)

            mu_range = ph%pexsi_options%muMax0-ph%pexsi_options%muMin0
         endif
      enddo

      call elsi_deallocate(bh,shifts,"shifts")
      call elsi_deallocate(bh,inertias,"inertias")
      call elsi_deallocate(bh,ne_lower,"ne_lower")
      call elsi_deallocate(bh,ne_upper,"ne_upper")

      call elsi_get_time(t1)

      write(info_str,"(2X,A)") "Finished inertia counting"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,info_str)

      if(ierr /= 0) then
         call elsi_stop(bh,"Inertia counting failed.",caller)
      endif
   endif

   ! Fermi operator expansion
   call elsi_get_time(t0)

   shift_width = mu_range/(ph%pexsi_options%nPoints+1)

   call elsi_allocate(bh,shifts,ph%pexsi_options%nPoints,"shifts",caller)

   do i = 1,ph%pexsi_options%nPoints
      shifts(i) = ph%pexsi_options%muMin0+i*shift_width
   enddo

   do i = 1,ph%pexsi_options%nPoints
      ph%mu = shifts(i)

      if(ph%pexsi_my_point == i-1) then
         call f_ppexsi_calculate_fermi_operator_real3(ph%pexsi_plan,&
                 ph%pexsi_options,ph%mu,ph%n_electrons,ph%pexsi_ne,ierr)
      endif
   enddo

   call elsi_allocate(bh,send_buf,ph%pexsi_options%nPoints,"send_buf",caller)

   send_buf(ph%pexsi_my_point+1) = ph%pexsi_ne*ph%i_weight

   call MPI_Allreduce(send_buf,ne_vec,ph%pexsi_options%nPoints,mpi_real8,&
           mpi_sum,ph%pexsi_comm_inter_point,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   ! Get global number of electrons
   if(ph%n_spins*ph%n_kpts > 1) then
      if(bh%myid == 0) then
         send_buf = ne_vec
      else
         send_buf = 0.0_r8
      endif

      call MPI_Allreduce(send_buf,ne_vec,ph%pexsi_options%nPoints,mpi_real8,&
              mpi_sum,bh%comm_all,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   endif

   call elsi_deallocate(bh,send_buf,"send_buf")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished Fermi operator calculation"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

   if(ierr /= 0) then
      call elsi_stop(bh,"Fermi operator calculation failed.",caller)
   endif

   ! Get density matrix
   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp_real,bh%nnz_l_sp1,"tmp_real",caller)

   call f_ppexsi_retrieve_real_dm(ph%pexsi_plan,tmp_real,local_energy,ierr)

   if(ierr /= 0) then
      call elsi_stop(bh,"Failed to get density matirx.",caller)
   endif

   ! Check convergence
   converged = .false.
   aux_min   = 0
   aux_max   = ph%pexsi_options%nPoints+1

   do i = 1,ph%pexsi_options%nPoints
      if(ne_vec(i) <&
         ph%n_electrons-ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMin0 = shifts(i)
         aux_min                 = i
      endif
   enddo

   do i = ph%pexsi_options%nPoints,1,-1
      if(ne_vec(i) >&
         ph%n_electrons+ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMax0 = shifts(i)
         aux_max                 = i
      endif
   enddo

   if(ph%pexsi_options%nPoints == 1) then
      ! Scale density matrix
      tmp_real  = (ph%n_electrons/ph%pexsi_ne)*tmp_real
      converged = .true.
      ph%mu     = shifts(1)
   else
      ! Safety check
      if(aux_min == 0) then
         aux_min = 1

         if(aux_max <= aux_min) then
            aux_max = 2
         endif
      endif

      if(aux_max == ph%pexsi_options%nPoints+1) then
         aux_max = ph%pexsi_options%nPoints

         if(aux_min >= aux_max) then
            aux_min = ph%pexsi_options%nPoints-1
         endif
      endif

      do i = aux_min,aux_max
         if(abs(ne_vec(i)-ph%n_electrons) <&
            ph%pexsi_options%numElectronPEXSITolerance) then
            ph%mu     = shifts(i)
            converged = .true.

            call MPI_Bcast(tmp_real,bh%nnz_l_sp1,mpi_real8,i-1,&
                    ph%pexsi_comm_inter_point,ierr)

            call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

            exit
         endif
      enddo
   endif

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Chemical potential
      ph%mu = shifts(aux_min)+(ph%n_electrons-ne_vec(aux_min))/&
                 (ne_vec(aux_max)-ne_vec(aux_min))*&
                 (shifts(aux_max)-shifts(aux_min))

      ! Density matrix
      factor_min = (ne_vec(aux_max)-ph%n_electrons)/&
                      (ne_vec(aux_max)-ne_vec(aux_min))
      factor_max = (ph%n_electrons-ne_vec(aux_min))/&
                      (ne_vec(aux_max)-ne_vec(aux_min))

      call elsi_allocate(bh,send_buf,bh%nnz_l_sp1,"send_buf",caller)

      if(ph%pexsi_my_point == aux_min-1) then
         send_buf = factor_min*tmp_real
      elseif(ph%pexsi_my_point == aux_max-1) then
         send_buf = factor_max*tmp_real
      endif

      call MPI_Allreduce(send_buf,tmp_real,bh%nnz_l_sp1,mpi_real8,mpi_sum,&
              ph%pexsi_comm_inter_point,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(bh,send_buf,"send_buf")
   endif

   if(.not. (ph%matrix_format == PEXSI_CSC .and. ph%pexsi_my_prow /= 0)) then
      dm = tmp_real
   endif

   call elsi_deallocate(bh,tmp_real,"tmp_real")
   call elsi_deallocate(bh,shifts,"shifts")

   ! Compute energy = Tr(H*DM)
   if(ph%pexsi_my_prow == 0) then
      local_energy = ddot(bh%nnz_l_sp1,ham,1,dm,1)

      call MPI_Reduce(local_energy,ph%ebs,1,mpi_real8,mpi_sum,0,&
              ph%pexsi_comm_intra_pole,ierr)

      call elsi_check_mpi(bh,"MPI_Reduce",ierr,caller)
   endif

   call MPI_Bcast(ph%ebs,1,mpi_real8,0,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished density matrix correction"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine computes the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_pexsi_real(ph,bh,ne_vec,edm)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   real(kind=r8),      intent(in)    :: ne_vec(ph%pexsi_options%nPoints)
   real(kind=r8),      intent(out)   :: edm(bh%nnz_l_sp1)

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

   character(len=*), parameter :: caller = "elsi_compute_edm_pexsi_real"

   call elsi_get_time(t0)

   call f_ppexsi_calculate_edm_correction_real(ph%pexsi_plan,ph%pexsi_options,&
           ierr)

   if(ierr /= 0) then
      call elsi_stop(bh,"Energy density matrix correction failed.",caller)
   endif

   ! Get energy density matrix
   call elsi_allocate(bh,tmp_real,bh%nnz_l_sp1,"tmp_real",caller)

   call f_ppexsi_retrieve_real_edm(ph%pexsi_plan,tmp_real,local_energy,ierr)

   if(ierr /= 0) then
      call elsi_stop(bh,"Failed to get energy density matirx.",caller)
   endif

   ! Check convergence
   mu_range    = ph%pexsi_options%muMax0-ph%pexsi_options%muMin0
   shift_width = mu_range/(ph%pexsi_options%nPoints+1)
   converged   = .false.
   aux_min     = 0
   aux_max     = ph%pexsi_options%nPoints+1

   call elsi_allocate(bh,shifts,ph%pexsi_options%nPoints,"shifts",caller)

   do i = 1,ph%pexsi_options%nPoints
      shifts(i) = ph%pexsi_options%muMin0+i*shift_width
   enddo

   do i = 1,ph%pexsi_options%nPoints
      if(ne_vec(i) <&
         ph%n_electrons-ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMin0 = shifts(i)
         aux_min                 = i
      endif
   enddo

   do i = ph%pexsi_options%nPoints,1,-1
      if(ne_vec(i) >&
         ph%n_electrons+ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMax0 = shifts(i)
         aux_max                 = i
      endif
   enddo

   if(ph%pexsi_options%nPoints == 1) then
      ! Scale energy density matrix
      tmp_real  = (ph%n_electrons/ph%pexsi_ne)*tmp_real
      converged = .true.
      ph%mu     = shifts(1)
   else
      ! Safety check
      if(aux_min == 0) then
         aux_min = 1

         if(aux_max <= aux_min) then
            aux_max = 2
         endif
      endif

      if(aux_max == ph%pexsi_options%nPoints+1) then
         aux_max = ph%pexsi_options%nPoints

         if(aux_min >= aux_max) then
            aux_min = ph%pexsi_options%nPoints-1
         endif
      endif

      do i = aux_min,aux_max
         if(abs(ne_vec(i)-ph%n_electrons) <&
            ph%pexsi_options%numElectronPEXSITolerance) then
            ph%mu     = shifts(i)
            converged = .true.

            call MPI_Bcast(tmp_real,bh%nnz_l_sp1,mpi_real8,i-1,&
                    ph%pexsi_comm_inter_point,ierr)

            call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

            exit
         endif
      enddo
   endif

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Energy density matrix
      factor_min = (ne_vec(aux_max)-ph%n_electrons)/&
                      (ne_vec(aux_max)-ne_vec(aux_min))
      factor_max = (ph%n_electrons-ne_vec(aux_min))/&
                      (ne_vec(aux_max)-ne_vec(aux_min))

      call elsi_allocate(bh,send_buf,bh%nnz_l_sp1,"send_buf",caller)

      if(ph%pexsi_my_point == aux_min-1) then
         send_buf = factor_min*tmp_real
      elseif(ph%pexsi_my_point == aux_max-1) then
         send_buf = factor_max*tmp_real
      endif

      call MPI_Allreduce(send_buf,tmp_real,bh%nnz_l_sp1,mpi_real8,mpi_sum,&
              ph%pexsi_comm_inter_point,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(bh,send_buf,"send_buf")
   endif

   if(.not. (ph%matrix_format == PEXSI_CSC .and. ph%pexsi_my_prow /= 0)) then
      edm = tmp_real
   endif

   call elsi_deallocate(bh,tmp_real,"tmp_real")
   call elsi_deallocate(bh,shifts,"shifts")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished energy density matrix calculation"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine interfaces to PEXSI.
!!
subroutine elsi_solve_pexsi_cmplx(ph,bh,row_ind,col_ptr,ne_vec,ham,ovlp,dm)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   integer(kind=i4),   intent(in)    :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4),   intent(in)    :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8),      intent(out)   :: ne_vec(ph%pexsi_options%nPoints)
   complex(kind=r8),   intent(in)    :: ham(bh%nnz_l_sp1)
   complex(kind=r8),   intent(inout) :: ovlp(bh%nnz_l_sp1)
   complex(kind=r8),   intent(out)   :: dm(bh%nnz_l_sp1)

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

   complex(kind=r8), external :: zdotc

   character(len=*), parameter :: caller = "elsi_solve_pexsi_cmplx"

   ! Load sparse matrices for PEXSI
   if(ph%ovlp_is_unit) then
      call f_ppexsi_load_complex_hs_matrix(ph%pexsi_plan,ph%pexsi_options,&
              ph%n_basis,bh%nnz_g,bh%nnz_l_sp1,bh%n_lcol_sp1,col_ptr,row_ind,&
              ham,1,ovlp,ierr)
   else
      call f_ppexsi_load_complex_hs_matrix(ph%pexsi_plan,ph%pexsi_options,&
              ph%n_basis,bh%nnz_g,bh%nnz_l_sp1,bh%n_lcol_sp1,col_ptr,row_ind,&
              ham,0,ovlp,ierr)
   endif

   if(ierr /= 0) then
      call elsi_stop(bh,"Failed to load matrices.",caller)
   endif

   write(info_str,"(2X,A)") "Starting PEXSI density matrix solver"
   call elsi_say(bh,info_str)

   ! Symbolic factorization
   if(ph%n_calls == 1) then
      call elsi_get_time(t0)

      call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(ph%pexsi_plan,&
              ph%pexsi_options,ierr)

      call f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix(&
              ph%pexsi_plan,ph%pexsi_options,ovlp,ierr)

      call elsi_get_time(t1)

      write(info_str,"(2X,A)") "Finished symbolic factorization"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,info_str)
   endif

   if(ierr /= 0) then
      call elsi_stop(bh,"Symbolic factorization failed.",caller)
   endif

   ! Inertia counting
   mu_range = ph%pexsi_options%muMax0-ph%pexsi_options%muMin0

   if(mu_range > ph%pexsi_options%muInertiaTolerance) then
      call elsi_get_time(t0)

      n_iner_steps = 0
      n_shift      = max(10,bh%n_procs/ph%pexsi_np_per_pole)

      call elsi_allocate(bh,shifts,n_shift,"shifts",caller)
      call elsi_allocate(bh,inertias,n_shift,"inertias",caller)
      call elsi_allocate(bh,ne_lower,n_shift,"ne_lower",caller)
      call elsi_allocate(bh,ne_upper,n_shift,"ne_upper",caller)

      do while(n_iner_steps < 10 .and.&
               mu_range > ph%pexsi_options%muInertiaTolerance)
         n_iner_steps = n_iner_steps+1
         shift_width  = mu_range/(n_shift-1)
         ne_lower     = 0.0_r8
         ne_upper     = ph%n_basis*ph%spin_degen

         do i = 1,n_shift
            shifts(i) = ph%pexsi_options%muMin0+(i-1)*shift_width
         enddo

         call f_ppexsi_inertia_count_complex_matrix(ph%pexsi_plan,&
                 ph%pexsi_options,n_shift,shifts,inertias,ierr)

         inertias = inertias*ph%spin_degen*ph%i_weight

         ! Get global inertias
         if(ph%n_spins*ph%n_kpts > 1) then
            call elsi_allocate(bh,send_buf,n_shift,"send_buf",caller)

            if(bh%myid == 0) then
               send_buf = inertias
            else
               send_buf = 0.0_r8
            endif

            call MPI_Allreduce(send_buf,inertias,n_shift,mpi_real8,mpi_sum,&
                    bh%comm_all,ierr)

            call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

            call elsi_deallocate(bh,send_buf,"send_buf")
         endif

         idx = ceiling(3*ph%pexsi_options%temperature/shift_width)

         do i = idx+1,n_shift
            ne_lower(i)     = 0.5_r8*(inertias(i-idx)+inertias(i))
            ne_upper(i-idx) = ne_lower(i)
         enddo

         aux_min = 1
         aux_max = n_shift

         do i = 2,n_shift-1
            if(ne_upper(i) < ph%n_electrons .and.&
               ne_upper(i+1) >= ph%n_electrons) then
               aux_min = i
            endif

            if(ne_lower(i) > ph%n_electrons .and.&
               ne_lower(i-1) <= ph%n_electrons) then
               aux_max = i
            endif
         enddo

         if(aux_min == 1 .and. aux_max == n_shift) then
            exit
         else
            ph%pexsi_options%muMin0 = shifts(aux_min)
            ph%pexsi_options%muMax0 = shifts(aux_max)

            mu_range = ph%pexsi_options%muMax0-ph%pexsi_options%muMin0
         endif
      enddo

      call elsi_deallocate(bh,shifts,"shifts")
      call elsi_deallocate(bh,inertias,"inertias")
      call elsi_deallocate(bh,ne_lower,"ne_lower")
      call elsi_deallocate(bh,ne_upper,"ne_upper")

      call elsi_get_time(t1)

      write(info_str,"(2X,A)") "Finished inertia counting"
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,info_str)

      if(ierr /= 0) then
         call elsi_stop(bh,"Inertia counting failed.",caller)
      endif
   endif

   ! Fermi operator expansion
   call elsi_get_time(t0)

   shift_width = mu_range/(ph%pexsi_options%nPoints+1)

   call elsi_allocate(bh,shifts,ph%pexsi_options%nPoints,"shifts",caller)

   do i = 1,ph%pexsi_options%nPoints
      shifts(i) = ph%pexsi_options%muMin0+i*shift_width
   enddo

   do i = 1,ph%pexsi_options%nPoints
      ph%mu = shifts(i)

      if(ph%pexsi_my_point == i-1) then
         call f_ppexsi_calculate_fermi_operator_complex(ph%pexsi_plan,&
                 ph%pexsi_options,ph%mu,ph%n_electrons,ph%pexsi_ne,ne_drv,ierr)
      endif
   enddo

   call elsi_allocate(bh,send_buf,ph%pexsi_options%nPoints,"send_buf",caller)

   send_buf(ph%pexsi_my_point+1) = ph%pexsi_ne*ph%i_weight

   call MPI_Allreduce(send_buf,ne_vec,ph%pexsi_options%nPoints,mpi_real8,&
           mpi_sum,ph%pexsi_comm_inter_point,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   ! Get global number of electrons
   if(ph%n_spins*ph%n_kpts > 1) then
      if(bh%myid == 0) then
         send_buf = ne_vec
      else
         send_buf = 0.0_r8
      endif

      call MPI_Allreduce(send_buf,ne_vec,ph%pexsi_options%nPoints,mpi_real8,&
              mpi_sum,bh%comm_all,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   endif

   call elsi_deallocate(bh,send_buf,"send_buf")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished Fermi operator calculation"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

   if(ierr /= 0) then
      call elsi_stop(bh,"Fermi operator calculation failed.",caller)
   endif

   ! Get density matrix
   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp_cmplx,bh%nnz_l_sp1,"tmp_cmplx",caller)

   call f_ppexsi_retrieve_complex_dm(ph%pexsi_plan,tmp_cmplx,local_energy,ierr)

   if(ierr /= 0) then
      call elsi_stop(bh,"Failed to get density matirx.",caller)
   endif

   ! Check convergence
   converged = .false.
   aux_min   = 0
   aux_max   = ph%pexsi_options%nPoints+1

   do i = 1,ph%pexsi_options%nPoints
      if(ne_vec(i) <&
         ph%n_electrons-ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMin0 = shifts(i)
         aux_min                 = i
      endif
   enddo

   do i = ph%pexsi_options%nPoints,1,-1
      if(ne_vec(i) >&
         ph%n_electrons+ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMax0 = shifts(i)
         aux_max                 = i
      endif
   enddo

   if(ph%pexsi_options%nPoints == 1) then
      ! Scale density matrix
      tmp_cmplx = (ph%n_electrons/ph%pexsi_ne)*tmp_cmplx
      converged = .true.
      ph%mu     = shifts(1)
   else
      ! Safety check
      if(aux_min == 0) then
         aux_min = 1

         if(aux_max <= aux_min) then
            aux_max = 2
         endif
      endif

      if(aux_max == ph%pexsi_options%nPoints+1) then
         aux_max = ph%pexsi_options%nPoints

         if(aux_min >= aux_max) then
            aux_min = ph%pexsi_options%nPoints-1
         endif
      endif

      do i = aux_min,aux_max
         if(abs(ne_vec(i)-ph%n_electrons) <&
            ph%pexsi_options%numElectronPEXSITolerance) then
            ph%mu     = shifts(i)
            converged = .true.

            call MPI_Bcast(tmp_cmplx,bh%nnz_l_sp1,mpi_complex16,i-1,&
                    ph%pexsi_comm_inter_point,ierr)

            call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

            exit
         endif
      enddo
   endif

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Chemical potential
      ph%mu = shifts(aux_min)+(ph%n_electrons-ne_vec(aux_min))/&
                 (ne_vec(aux_max)-ne_vec(aux_min))*&
                 (shifts(aux_max)-shifts(aux_min))

      ! Density matrix
      factor_min = (ne_vec(aux_max)-ph%n_electrons)/&
                      (ne_vec(aux_max)-ne_vec(aux_min))
      factor_max = (ph%n_electrons-ne_vec(aux_min))/&
                      (ne_vec(aux_max)-ne_vec(aux_min))

      call elsi_allocate(bh,send_buf_cmplx,bh%nnz_l_sp1,"send_buf_cmplx",caller)

      if(ph%pexsi_my_point == aux_min-1) then
         send_buf_cmplx = factor_min*tmp_cmplx
      elseif(ph%pexsi_my_point == aux_max-1) then
         send_buf_cmplx = factor_max*tmp_cmplx
      endif

      call MPI_Allreduce(send_buf_cmplx,tmp_cmplx,bh%nnz_l_sp1,mpi_complex16,&
              mpi_sum,ph%pexsi_comm_inter_point,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(bh,send_buf_cmplx,"send_buf_cmplx")
   endif

   if(.not. (ph%matrix_format == PEXSI_CSC .and. ph%pexsi_my_prow /= 0)) then
      dm = tmp_cmplx
   endif

   call elsi_deallocate(bh,tmp_cmplx,"tmp_cmplx")
   call elsi_deallocate(bh,shifts,"shifts")

   ! Compute energy = Tr(H*DM)
   if(ph%pexsi_my_prow == 0) then
      local_cmplx  = zdotc(bh%nnz_l_sp1,ham,1,dm,1)
      local_energy = real(local_cmplx,kind=r8)

      call MPI_Reduce(local_energy,ph%ebs,1,mpi_real8,mpi_sum,0,&
              ph%pexsi_comm_intra_pole,ierr)

      call elsi_check_mpi(bh,"MPI_Reduce",ierr,caller)
   endif

   call MPI_Bcast(ph%ebs,1,mpi_real8,0,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished density matrix correction"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine computes the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_pexsi_cmplx(ph,bh,ne_vec,edm)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in)    :: bh
   real(kind=r8),      intent(in)    :: ne_vec(ph%pexsi_options%nPoints)
   complex(kind=r8),   intent(out)   :: edm(bh%nnz_l_sp1)

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

   character(len=*), parameter :: caller = "elsi_compute_edm_pexsi_cmplx"

   call elsi_get_time(t0)

   call f_ppexsi_calculate_edm_correction_complex(ph%pexsi_plan,&
           ph%pexsi_options,ierr)

   if(ierr /= 0) then
      call elsi_stop(bh,"Energy density matrix correction failed.",caller)
   endif

   ! Get energy density matrix
   call elsi_allocate(bh,tmp_cmplx,bh%nnz_l_sp1,"tmp_cmplx",caller)

   call f_ppexsi_retrieve_complex_edm(ph%pexsi_plan,tmp_cmplx,local_energy,ierr)

   if(ierr /= 0) then
      call elsi_stop(bh,"Failed to get energy density matirx.",caller)
   endif

   ! Check convergence
   mu_range    = ph%pexsi_options%muMax0-ph%pexsi_options%muMin0
   shift_width = mu_range/(ph%pexsi_options%nPoints+1)
   converged   = .false.
   aux_min     = 0
   aux_max     = ph%pexsi_options%nPoints+1

   call elsi_allocate(bh,shifts,ph%pexsi_options%nPoints,"shifts",caller)

   do i = 1,ph%pexsi_options%nPoints
      shifts(i) = ph%pexsi_options%muMin0+i*shift_width
   enddo

   do i = 1,ph%pexsi_options%nPoints
      if(ne_vec(i) <&
         ph%n_electrons-ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMin0 = shifts(i)
         aux_min                 = i
      endif
   enddo

   do i = ph%pexsi_options%nPoints,1,-1
      if(ne_vec(i) >&
         ph%n_electrons+ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMax0 = shifts(i)
         aux_max                 = i
      endif
   enddo

   if(ph%pexsi_options%nPoints == 1) then
      ! Scale energy density matrix
      tmp_cmplx = (ph%n_electrons/ph%pexsi_ne)*tmp_cmplx
      converged = .true.
      ph%mu     = shifts(1)
   else
      ! Safety check
      if(aux_min == 0) then
         aux_min = 1

         if(aux_max <= aux_min) then
            aux_max = 2
         endif
      endif

      if(aux_max == ph%pexsi_options%nPoints+1) then
         aux_max = ph%pexsi_options%nPoints

         if(aux_min >= aux_max) then
            aux_min = ph%pexsi_options%nPoints-1
         endif
      endif

      do i = aux_min,aux_max
         if(abs(ne_vec(i)-ph%n_electrons) <&
            ph%pexsi_options%numElectronPEXSITolerance) then
            ph%mu     = shifts(i)
            converged = .true.

            call MPI_Bcast(tmp_cmplx,bh%nnz_l_sp1,mpi_complex16,i-1,&
                    ph%pexsi_comm_inter_point,ierr)

            call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

            exit
         endif
      enddo
   endif

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Energy density matrix
      factor_min = (ne_vec(aux_max)-ph%n_electrons)/&
                      (ne_vec(aux_max)-ne_vec(aux_min))
      factor_max = (ph%n_electrons-ne_vec(aux_min))/&
                      (ne_vec(aux_max)-ne_vec(aux_min))

      call elsi_allocate(bh,send_buf_cmplx,bh%nnz_l_sp1,"send_buf_cmplx",caller)

      if(ph%pexsi_my_point == aux_min-1) then
         send_buf_cmplx = factor_min*tmp_cmplx
      elseif(ph%pexsi_my_point == aux_max-1) then
         send_buf_cmplx = factor_max*tmp_cmplx
      endif

      call MPI_Allreduce(send_buf_cmplx,tmp_cmplx,bh%nnz_l_sp1,mpi_complex16,&
              mpi_sum,ph%pexsi_comm_inter_point,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(bh,send_buf_cmplx,"send_buf_cmplx")
   endif

   if(.not. (ph%matrix_format == PEXSI_CSC .and. ph%pexsi_my_prow /= 0)) then
      edm = tmp_cmplx
   endif

   call elsi_deallocate(bh,tmp_cmplx,"tmp_cmplx")
   call elsi_deallocate(bh,shifts,"shifts")

   call elsi_get_time(t1)

   write(info_str,"(2X,A)") "Finished energy density matrix calculation"
   call elsi_say(bh,info_str)
   write(info_str,"(2X,A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,info_str)

end subroutine

!>
!! This routine sets default PEXSI parameters.
!!
subroutine elsi_set_pexsi_default(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   character(len=*), parameter :: caller = "elsi_set_pexsi_default"

   ! Use the PEXSI Default options
   call f_ppexsi_set_default_options(ph%pexsi_options)

   ! Use 1 process in symbolic factorization
   ph%pexsi_options%npSymbFact = 1

   ! Use parallel matrix ordering by ParMETIS/PT-SCOTCH
   ! Note: must use serial ordering on some platform (segfault otherwise)
   ph%pexsi_options%ordering = 0

end subroutine

!>
!! This routine cleans up PEXSI.
!!
subroutine elsi_cleanup_pexsi(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_cleanup_pexsi"

   if(ph%pexsi_started) then
      call f_ppexsi_plan_finalize(ph%pexsi_plan,ierr)
      call MPI_Comm_free(ph%pexsi_comm_intra_pole,ierr)
      call MPI_Comm_free(ph%pexsi_comm_inter_pole,ierr)
      call MPI_Comm_free(ph%pexsi_comm_inter_point,ierr)
   endif

   ph%pexsi_started = .false.

end subroutine

end module ELSI_PEXSI
