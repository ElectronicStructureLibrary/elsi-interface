! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Interface to PEXSI.
!!
module ELSI_PEXSI

   use ELSI_CONSTANT, only: UNSET,PEXSI_SOLVER,PEXSI_CSC,DECISION_WIP
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI, only: elsi_stop,elsi_check_mpi,mpi_sum,mpi_real8,mpi_complex16
   use ELSI_OUTPUT, only: elsi_say,elsi_get_time
   use ELSI_PRECISION, only: r8,i4
   use F_PPEXSI_INTERFACE, only: f_ppexsi_plan_initialize,&
       f_ppexsi_plan_finalize,f_ppexsi_set_default_options,&
       f_ppexsi_load_real_hs_matrix,f_ppexsi_load_complex_hs_matrix,&
       f_ppexsi_symbolic_factorize_real_symmetric_matrix,&
       f_ppexsi_symbolic_factorize_complex_symmetric_matrix,&
       f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix,&
       f_ppexsi_inertia_count_real_matrix,&
       f_ppexsi_inertia_count_complex_matrix,&
       f_ppexsi_calculate_fermi_operator_real3,&
       f_ppexsi_calculate_fermi_operator_complex,&
       f_ppexsi_retrieve_real_dm,f_ppexsi_retrieve_complex_dm,&
       f_ppexsi_retrieve_real_edm,f_ppexsi_retrieve_complex_edm

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
!! Initialize PEXSI.
!!
subroutine elsi_init_pexsi(ph,bh)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(inout) :: bh

   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: log_id
   integer(kind=i4) :: ierr
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_init_pexsi"

   if(.not. ph%pexsi_started) then
      ph%pexsi_options%spin = ph%spin_degen

      ! Point parallelization
      ph%pexsi_np_per_point = bh%n_procs/ph%pexsi_options%nPoints
      ph%pexsi_my_point = bh%myid/ph%pexsi_np_per_point
      ph%pexsi_myid_point = mod(bh%myid,ph%pexsi_np_per_point)

      if(ph%pexsi_np_per_pole == UNSET) then
         j = nint(sqrt(real(ph%pexsi_np_per_point,kind=r8)),kind=i4)

         do i = ph%pexsi_np_per_point,j,-1
            if(mod(ph%pexsi_np_per_point,i) == 0) then
               if((ph%pexsi_np_per_point/i) > ph%pexsi_options%numPole) then
                  exit
               end if

               ph%pexsi_np_per_pole = i
            end if
         end do
      end if

      ! Set square-like process grid for selected inversion of each pole
      do i = nint(sqrt(real(ph%pexsi_np_per_pole,kind=r8)),kind=i4),2,-1
         if(mod(ph%pexsi_np_per_pole,i) == 0) then
            exit
         end if
      end do

      ! PEXSI process grid
      ph%pexsi_n_prow = i
      ph%pexsi_n_pcol = ph%pexsi_np_per_pole/i
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
         end if
      end if

      ! Only one process outputs
      if(bh%myid_all == 0) then
         log_id = 0
      else
         log_id = -1
      end if

      ph%pexsi_plan = f_ppexsi_plan_initialize(bh%comm,ph%pexsi_n_prow,&
         ph%pexsi_n_pcol,log_id,ierr)

      if(ierr /= 0) then
         write(msg,"(A)") "Initialization failed"
         call elsi_stop(bh,msg,caller)
      end if

      if(bh%n_lcol_sp == UNSET) then
         bh%n_lcol_sp = bh%n_lcol_sp1
      end if

      ph%pexsi_started = .true.
   end if

end subroutine

!>
!! Interface to PEXSI.
!!
subroutine elsi_solve_pexsi_real(ph,bh,row_ind,col_ptr,ne_vec,ham,ovlp,dm)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8), intent(out) :: ne_vec(ph%pexsi_options%nPoints)
   real(kind=r8), intent(in) :: ham(bh%nnz_l_sp1)
   real(kind=r8), intent(in) :: ovlp(bh%nnz_l_sp1)
   real(kind=r8), intent(out) :: dm(bh%nnz_l_sp1)

   real(kind=r8) :: mu_range
   real(kind=r8) :: shift_width
   real(kind=r8) :: local_energy
   real(kind=r8) :: factor_min
   real(kind=r8) :: factor_max
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: i_step
   integer(kind=i4) :: n_shift
   integer(kind=i4) :: aux_min
   integer(kind=i4) :: aux_max
   integer(kind=i4) :: i
   integer(kind=i4) :: idx
   integer(kind=i4) :: ierr
   logical :: converged
   character(len=200) :: msg

   real(kind=r8), allocatable :: shifts(:)
   real(kind=r8), allocatable :: inertias(:)
   real(kind=r8), allocatable :: ne_lower(:)
   real(kind=r8), allocatable :: ne_upper(:)
   real(kind=r8), allocatable :: tmp(:)
   real(kind=r8), allocatable :: send_buf(:)

   real(kind=r8), external :: ddot

   character(len=*), parameter :: caller = "elsi_solve_pexsi_real"

   write(msg,"(A)") "Starting PEXSI density matrix solver"
   call elsi_say(bh,msg)

   call elsi_get_time(t0)

   ! Load sparse matrices for PEXSI
   if(ph%unit_ovlp) then
      call f_ppexsi_load_real_hs_matrix(ph%pexsi_plan,ph%pexsi_options,&
           ph%n_basis,bh%nnz_g,bh%nnz_l_sp1,bh%n_lcol_sp1,col_ptr,row_ind,ham,&
           1,ovlp,ierr)
   else
      call f_ppexsi_load_real_hs_matrix(ph%pexsi_plan,ph%pexsi_options,&
           ph%n_basis,bh%nnz_g,bh%nnz_l_sp1,bh%n_lcol_sp1,col_ptr,row_ind,ham,&
           0,ovlp,ierr)
   end if

   if(ierr /= 0) then
      write(msg,"(A)") "Failed to load matrices"
      call elsi_stop(bh,msg,caller)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished loading matrices"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ! Symbolic factorization
   if(ph%pexsi_first) then
      call elsi_get_time(t0)

      call f_ppexsi_symbolic_factorize_real_symmetric_matrix(ph%pexsi_plan,&
           ph%pexsi_options,ierr)

      call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(ph%pexsi_plan,&
           ph%pexsi_options,ierr)

      call elsi_get_time(t1)

      write(msg,"(A)") "Finished symbolic factorization"
      call elsi_say(bh,msg)
      write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,msg)
   end if

   if(ierr /= 0) then
      write(msg,"(A)") "Symbolic factorization failed"
      call elsi_stop(bh,msg,caller)
   end if

   ! Inertia counting
   mu_range = ph%pexsi_options%muMax0-ph%pexsi_options%muMin0

   if(mu_range > ph%pexsi_options%muInertiaTolerance) then
      call elsi_get_time(t0)

      i_step = 0
      n_shift = max(10,bh%n_procs/ph%pexsi_np_per_pole)

      call elsi_allocate(bh,shifts,n_shift,"shifts",caller)
      call elsi_allocate(bh,inertias,n_shift,"inertias",caller)
      call elsi_allocate(bh,ne_lower,n_shift,"ne_lower",caller)
      call elsi_allocate(bh,ne_upper,n_shift,"ne_upper",caller)

      do while(i_step < 10 .and. mu_range > ph%pexsi_options%muInertiaTolerance)
         i_step = i_step+1
         shift_width = mu_range/(n_shift-1)
         ne_lower = 0.0_r8
         ne_upper = ph%n_basis*ph%spin_degen

         do i = 1,n_shift
            shifts(i) = ph%pexsi_options%muMin0+(i-1)*shift_width
         end do

         call f_ppexsi_inertia_count_real_matrix(ph%pexsi_plan,&
              ph%pexsi_options,n_shift,shifts,inertias,ierr)

         inertias = inertias*ph%spin_degen*ph%i_wt

         ! Get global inertias
         if(ph%n_spins*ph%n_kpts > 1) then
            call elsi_allocate(bh,send_buf,n_shift,"send_buf",caller)

            if(bh%myid == 0) then
               send_buf = inertias
            else
               send_buf = 0.0_r8
            end if

            call MPI_Allreduce(send_buf,inertias,n_shift,mpi_real8,mpi_sum,&
                 bh%comm_all,ierr)

            call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

            call elsi_deallocate(bh,send_buf,"send_buf")
         end if

         idx = ceiling(3*ph%pexsi_options%temperature/shift_width)

         do i = idx+1,n_shift
            ne_lower(i) = 0.5_r8*(inertias(i-idx)+inertias(i))
            ne_upper(i-idx) = ne_lower(i)
         end do

         aux_min = 1
         aux_max = n_shift

         do i = 2,n_shift-1
            if(ne_upper(i) < ph%n_electrons&
               .and. ne_upper(i+1) >= ph%n_electrons) then
               aux_min = i
            end if

            if(ne_lower(i) > ph%n_electrons&
               .and. ne_lower(i-1) <= ph%n_electrons) then
               aux_max = i
            end if
         end do

         if(aux_min == 1 .and. aux_max == n_shift) then
            exit
         else
            ph%pexsi_options%muMin0 = shifts(aux_min)
            ph%pexsi_options%muMax0 = shifts(aux_max)

            mu_range = ph%pexsi_options%muMax0-ph%pexsi_options%muMin0
         end if
      end do

      call elsi_deallocate(bh,shifts,"shifts")
      call elsi_deallocate(bh,inertias,"inertias")
      call elsi_deallocate(bh,ne_lower,"ne_lower")
      call elsi_deallocate(bh,ne_upper,"ne_upper")

      call elsi_get_time(t1)

      write(msg,"(A)") "Finished inertia counting"
      call elsi_say(bh,msg)
      write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,msg)

      if(ierr /= 0) then
         write(msg,"(A)") "Inertia counting failed"
         call elsi_stop(bh,msg,caller)
      end if
   end if

   ! Fermi operator expansion
   call elsi_get_time(t0)

   shift_width = mu_range/(ph%pexsi_options%nPoints+1)

   call elsi_allocate(bh,shifts,ph%pexsi_options%nPoints,"shifts",caller)

   do i = 1,ph%pexsi_options%nPoints
      shifts(i) = ph%pexsi_options%muMin0+i*shift_width
   end do

   do i = 1,ph%pexsi_options%nPoints
      ph%mu = shifts(i)

      if(ph%pexsi_my_point == i-1) then
         call f_ppexsi_calculate_fermi_operator_real3(ph%pexsi_plan,&
              ph%pexsi_options,ph%n_electrons,ph%mu,ph%pexsi_ne,ierr)
      end if
   end do

   call elsi_allocate(bh,send_buf,ph%pexsi_options%nPoints,"send_buf",caller)

   send_buf(ph%pexsi_my_point+1) = ph%pexsi_ne*ph%i_wt

   call MPI_Allreduce(send_buf,ne_vec,ph%pexsi_options%nPoints,mpi_real8,&
        mpi_sum,ph%pexsi_comm_inter_point,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   ! Get global number of electrons
   if(ph%n_spins*ph%n_kpts > 1) then
      if(bh%myid == 0) then
         send_buf = ne_vec
      else
         send_buf = 0.0_r8
      end if

      call MPI_Allreduce(send_buf,ne_vec,ph%pexsi_options%nPoints,mpi_real8,&
           mpi_sum,bh%comm_all,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   end if

   call elsi_deallocate(bh,send_buf,"send_buf")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished Fermi operator calculation"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   if(ierr /= 0) then
      write(msg,"(A)") "Fermi operator calculation failed"
      call elsi_stop(bh,msg,caller)
   end if

   if(ph%decision_stage == DECISION_WIP) then
      ph%decision_data(PEXSI_SOLVER) = t1-t0
   end if

   ! Get density matrix
   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp,bh%nnz_l_sp1,"tmp",caller)

   call f_ppexsi_retrieve_real_dm(ph%pexsi_plan,tmp,local_energy,ierr)

   if(ierr /= 0) then
      write(msg,"(A)") "Failed to get density matirx"
      call elsi_stop(bh,msg,caller)
   end if

   ! Check convergence
   converged = .false.
   aux_min = 0
   aux_max = ph%pexsi_options%nPoints+1

   do i = 1,ph%pexsi_options%nPoints
      if(ne_vec(i)&
         < ph%n_electrons-ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMin0 = shifts(i)
         aux_min = i
      end if
   end do

   do i = ph%pexsi_options%nPoints,1,-1
      if(ne_vec(i)&
         > ph%n_electrons+ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMax0 = shifts(i)
         aux_max = i
      end if
   end do

   if(ph%pexsi_options%nPoints == 1) then
      ! Scale density matrix
      tmp = (ph%n_electrons/ph%pexsi_ne)*tmp
      converged = .true.
      ph%mu = shifts(1)
   else
      ! Safety check
      if(aux_min == 0) then
         aux_min = 1

         if(aux_max <= aux_min) then
            aux_max = 2
         end if
      end if

      if(aux_max == ph%pexsi_options%nPoints+1) then
         aux_max = ph%pexsi_options%nPoints

         if(aux_min >= aux_max) then
            aux_min = ph%pexsi_options%nPoints-1
         end if
      end if

      do i = aux_min,aux_max
         if(abs(ne_vec(i)-ph%n_electrons)&
            < ph%pexsi_options%numElectronPEXSITolerance) then
            ph%mu = shifts(i)
            converged = .true.

            call MPI_Bcast(tmp,bh%nnz_l_sp1,mpi_real8,i-1,&
                 ph%pexsi_comm_inter_point,ierr)

            call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

            exit
         end if
      end do
   end if

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Chemical potential
      ph%mu = shifts(aux_min)+(ph%n_electrons-ne_vec(aux_min))&
         /(ne_vec(aux_max)-ne_vec(aux_min))*(shifts(aux_max)-shifts(aux_min))

      ! Density matrix
      factor_min = (ne_vec(aux_max)-ph%n_electrons)&
         /(ne_vec(aux_max)-ne_vec(aux_min))
      factor_max = (ph%n_electrons-ne_vec(aux_min))&
         /(ne_vec(aux_max)-ne_vec(aux_min))

      call elsi_allocate(bh,send_buf,bh%nnz_l_sp1,"send_buf",caller)

      if(ph%pexsi_my_point == aux_min-1) then
         send_buf = factor_min*tmp
      else if(ph%pexsi_my_point == aux_max-1) then
         send_buf = factor_max*tmp
      end if

      call MPI_Allreduce(send_buf,tmp,bh%nnz_l_sp1,mpi_real8,mpi_sum,&
           ph%pexsi_comm_inter_point,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(bh,send_buf,"send_buf")
   end if

   if(.not. (ph%matrix_format == PEXSI_CSC .and. ph%pexsi_my_prow /= 0)) then
      dm = tmp
   end if

   call elsi_deallocate(bh,tmp,"tmp")
   call elsi_deallocate(bh,shifts,"shifts")

   ! Compute energy = Tr(H*DM)
   if(ph%pexsi_my_prow == 0) then
      local_energy = ddot(bh%nnz_l_sp1,ham,1,dm,1)

      call MPI_Reduce(local_energy,ph%ebs,1,mpi_real8,mpi_sum,0,&
           ph%pexsi_comm_intra_pole,ierr)

      call elsi_check_mpi(bh,"MPI_Reduce",ierr,caller)
   end if

   call MPI_Bcast(ph%ebs,1,mpi_real8,0,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished density matrix correction"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%pexsi_first = .false.

end subroutine

!>
!! Compute the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_pexsi_real(ph,bh,ne_vec,edm)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: ne_vec(ph%pexsi_options%nPoints)
   real(kind=r8), intent(out) :: edm(bh%nnz_l_sp1)

   real(kind=r8) :: mu_range
   real(kind=r8) :: shift_width
   real(kind=r8) :: local_energy
   real(kind=r8) :: factor_min
   real(kind=r8) :: factor_max
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: aux_min
   integer(kind=i4) :: aux_max
   integer(kind=i4) :: i
   integer(kind=i4) :: ierr
   logical :: converged
   character(len=200) :: msg

   real(kind=r8), allocatable :: shifts(:)
   real(kind=r8), allocatable :: tmp(:)
   real(kind=r8), allocatable :: send_buf(:)

   character(len=*), parameter :: caller = "elsi_compute_edm_pexsi_real"

   call elsi_get_time(t0)

   ! Get energy density matrix
   call elsi_allocate(bh,tmp,bh%nnz_l_sp1,"tmp",caller)

   call f_ppexsi_retrieve_real_edm(ph%pexsi_plan,ph%pexsi_options,tmp,&
        local_energy,ierr)

   if(ierr /= 0) then
      write(msg,"(A)") "Failed to get energy density matirx"
      call elsi_stop(bh,msg,caller)
   end if

   ! Check convergence
   mu_range = ph%pexsi_options%muMax0-ph%pexsi_options%muMin0
   shift_width = mu_range/(ph%pexsi_options%nPoints+1)
   converged = .false.
   aux_min = 0
   aux_max = ph%pexsi_options%nPoints+1

   call elsi_allocate(bh,shifts,ph%pexsi_options%nPoints,"shifts",caller)

   do i = 1,ph%pexsi_options%nPoints
      shifts(i) = ph%pexsi_options%muMin0+i*shift_width
   end do

   do i = 1,ph%pexsi_options%nPoints
      if(ne_vec(i)&
         < ph%n_electrons-ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMin0 = shifts(i)
         aux_min = i
      end if
   end do

   do i = ph%pexsi_options%nPoints,1,-1
      if(ne_vec(i)&
         > ph%n_electrons+ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMax0 = shifts(i)
         aux_max = i
      end if
   end do

   if(ph%pexsi_options%nPoints == 1) then
      ! Scale energy density matrix
      tmp = (ph%n_electrons/ph%pexsi_ne)*tmp
      converged = .true.
      ph%mu = shifts(1)
   else
      ! Safety check
      if(aux_min == 0) then
         aux_min = 1

         if(aux_max <= aux_min) then
            aux_max = 2
         end if
      end if

      if(aux_max == ph%pexsi_options%nPoints+1) then
         aux_max = ph%pexsi_options%nPoints

         if(aux_min >= aux_max) then
            aux_min = ph%pexsi_options%nPoints-1
         end if
      end if

      do i = aux_min,aux_max
         if(abs(ne_vec(i)-ph%n_electrons)&
            < ph%pexsi_options%numElectronPEXSITolerance) then
            ph%mu = shifts(i)
            converged = .true.

            call MPI_Bcast(tmp,bh%nnz_l_sp1,mpi_real8,i-1,&
                 ph%pexsi_comm_inter_point,ierr)

            call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

            exit
         end if
      end do
   end if

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Energy density matrix
      factor_min = (ne_vec(aux_max)-ph%n_electrons)&
         /(ne_vec(aux_max)-ne_vec(aux_min))
      factor_max = (ph%n_electrons-ne_vec(aux_min))&
         /(ne_vec(aux_max)-ne_vec(aux_min))

      call elsi_allocate(bh,send_buf,bh%nnz_l_sp1,"send_buf",caller)

      if(ph%pexsi_my_point == aux_min-1) then
         send_buf = factor_min*tmp
      else if(ph%pexsi_my_point == aux_max-1) then
         send_buf = factor_max*tmp
      end if

      call MPI_Allreduce(send_buf,tmp,bh%nnz_l_sp1,mpi_real8,mpi_sum,&
           ph%pexsi_comm_inter_point,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(bh,send_buf,"send_buf")
   end if

   if(.not. (ph%matrix_format == PEXSI_CSC .and. ph%pexsi_my_prow /= 0)) then
      edm = tmp
   end if

   call elsi_deallocate(bh,tmp,"tmp")
   call elsi_deallocate(bh,shifts,"shifts")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished energy density matrix calculation"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Interface to PEXSI.
!!
subroutine elsi_solve_pexsi_cmplx(ph,bh,row_ind,col_ptr,ne_vec,ham,ovlp,dm)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: row_ind(bh%nnz_l_sp1)
   integer(kind=i4), intent(in) :: col_ptr(bh%n_lcol_sp1+1)
   real(kind=r8), intent(out) :: ne_vec(ph%pexsi_options%nPoints)
   complex(kind=r8), intent(in) :: ham(bh%nnz_l_sp1)
   complex(kind=r8), intent(inout) :: ovlp(bh%nnz_l_sp1)
   complex(kind=r8), intent(out) :: dm(bh%nnz_l_sp1)

   complex(kind=r8) :: local_cmplx
   real(kind=r8) :: ne_drv
   real(kind=r8) :: mu_range
   real(kind=r8) :: shift_width
   real(kind=r8) :: local_energy
   real(kind=r8) :: factor_min
   real(kind=r8) :: factor_max
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: i_step
   integer(kind=i4) :: n_shift
   integer(kind=i4) :: aux_min
   integer(kind=i4) :: aux_max
   integer(kind=i4) :: i
   integer(kind=i4) :: idx
   integer(kind=i4) :: ierr
   logical :: converged
   character(len=200) :: msg

   complex(kind=r8), allocatable :: tmp(:)
   complex(kind=r8), allocatable :: send_buf_cmplx(:)
   real(kind=r8), allocatable :: shifts(:)
   real(kind=r8), allocatable :: inertias(:)
   real(kind=r8), allocatable :: ne_lower(:)
   real(kind=r8), allocatable :: ne_upper(:)
   real(kind=r8), allocatable :: send_buf(:)

   complex(kind=r8), external :: zdotc

   character(len=*), parameter :: caller = "elsi_solve_pexsi_cmplx"

   write(msg,"(A)") "Starting PEXSI density matrix solver"
   call elsi_say(bh,msg)

   call elsi_get_time(t0)

   ! Load sparse matrices for PEXSI
   if(ph%unit_ovlp) then
      call f_ppexsi_load_complex_hs_matrix(ph%pexsi_plan,ph%pexsi_options,&
           ph%n_basis,bh%nnz_g,bh%nnz_l_sp1,bh%n_lcol_sp1,col_ptr,row_ind,ham,&
           1,ovlp,ierr)
   else
      call f_ppexsi_load_complex_hs_matrix(ph%pexsi_plan,ph%pexsi_options,&
           ph%n_basis,bh%nnz_g,bh%nnz_l_sp1,bh%n_lcol_sp1,col_ptr,row_ind,ham,&
           0,ovlp,ierr)
   end if

   if(ierr /= 0) then
      write(msg,"(A)") "Failed to load matirces"
      call elsi_stop(bh,msg,caller)
   end if

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished loading matrices"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ! Symbolic factorization
   if(ph%pexsi_first) then
      call elsi_get_time(t0)

      call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(ph%pexsi_plan,&
           ph%pexsi_options,ierr)

      call f_ppexsi_symbolic_factorize_complex_unsymmetric_matrix(&
           ph%pexsi_plan,ph%pexsi_options,ovlp,ierr)

      call elsi_get_time(t1)

      write(msg,"(A)") "Finished symbolic factorization"
      call elsi_say(bh,msg)
      write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,msg)
   end if

   if(ierr /= 0) then
      write(msg,"(A)") "Symbolic factorization failed"
      call elsi_stop(bh,msg,caller)
   end if

   ! Inertia counting
   mu_range = ph%pexsi_options%muMax0-ph%pexsi_options%muMin0

   if(mu_range > ph%pexsi_options%muInertiaTolerance) then
      call elsi_get_time(t0)

      i_step = 0
      n_shift = max(10,bh%n_procs/ph%pexsi_np_per_pole)

      call elsi_allocate(bh,shifts,n_shift,"shifts",caller)
      call elsi_allocate(bh,inertias,n_shift,"inertias",caller)
      call elsi_allocate(bh,ne_lower,n_shift,"ne_lower",caller)
      call elsi_allocate(bh,ne_upper,n_shift,"ne_upper",caller)

      do while(i_step < 10 .and. mu_range > ph%pexsi_options%muInertiaTolerance)
         i_step = i_step+1
         shift_width = mu_range/(n_shift-1)
         ne_lower = 0.0_r8
         ne_upper = ph%n_basis*ph%spin_degen

         do i = 1,n_shift
            shifts(i) = ph%pexsi_options%muMin0+(i-1)*shift_width
         end do

         call f_ppexsi_inertia_count_complex_matrix(ph%pexsi_plan,&
              ph%pexsi_options,n_shift,shifts,inertias,ierr)

         inertias = inertias*ph%spin_degen*ph%i_wt

         ! Get global inertias
         if(ph%n_spins*ph%n_kpts > 1) then
            call elsi_allocate(bh,send_buf,n_shift,"send_buf",caller)

            if(bh%myid == 0) then
               send_buf = inertias
            else
               send_buf = 0.0_r8
            end if

            call MPI_Allreduce(send_buf,inertias,n_shift,mpi_real8,mpi_sum,&
                 bh%comm_all,ierr)

            call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

            call elsi_deallocate(bh,send_buf,"send_buf")
         end if

         idx = ceiling(3*ph%pexsi_options%temperature/shift_width)

         do i = idx+1,n_shift
            ne_lower(i) = 0.5_r8*(inertias(i-idx)+inertias(i))
            ne_upper(i-idx) = ne_lower(i)
         end do

         aux_min = 1
         aux_max = n_shift

         do i = 2,n_shift-1
            if(ne_upper(i) < ph%n_electrons&
               .and. ne_upper(i+1) >= ph%n_electrons) then
               aux_min = i
            end if

            if(ne_lower(i) > ph%n_electrons&
               .and. ne_lower(i-1) <= ph%n_electrons) then
               aux_max = i
            end if
         end do

         if(aux_min == 1 .and. aux_max == n_shift) then
            exit
         else
            ph%pexsi_options%muMin0 = shifts(aux_min)
            ph%pexsi_options%muMax0 = shifts(aux_max)

            mu_range = ph%pexsi_options%muMax0-ph%pexsi_options%muMin0
         end if
      end do

      call elsi_deallocate(bh,shifts,"shifts")
      call elsi_deallocate(bh,inertias,"inertias")
      call elsi_deallocate(bh,ne_lower,"ne_lower")
      call elsi_deallocate(bh,ne_upper,"ne_upper")

      call elsi_get_time(t1)

      write(msg,"(A)") "Finished inertia counting"
      call elsi_say(bh,msg)
      write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
      call elsi_say(bh,msg)

      if(ierr /= 0) then
         write(msg,"(A)") "Inertia counting failed"
         call elsi_stop(bh,msg,caller)
      end if
   end if

   ! Fermi operator expansion
   call elsi_get_time(t0)

   shift_width = mu_range/(ph%pexsi_options%nPoints+1)

   call elsi_allocate(bh,shifts,ph%pexsi_options%nPoints,"shifts",caller)

   do i = 1,ph%pexsi_options%nPoints
      shifts(i) = ph%pexsi_options%muMin0+i*shift_width
   end do

   do i = 1,ph%pexsi_options%nPoints
      ph%mu = shifts(i)

      if(ph%pexsi_my_point == i-1) then
         call f_ppexsi_calculate_fermi_operator_complex(ph%pexsi_plan,&
              ph%pexsi_options,ph%mu,ph%n_electrons,ph%pexsi_ne,ne_drv,ierr)
      end if
   end do

   call elsi_allocate(bh,send_buf,ph%pexsi_options%nPoints,"send_buf",caller)

   send_buf(ph%pexsi_my_point+1) = ph%pexsi_ne*ph%i_wt

   call MPI_Allreduce(send_buf,ne_vec,ph%pexsi_options%nPoints,mpi_real8,&
        mpi_sum,ph%pexsi_comm_inter_point,ierr)

   call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

   ! Get global number of electrons
   if(ph%n_spins*ph%n_kpts > 1) then
      if(bh%myid == 0) then
         send_buf = ne_vec
      else
         send_buf = 0.0_r8
      end if

      call MPI_Allreduce(send_buf,ne_vec,ph%pexsi_options%nPoints,mpi_real8,&
           mpi_sum,bh%comm_all,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)
   end if

   call elsi_deallocate(bh,send_buf,"send_buf")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished Fermi operator calculation"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   if(ierr /= 0) then
      write(msg,"(A)") "Fermi operator calculation failed"
      call elsi_stop(bh,msg,caller)
   end if

   if(ph%decision_stage == DECISION_WIP) then
      ph%decision_data(PEXSI_SOLVER) = t1-t0
   end if

   ! Get density matrix
   call elsi_get_time(t0)

   call elsi_allocate(bh,tmp,bh%nnz_l_sp1,"tmp",caller)

   call f_ppexsi_retrieve_complex_dm(ph%pexsi_plan,tmp,local_energy,ierr)

   if(ierr /= 0) then
      write(msg,"(A)") "Failed to get density matirx"
      call elsi_stop(bh,msg,caller)
   end if

   ! Check convergence
   converged = .false.
   aux_min = 0
   aux_max = ph%pexsi_options%nPoints+1

   do i = 1,ph%pexsi_options%nPoints
      if(ne_vec(i)&
         < ph%n_electrons-ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMin0 = shifts(i)
         aux_min = i
      end if
   end do

   do i = ph%pexsi_options%nPoints,1,-1
      if(ne_vec(i)&
         > ph%n_electrons+ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMax0 = shifts(i)
         aux_max = i
      end if
   end do

   if(ph%pexsi_options%nPoints == 1) then
      ! Scale density matrix
      tmp = (ph%n_electrons/ph%pexsi_ne)*tmp
      converged = .true.
      ph%mu = shifts(1)
   else
      ! Safety check
      if(aux_min == 0) then
         aux_min = 1

         if(aux_max <= aux_min) then
            aux_max = 2
         end if
      end if

      if(aux_max == ph%pexsi_options%nPoints+1) then
         aux_max = ph%pexsi_options%nPoints

         if(aux_min >= aux_max) then
            aux_min = ph%pexsi_options%nPoints-1
         end if
      end if

      do i = aux_min,aux_max
         if(abs(ne_vec(i)-ph%n_electrons)&
            < ph%pexsi_options%numElectronPEXSITolerance) then
            ph%mu = shifts(i)
            converged = .true.

            call MPI_Bcast(tmp,bh%nnz_l_sp1,mpi_complex16,i-1,&
                 ph%pexsi_comm_inter_point,ierr)

            call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

            exit
         end if
      end do
   end if

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Chemical potential
      ph%mu = shifts(aux_min)+(ph%n_electrons-ne_vec(aux_min))&
         /(ne_vec(aux_max)-ne_vec(aux_min))*(shifts(aux_max)-shifts(aux_min))

      ! Density matrix
      factor_min = (ne_vec(aux_max)-ph%n_electrons)&
         /(ne_vec(aux_max)-ne_vec(aux_min))
      factor_max = (ph%n_electrons-ne_vec(aux_min))&
         /(ne_vec(aux_max)-ne_vec(aux_min))

      call elsi_allocate(bh,send_buf_cmplx,bh%nnz_l_sp1,"send_buf_cmplx",caller)

      if(ph%pexsi_my_point == aux_min-1) then
         send_buf_cmplx = factor_min*tmp
      else if(ph%pexsi_my_point == aux_max-1) then
         send_buf_cmplx = factor_max*tmp
      end if

      call MPI_Allreduce(send_buf_cmplx,tmp,bh%nnz_l_sp1,mpi_complex16,mpi_sum,&
           ph%pexsi_comm_inter_point,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(bh,send_buf_cmplx,"send_buf_cmplx")
   end if

   if(.not. (ph%matrix_format == PEXSI_CSC .and. ph%pexsi_my_prow /= 0)) then
      dm = tmp
   end if

   call elsi_deallocate(bh,tmp,"tmp")
   call elsi_deallocate(bh,shifts,"shifts")

   ! Compute energy = Tr(H*DM)
   if(ph%pexsi_my_prow == 0) then
      local_cmplx = zdotc(bh%nnz_l_sp1,ham,1,dm,1)
      local_energy = real(local_cmplx,kind=r8)

      call MPI_Reduce(local_energy,ph%ebs,1,mpi_real8,mpi_sum,0,&
           ph%pexsi_comm_intra_pole,ierr)

      call elsi_check_mpi(bh,"MPI_Reduce",ierr,caller)
   end if

   call MPI_Bcast(ph%ebs,1,mpi_real8,0,bh%comm,ierr)

   call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished density matrix correction"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

   ph%pexsi_first = .false.

end subroutine

!>
!! Compute the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_pexsi_cmplx(ph,bh,ne_vec,edm)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: ne_vec(ph%pexsi_options%nPoints)
   complex(kind=r8), intent(out) :: edm(bh%nnz_l_sp1)

   real(kind=r8) :: mu_range
   real(kind=r8) :: shift_width
   real(kind=r8) :: local_energy
   real(kind=r8) :: factor_min
   real(kind=r8) :: factor_max
   real(kind=r8) :: t0
   real(kind=r8) :: t1
   integer(kind=i4) :: aux_min
   integer(kind=i4) :: aux_max
   integer(kind=i4) :: i
   integer(kind=i4) :: ierr
   logical :: converged
   character(len=200) :: msg

   complex(kind=r8), allocatable :: tmp(:)
   complex(kind=r8), allocatable :: send_buf(:)
   real(kind=r8), allocatable :: shifts(:)

   character(len=*), parameter :: caller = "elsi_compute_edm_pexsi_cmplx"

   call elsi_get_time(t0)

   ! Get energy density matrix
   call elsi_allocate(bh,tmp,bh%nnz_l_sp1,"tmp",caller)

   call f_ppexsi_retrieve_complex_edm(ph%pexsi_plan,ph%pexsi_options,tmp,&
        local_energy,ierr)

   if(ierr /= 0) then
      write(msg,"(A)") "Failed to get energy density matirx"
      call elsi_stop(bh,msg,caller)
   end if

   ! Check convergence
   mu_range = ph%pexsi_options%muMax0-ph%pexsi_options%muMin0
   shift_width = mu_range/(ph%pexsi_options%nPoints+1)
   converged = .false.
   aux_min = 0
   aux_max = ph%pexsi_options%nPoints+1

   call elsi_allocate(bh,shifts,ph%pexsi_options%nPoints,"shifts",caller)

   do i = 1,ph%pexsi_options%nPoints
      shifts(i) = ph%pexsi_options%muMin0+i*shift_width
   end do

   do i = 1,ph%pexsi_options%nPoints
      if(ne_vec(i)&
         < ph%n_electrons-ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMin0 = shifts(i)
         aux_min = i
      end if
   end do

   do i = ph%pexsi_options%nPoints,1,-1
      if(ne_vec(i)&
         > ph%n_electrons+ph%pexsi_options%numElectronPEXSITolerance) then
         ph%pexsi_options%muMax0 = shifts(i)
         aux_max = i
      end if
   end do

   if(ph%pexsi_options%nPoints == 1) then
      ! Scale energy density matrix
      tmp = (ph%n_electrons/ph%pexsi_ne)*tmp
      converged = .true.
      ph%mu = shifts(1)
   else
      ! Safety check
      if(aux_min == 0) then
         aux_min = 1

         if(aux_max <= aux_min) then
            aux_max = 2
         end if
      end if

      if(aux_max == ph%pexsi_options%nPoints+1) then
         aux_max = ph%pexsi_options%nPoints

         if(aux_min >= aux_max) then
            aux_min = ph%pexsi_options%nPoints-1
         end if
      end if

      do i = aux_min,aux_max
         if(abs(ne_vec(i)-ph%n_electrons)&
            < ph%pexsi_options%numElectronPEXSITolerance) then
            ph%mu = shifts(i)
            converged = .true.

            call MPI_Bcast(tmp,bh%nnz_l_sp1,mpi_complex16,i-1,&
                 ph%pexsi_comm_inter_point,ierr)

            call elsi_check_mpi(bh,"MPI_Bcast",ierr,caller)

            exit
         end if
      end do
   end if

   ! Adjust to exact number of electrons
   if(.not. converged) then
      ! Energy density matrix
      factor_min = (ne_vec(aux_max)-ph%n_electrons)&
         /(ne_vec(aux_max)-ne_vec(aux_min))
      factor_max = (ph%n_electrons-ne_vec(aux_min))&
         /(ne_vec(aux_max)-ne_vec(aux_min))

      call elsi_allocate(bh,send_buf,bh%nnz_l_sp1,"send_buf",caller)

      if(ph%pexsi_my_point == aux_min-1) then
         send_buf = factor_min*tmp
      else if(ph%pexsi_my_point == aux_max-1) then
         send_buf = factor_max*tmp
      end if

      call MPI_Allreduce(send_buf,tmp,bh%nnz_l_sp1,mpi_complex16,mpi_sum,&
           ph%pexsi_comm_inter_point,ierr)

      call elsi_check_mpi(bh,"MPI_Allreduce",ierr,caller)

      call elsi_deallocate(bh,send_buf,"send_buf")
   end if

   if(.not. (ph%matrix_format == PEXSI_CSC .and. ph%pexsi_my_prow /= 0)) then
      edm = tmp
   end if

   call elsi_deallocate(bh,tmp,"tmp")
   call elsi_deallocate(bh,shifts,"shifts")

   call elsi_get_time(t1)

   write(msg,"(A)") "Finished energy density matrix calculation"
   call elsi_say(bh,msg)
   write(msg,"(A,F10.3,A)") "| Time :",t1-t0," s"
   call elsi_say(bh,msg)

end subroutine

!>
!! Set default PEXSI parameters.
!!
subroutine elsi_set_pexsi_default(ph)

   implicit none

   type(elsi_param_t), intent(inout) :: ph

   character(len=*), parameter :: caller = "elsi_set_pexsi_default"

   call f_ppexsi_set_default_options(ph%pexsi_options)

   ! Pole expansion method 2 for now (PEXSI default 3)
   ph%pexsi_options%method = 2

end subroutine

!>
!! Clean up PEXSI.
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
   end if

   ph%pexsi_first = .true.
   ph%pexsi_started = .false.

end subroutine

end module ELSI_PEXSI
