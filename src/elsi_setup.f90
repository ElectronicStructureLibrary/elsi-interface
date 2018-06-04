! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides routines for setting up an ELSI instance.
!!
module ELSI_SETUP

   use ELSI_CONSTANTS,     only: ELPA_SOLVER,OMM_SOLVER,SIPS_SOLVER,DMP_SOLVER,&
                                 PEXSI_SOLVER,UNSET,SINGLE_PROC,MULTI_PROC,&
                                 BLACS_DENSE,PEXSI_CSC,SIESTA_CSC
   use ELSI_DATATYPE,      only: elsi_handle,elsi_param_t,elsi_basic_t
   use ELSI_DMP,           only: elsi_cleanup_dmp
   use ELSI_ELPA,          only: elsi_cleanup_elpa
   use ELSI_IO,            only: elsi_say,elsi_final_print,fjson_close_file,&
                                 fjson_finish_array,fjson_reset_fj_handle
   use ELSI_MALLOC,        only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,           only: elsi_stop
   use ELSI_OMM,           only: elsi_cleanup_omm
   use ELSI_PEXSI,         only: elsi_set_pexsi_default,elsi_cleanup_pexsi
   use ELSI_PRECISION,     only: r8,i4
   use ELSI_SIPS,          only: elsi_cleanup_sips
   use ELSI_UTILS,         only: elsi_check_init,elsi_reset_param,&
                                 elsi_reset_basic

   implicit none

   private

   public :: elsi_init
   public :: elsi_finalize
   public :: elsi_set_mpi
   public :: elsi_set_mpi_global
   public :: elsi_set_spin
   public :: elsi_set_kpoint
   public :: elsi_set_blacs
   public :: elsi_set_csc
   public :: elsi_cleanup

contains

!>
!! This routine initializes ELSI with the solver, parallel mode, matrix storage
!! format, number of basis functions (global size of the Hamiltonian matrix),
!! number of electrons, and number of states.
!!
subroutine elsi_init(eh,solver,parallel_mode,matrix_format,n_basis,n_electron,&
              n_state)

   implicit none

   type(elsi_handle), intent(out) :: eh            !< Handle
   integer(kind=i4),  intent(in)  :: solver        !< AUTO,ELPA,LIBOMM,PEXSI,CHESS,SIPS,DMP
   integer(kind=i4),  intent(in)  :: parallel_mode !< SINGLE_PROC,MULTI_PROC
   integer(kind=i4),  intent(in)  :: matrix_format !< BLACS_DENSE,PEXSI_CSC,SIESTA_CSC
   integer(kind=i4),  intent(in)  :: n_basis       !< Number of basis functions
   real(kind=r8),     intent(in)  :: n_electron    !< Number of electrons
   integer(kind=i4),  intent(in)  :: n_state       !< Number of states

   character(len=40), parameter :: caller = "elsi_init"

   ! For safety
   call elsi_cleanup(eh)

   eh%handle_init       = .true.
   eh%ph%n_basis        = n_basis
   eh%ph%n_good         = n_basis
   eh%ph%n_electrons    = n_electron
   eh%ph%n_states       = n_state
   eh%ph%n_states_solve = n_state
   eh%ph%omm_n_states   = nint(n_electron/2.0_r8)
   eh%ph%dmp_n_states   = nint(n_electron/2.0_r8)
   eh%ph%solver         = solver
   eh%ph%matrix_format  = matrix_format
   eh%ph%parallel_mode  = parallel_mode
   eh%ph%n_calls        = 0

   if(parallel_mode == SINGLE_PROC) then
      eh%bh%n_lrow      = n_basis
      eh%bh%n_lcol      = n_basis
      eh%bh%blk         = n_basis
      eh%bh%n_prow      = 1
      eh%bh%n_pcol      = 1
      eh%bh%myid        = 0
      eh%bh%n_procs     = 1
      eh%bh%myid_all    = 0
      eh%bh%n_procs_all = 1
   endif

   if(solver == PEXSI_SOLVER) then
      call elsi_set_pexsi_default(eh%ph)
   endif

end subroutine

!>
!! This routine sets the unit MPI communicator.
!!
subroutine elsi_set_mpi(eh,comm)

   implicit none

   type(elsi_handle), intent(inout) :: eh   !< Handle
   integer(kind=i4),  intent(in)    :: comm !< Unit communicator

   integer(kind=i4) :: ierr

   character(len=40), parameter :: caller = "elsi_set_mpi"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(eh%ph%parallel_mode == MULTI_PROC) then
      eh%bh%comm = comm

      call MPI_Comm_rank(comm,eh%bh%myid,ierr)
      call MPI_Comm_size(comm,eh%bh%n_procs,ierr)

      eh%bh%mpi_ready = .true.
   endif

end subroutine

!>
!! This routine sets the global MPI communicator.
!!
subroutine elsi_set_mpi_global(eh,comm_all)

   implicit none

   type(elsi_handle), intent(inout) :: eh       !< Handle
   integer(kind=i4),  intent(in)    :: comm_all !< Global communicator

   integer(kind=i4) :: ierr

   character(len=40), parameter :: caller = "elsi_set_mpi_global"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(eh%ph%parallel_mode == MULTI_PROC) then
      eh%bh%comm_all = comm_all

      call MPI_Comm_rank(comm_all,eh%bh%myid_all,ierr)
      call MPI_Comm_size(comm_all,eh%bh%n_procs_all,ierr)

      eh%bh%mpi_all_ready = .true.
   endif

end subroutine

!>
!! This routine sets the spin information.
!!
subroutine elsi_set_spin(eh,n_spin,i_spin)

   implicit none

   type(elsi_handle), intent(inout) :: eh     !< Handle
   integer(kind=i4),  intent(in)    :: n_spin !< Number of spin channels
   integer(kind=i4),  intent(in)    :: i_spin !< Spin index

   eh%ph%n_spins = n_spin
   eh%ph%i_spin  = i_spin

end subroutine

!>
!! This routine sets the k-point information.
!!
subroutine elsi_set_kpoint(eh,n_kpt,i_kpt,weight)

   implicit none

   type(elsi_handle), intent(inout) :: eh     !< Handle
   integer(kind=i4),  intent(in)    :: n_kpt  !< Number of k-points
   integer(kind=i4),  intent(in)    :: i_kpt  !< K-point index
   real(kind=r8),     intent(in)    :: weight !< Weight

   eh%ph%n_kpts   = n_kpt
   eh%ph%i_kpt    = i_kpt
   eh%ph%i_weight = weight

end subroutine

!>
!! This routine sets the BLACS context and the block size.
!!
subroutine elsi_set_blacs(eh,blacs_ctxt,block_size)

   implicit none

   type(elsi_handle), intent(inout) :: eh         !< Handle
   integer(kind=i4),  intent(in)    :: blacs_ctxt !< BLACS context
   integer(kind=i4),  intent(in)    :: block_size !< Block size

   integer(kind=i4) :: i
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: ierr

   integer(kind=i4), external :: numroc

   character(len=40), parameter :: caller = "elsi_set_blacs"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(eh%ph%parallel_mode == MULTI_PROC) then
      eh%bh%blacs_ctxt = blacs_ctxt
      eh%bh%blk        = block_size

      ! Get processor grid information
      call BLACS_Gridinfo(blacs_ctxt,eh%bh%n_prow,eh%bh%n_pcol,eh%bh%my_prow,&
              eh%bh%my_pcol)

      ! Get local size of matrix
      eh%bh%n_lrow = numroc(eh%ph%n_basis,eh%bh%blk,eh%bh%my_prow,0,&
                        eh%bh%n_prow)
      eh%bh%n_lcol = numroc(eh%ph%n_basis,eh%bh%blk,eh%bh%my_pcol,0,&
                        eh%bh%n_pcol)

      ! Get BLACS descriptor
      call descinit(eh%bh%desc,eh%ph%n_basis,eh%ph%n_basis,eh%bh%blk,eh%bh%blk,&
              0,0,eh%bh%blacs_ctxt,max(1,eh%bh%n_lrow),ierr)

      ! Create global-local mapping
      call elsi_allocate(eh%bh,eh%row_map,eh%ph%n_basis,"row_map",caller)
      call elsi_allocate(eh%bh,eh%col_map,eh%ph%n_basis,"col_map",caller)

      i_row = 0
      i_col = 0

      do i = 1,eh%ph%n_basis
         if(mod((i-1)/eh%bh%blk,eh%bh%n_prow) == eh%bh%my_prow) then
            i_row         = i_row+1
            eh%row_map(i) = i_row
         endif
         if(mod((i-1)/eh%bh%blk,eh%bh%n_pcol) == eh%bh%my_pcol) then
            i_col         = i_col+1
            eh%col_map(i) = i_col
         endif
      enddo

      eh%bh%blacs_ready = .true.
   endif

end subroutine

!>
!! This routine sets the sparsity pattern.
!!
subroutine elsi_set_csc(eh,nnz_g,nnz_l,n_lcol,row_ind,col_ptr)

   implicit none

   type(elsi_handle), intent(inout) :: eh                !< Handle
   integer(kind=i4),  intent(in)    :: nnz_g             !< Global number of nonzeros
   integer(kind=i4),  intent(in)    :: nnz_l             !< Local number of nonzeros
   integer(kind=i4),  intent(in)    :: n_lcol            !< Local number of columns
   integer(kind=i4),  intent(in)    :: row_ind(nnz_l)    !< Row index
   integer(kind=i4),  intent(in)    :: col_ptr(n_lcol+1) !< Column pointer

   character(len=40), parameter :: caller = "elsi_set_csc"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%bh%nnz_g     = nnz_g
   eh%bh%nnz_l_sp  = nnz_l
   eh%bh%n_lcol_sp = n_lcol

   select case(eh%ph%matrix_format)
   case(PEXSI_CSC)
      eh%bh%nnz_l_sp1  = nnz_l
      eh%bh%n_lcol_sp1 = n_lcol

      if(allocated(eh%row_ind_sp1)) then
         call elsi_deallocate(eh%bh,eh%row_ind_sp1,"row_ind_sp1")
      endif
      if(allocated(eh%col_ptr_sp1)) then
         call elsi_deallocate(eh%bh,eh%col_ptr_sp1,"col_ptr_sp1")
      endif

      call elsi_allocate(eh%bh,eh%row_ind_sp1,nnz_l,"row_ind_sp1",caller)
      call elsi_allocate(eh%bh,eh%col_ptr_sp1,n_lcol+1,"col_ptr_sp1",caller)

      eh%row_ind_sp1 = row_ind
      eh%col_ptr_sp1 = col_ptr

      eh%bh%pexsi_csc_ready = .true.
   case(SIESTA_CSC)
      eh%bh%nnz_l_sp2  = nnz_l
      eh%bh%n_lcol_sp2 = n_lcol

      if(allocated(eh%row_ind_sp2)) then
         call elsi_deallocate(eh%bh,eh%row_ind_sp2,"row_ind_sp2")
      endif
      if(allocated(eh%col_ptr_sp2)) then
         call elsi_deallocate(eh%bh,eh%col_ptr_sp2,"col_ptr_sp2")
      endif

      call elsi_allocate(eh%bh,eh%row_ind_sp2,nnz_l,"row_ind_sp2",caller)
      call elsi_allocate(eh%bh,eh%col_ptr_sp2,n_lcol+1,"col_ptr_sp2",caller)

      eh%row_ind_sp2 = row_ind
      eh%col_ptr_sp2 = col_ptr

      eh%bh%siesta_csc_ready = .true.
   end select

end subroutine

!>
!! This routine finalizes ELSI.
!!
subroutine elsi_finalize(eh)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle

   character(len=40), parameter :: caller = "elsi_finalize"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_final_print(eh%ph,eh%bh)
   call elsi_cleanup(eh)

end subroutine

!>
!! This routine frees memory.
!!
subroutine elsi_cleanup(eh)

   implicit none

   type(elsi_handle), intent(inout) :: eh

   character(len=40), parameter :: caller = "elsi_cleanup"

   call elsi_cleanup_dmp(eh%ph)
   call elsi_cleanup_elpa(eh%ph)
   call elsi_cleanup_omm(eh%ph)
   call elsi_cleanup_pexsi(eh%ph)
   call elsi_cleanup_sips(eh%ph)

   ! Dense arrays
   if(allocated(eh%ham_real_den)) then
      call elsi_deallocate(eh%bh,eh%ham_real_den,"ham_real_den")
   endif
   if(allocated(eh%ham_cmplx_den)) then
      call elsi_deallocate(eh%bh,eh%ham_cmplx_den,"ham_cmplx_den")
   endif
   if(allocated(eh%ovlp_real_den)) then
      call elsi_deallocate(eh%bh,eh%ovlp_real_den,"ovlp_real_den")
   endif
   if(allocated(eh%ovlp_cmplx_den)) then
      call elsi_deallocate(eh%bh,eh%ovlp_cmplx_den,"ovlp_cmplx_den")
   endif
   if(allocated(eh%eval)) then
      call elsi_deallocate(eh%bh,eh%eval,"eval")
   endif
   if(allocated(eh%evec_real)) then
      call elsi_deallocate(eh%bh,eh%evec_real,"evec_real")
   endif
   if(allocated(eh%evec_cmplx)) then
      call elsi_deallocate(eh%bh,eh%evec_cmplx,"evec_cmplx")
   endif
   if(allocated(eh%dm_real_den)) then
      call elsi_deallocate(eh%bh,eh%dm_real_den,"dm_real_den")
   endif
   if(allocated(eh%dm_cmplx_den)) then
      call elsi_deallocate(eh%bh,eh%dm_cmplx_den,"dm_cmplx_den")
   endif

   ! Sparse arrays
   if(allocated(eh%ham_real_csc)) then
      call elsi_deallocate(eh%bh,eh%ham_real_csc,"ham_real_csc")
   endif
   if(allocated(eh%ham_cmplx_csc)) then
      call elsi_deallocate(eh%bh,eh%ham_cmplx_csc,"ham_cmplx_csc")
   endif
   if(allocated(eh%ovlp_real_csc)) then
      call elsi_deallocate(eh%bh,eh%ovlp_real_csc,"ovlp_real_csc")
   endif
   if(allocated(eh%ovlp_cmplx_csc)) then
      call elsi_deallocate(eh%bh,eh%ovlp_cmplx_csc,"ovlp_cmplx_csc")
   endif
   if(allocated(eh%dm_real_csc)) then
      call elsi_deallocate(eh%bh,eh%dm_real_csc,"dm_real_csc")
   endif
   if(allocated(eh%dm_cmplx_csc)) then
      call elsi_deallocate(eh%bh,eh%dm_cmplx_csc,"dm_cmplx_csc")
   endif
   if(allocated(eh%row_ind_sp1)) then
      call elsi_deallocate(eh%bh,eh%row_ind_sp1,"row_ind_sp1")
   endif
   if(allocated(eh%col_ptr_sp1)) then
      call elsi_deallocate(eh%bh,eh%col_ptr_sp1,"col_ptr_sp1")
   endif
   if(allocated(eh%row_ind_sp2)) then
      call elsi_deallocate(eh%bh,eh%row_ind_sp2,"row_ind_sp2")
   endif
   if(allocated(eh%col_ptr_sp2)) then
      call elsi_deallocate(eh%bh,eh%col_ptr_sp2,"col_ptr_sp2")
   endif

   ! Auxiliary arrays
   if(allocated(eh%ham_real_copy)) then
      call elsi_deallocate(eh%bh,eh%ham_real_copy,"ham_real_copy")
   endif
   if(allocated(eh%ovlp_real_copy)) then
      call elsi_deallocate(eh%bh,eh%ovlp_real_copy,"ovlp_real_copy")
   endif
   if(allocated(eh%ovlp_cmplx_copy)) then
      call elsi_deallocate(eh%bh,eh%ovlp_cmplx_copy,"ovlp_cmplx_copy")
   endif
   if(allocated(eh%ovlp_real_inv)) then
      call elsi_deallocate(eh%bh,eh%ovlp_real_inv,"ovlp_real_inv")
   endif
   if(allocated(eh%occ)) then
      call elsi_deallocate(eh%bh,eh%occ,"occ")
   endif
   if(allocated(eh%row_map)) then
      call elsi_deallocate(eh%bh,eh%row_map,"row_map")
   endif
   if(allocated(eh%col_map)) then
      call elsi_deallocate(eh%bh,eh%col_map,"col_map")
   endif
   if(allocated(eh%omm_c_real)) then
      call elsi_deallocate(eh%bh,eh%omm_c_real,"omm_c_real")
   endif
   if(allocated(eh%omm_c_cmplx)) then
      call elsi_deallocate(eh%bh,eh%omm_c_cmplx,"omm_c_cmplx")
   endif
   if(allocated(eh%pexsi_ne_vec)) then
      call elsi_deallocate(eh%bh,eh%pexsi_ne_vec,"pexsi_ne_vec")
   endif
   if(allocated(eh%dmp_vec1)) then
      call elsi_deallocate(eh%bh,eh%dmp_vec1,"dmp_vec1")
   endif
   if(allocated(eh%dmp_vec2)) then
      call elsi_deallocate(eh%bh,eh%dmp_vec2,"dmp_vec2")
   endif

   if(eh%bh%json_init) then
      call fjson_finish_array(eh%jh)
      call fjson_close_file(eh%jh)
      call fjson_reset_fj_handle(eh%jh)
   endif

   ! Reset handle
   call elsi_reset_param(eh%ph)
   call elsi_reset_basic(eh%bh)

   eh%handle_init = .false.

end subroutine

end module ELSI_SETUP
