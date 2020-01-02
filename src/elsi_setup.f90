! Copyright (c) 2015-2020, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide routines to set up an ELSI instance.
!!
module ELSI_SETUP

   use ELSI_CONSTANT, only: AUTO_SOLVER,ELPA_SOLVER,PEXSI_SOLVER,SINGLE_PROC,&
       MULTI_PROC,PEXSI_CSC,SIESTA_CSC,UNSET,DECISION_INIT
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_EIGENEXA, only: elsi_cleanup_eigenexa
   use ELSI_ELPA, only: elsi_cleanup_elpa
   use ELSI_MAGMA, only: elsi_cleanup_magma
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_NTPOLY, only: elsi_cleanup_ntpoly
   use ELSI_OMM, only: elsi_cleanup_omm
   use ELSI_OUTPUT, only: fjson_close_file,fjson_finish_array,&
       fjson_reset_fj_handle
   use ELSI_PEXSI, only: elsi_set_pexsi_default,elsi_cleanup_pexsi
   use ELSI_PRECISION, only: r8,i4,i8
   use ELSI_SIPS, only: elsi_cleanup_sips
   use ELSI_SORT, only: elsi_heapsort
   use ELSI_UTIL, only: elsi_check_init,elsi_reset_param,elsi_reset_basic

   implicit none

   private

   public :: elsi_init
   public :: elsi_set_mpi
   public :: elsi_set_mpi_global
   public :: elsi_set_spin
   public :: elsi_set_kpoint
   public :: elsi_set_blacs
   public :: elsi_set_csc
   public :: elsi_set_csc_blk
   public :: elsi_set_coo
   public :: elsi_reinit
   public :: elsi_finalize
   public :: elsi_cleanup

contains

!>
!! Initialize ELSI with user's choice of solver, parallel mode, matrix format,
!! number of basis functions (global size of Hamiltonian), number of electrons,
!! and number of states.
!!
subroutine elsi_init(eh,solver,parallel_mode,matrix_format,n_basis,n_electron,&
   n_state)

   implicit none

   type(elsi_handle), intent(out) :: eh !< Handle
   integer(kind=i4), intent(in) :: solver !< Solver
   integer(kind=i4), intent(in) :: parallel_mode !< Parallel mode
   integer(kind=i4), intent(in) :: matrix_format !< Matrix format
   integer(kind=i4), intent(in) :: n_basis !< Number of basis functions
   real(kind=r8), intent(in) :: n_electron !< Number of electrons
   integer(kind=i4), intent(in) :: n_state !< Number of states

   character(len=*), parameter :: caller = "elsi_init"

   ! For safety
   call elsi_cleanup(eh)

   eh%handle_init = .true.
   eh%ph%n_basis = n_basis
   eh%ph%n_good = n_basis
   eh%ph%n_electrons = n_electron
   eh%ph%n_states = min(n_basis,n_state)
   eh%ph%n_states_solve = min(n_basis,n_state)
   eh%ph%solver = solver
   eh%ph%matrix_format = matrix_format
   eh%ph%parallel_mode = parallel_mode

   if(parallel_mode == SINGLE_PROC) then
      eh%bh%n_lrow = n_basis
      eh%bh%n_lcol = n_basis
      eh%bh%blk = n_basis
      eh%bh%n_prow = 1
      eh%bh%n_pcol = 1
      eh%bh%myid = 0
      eh%bh%n_procs = 1
      eh%bh%myid_all = 0
      eh%bh%n_procs_all = 1
   end if

   if(solver == PEXSI_SOLVER .or. solver == AUTO_SOLVER) then
      ! This overrides user settings, so must call it here
      call elsi_set_pexsi_default(eh%ph)
   end if

   if(solver == AUTO_SOLVER) then
      eh%ph%decision_stage = DECISION_INIT
      eh%bh%print_json = 1
   end if

end subroutine

!>
!! Set the MPI communicator to be used to solve one eigenproblem (one spin
!! channel, one k-point).
!!
subroutine elsi_set_mpi(eh,comm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: comm !< Unit communicator

   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_set_mpi"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(eh%ph%parallel_mode == MULTI_PROC) then
      eh%bh%comm = comm

      call MPI_Comm_rank(comm,eh%bh%myid,ierr)
      call MPI_Comm_size(comm,eh%bh%n_procs,ierr)

      eh%bh%mpi_ready = .true.
   end if

end subroutine

!>
!! Set the global MPI communicator to be used to exchange information across all
!! spin channels and k-points.
!!
subroutine elsi_set_mpi_global(eh,comm_all)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: comm_all !< Global communicator

   integer(kind=i4) :: ierr

   character(len=*), parameter :: caller = "elsi_set_mpi_global"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(eh%ph%parallel_mode == MULTI_PROC) then
      eh%bh%comm_all = comm_all

      call MPI_Comm_rank(comm_all,eh%bh%myid_all,ierr)
      call MPI_Comm_size(comm_all,eh%bh%n_procs_all,ierr)

      eh%bh%mpi_all_ready = .true.
   end if

end subroutine

!>
!! Set the number of spin channels and the index of the spin channel the calling
!! process is solving.
!!
subroutine elsi_set_spin(eh,n_spin,i_spin)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_spin !< Number of spin channels
   integer(kind=i4), intent(in) :: i_spin !< Spin index

   eh%ph%n_spins = n_spin
   eh%ph%i_spin = i_spin

end subroutine

!>
!! Set the number of k-points, and the index and weight of the k-point the
!! calling process is solving.
!!
subroutine elsi_set_kpoint(eh,n_kpt,i_kpt,i_wt)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_kpt !< Number of k-points
   integer(kind=i4), intent(in) :: i_kpt !< K-point index
   real(kind=r8), intent(in) :: i_wt !< Weight

   eh%ph%n_kpts = n_kpt
   eh%ph%i_kpt = i_kpt
   eh%ph%i_wt = i_wt

end subroutine

!>
!! Set the BLACS context and block size used by the BLACS_DENSE matrix format.
!!
subroutine elsi_set_blacs(eh,blacs_ctxt,block_size)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: blacs_ctxt !< BLACS context
   integer(kind=i4), intent(in) :: block_size !< Block size

   integer(kind=i4) :: ierr

   integer(kind=i4), external :: numroc

   character(len=*), parameter :: caller = "elsi_set_blacs"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(eh%ph%parallel_mode == MULTI_PROC) then
      eh%bh%blacs_ctxt = blacs_ctxt
      eh%bh%blk = block_size

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

      eh%bh%blacs_ready = .true.
   end if

end subroutine

!>
!! Set the global number of nonzero matrix elements, local number of nonzero
!! matrix elements, local number of matrix columns, row index array, and column
!! pointer array. These variables are collectively referred to as the CSC
!! sparsity pattern, used by the PEXSI_CSC and SIESTA_CSC matrix formats.
!!
subroutine elsi_set_csc(eh,nnz_g,nnz_l,n_lcol,row_ind,col_ptr)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: nnz_g !< Global number of nonzeros
   integer(kind=i4), intent(in) :: nnz_l !< Local number of nonzeros
   integer(kind=i4), intent(in) :: n_lcol !< Local number of columns
   integer(kind=i4), intent(in) :: row_ind(nnz_l) !< Row index
   integer(kind=i4), intent(in) :: col_ptr(n_lcol+1) !< Column pointer

   character(len=*), parameter :: caller = "elsi_set_csc"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%bh%nnz_g = nnz_g
   eh%bh%nnz_l_sp = nnz_l
   eh%bh%n_lcol_sp = n_lcol

   select case(eh%ph%matrix_format)
   case(PEXSI_CSC)
      eh%bh%nnz_l_sp1 = nnz_l
      eh%bh%n_lcol_sp1 = n_lcol

      if(allocated(eh%row_ind_sp1)) then
         call elsi_deallocate(eh%bh,eh%row_ind_sp1,"row_ind_sp1")
      end if

      if(allocated(eh%col_ptr_sp1)) then
         call elsi_deallocate(eh%bh,eh%col_ptr_sp1,"col_ptr_sp1")
      end if

      call elsi_allocate(eh%bh,eh%row_ind_sp1,nnz_l,"row_ind_sp1",caller)
      call elsi_allocate(eh%bh,eh%col_ptr_sp1,n_lcol+1,"col_ptr_sp1",caller)

      eh%row_ind_sp1 = row_ind
      eh%col_ptr_sp1 = col_ptr

      eh%bh%pexsi_csc_ready = .true.
   case(SIESTA_CSC)
      eh%bh%nnz_l_sp2 = nnz_l
      eh%bh%n_lcol_sp2 = n_lcol

      if(allocated(eh%row_ind_sp2)) then
         call elsi_deallocate(eh%bh,eh%row_ind_sp2,"row_ind_sp2")
      end if

      if(allocated(eh%col_ptr_sp2)) then
         call elsi_deallocate(eh%bh,eh%col_ptr_sp2,"col_ptr_sp2")
      end if

      call elsi_allocate(eh%bh,eh%row_ind_sp2,nnz_l,"row_ind_sp2",caller)
      call elsi_allocate(eh%bh,eh%col_ptr_sp2,n_lcol+1,"col_ptr_sp2",caller)

      eh%row_ind_sp2 = row_ind
      eh%col_ptr_sp2 = col_ptr

      eh%bh%siesta_csc_ready = .true.
   end select

end subroutine

!>
!! Set the block size in 1D block-cyclic distributed CSC format.
!!
subroutine elsi_set_csc_blk(eh,blk)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: blk !< Block size

   character(len=*), parameter :: caller = "elsi_set_csc_blk"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%bh%blk_sp2 = blk

end subroutine

!>
!! Set the global number of nonzero matrix elements, local number of nonzero
!! matrix elements, row index array, and column index array. These variables are
!! collectively referred to as COO sparsity pattern, used by the GENERIC_COO
!! matrix format.
!!
subroutine elsi_set_coo(eh,nnz_g,nnz_l,row_ind,col_ind)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   integer(kind=i4), intent(in) :: nnz_g !< Global number of nonzeros
   integer(kind=i4), intent(in) :: nnz_l !< Local number of nonzeros
   integer(kind=i4), intent(in) :: row_ind(nnz_l) !< Row index
   integer(kind=i4), intent(in) :: col_ind(nnz_l) !< Column index

   integer(kind=i8), allocatable :: gid(:) ! Global 1D id

   character(len=*), parameter :: caller = "elsi_set_coo"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   eh%bh%nnz_g = nnz_g
   eh%bh%nnz_l_sp = nnz_l
   eh%bh%nnz_l_sp3 = nnz_l

   if(allocated(eh%row_ind_sp3)) then
      call elsi_deallocate(eh%bh,eh%row_ind_sp3,"row_ind_sp3")
   end if

   if(allocated(eh%col_ind_sp3)) then
      call elsi_deallocate(eh%bh,eh%col_ind_sp3,"col_ind_sp3")
   end if

   if(allocated(eh%perm_sp3)) then
      call elsi_deallocate(eh%bh,eh%perm_sp3,"perm_sp3")
   end if

   call elsi_allocate(eh%bh,eh%row_ind_sp3,nnz_l,"row_ind_sp3",caller)
   call elsi_allocate(eh%bh,eh%col_ind_sp3,nnz_l,"col_ind_sp3",caller)
   call elsi_allocate(eh%bh,eh%perm_sp3,nnz_l,"perm_sp3",caller)
   call elsi_allocate(eh%bh,gid,nnz_l,"gid",caller)

   eh%row_ind_sp3 = row_ind
   eh%col_ind_sp3 = col_ind

   ! Compute global 1D id
   gid = int(col_ind-1,kind=i8)*int(eh%ph%n_basis,kind=i8)+int(row_ind,kind=i8)

   ! Sort
   call elsi_heapsort(nnz_l,gid,eh%perm_sp3)

   call elsi_deallocate(eh%bh,gid,"gid")

   eh%bh%generic_coo_ready = .true.

end subroutine

!>
!! Start a new geometry step, which usually means a new overlap.
!!
subroutine elsi_reinit(eh)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle

   character(len=*), parameter :: caller = "elsi_reinit"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(eh%ph%n_calls > 0) then
      eh%bh%nnz_l = UNSET
      eh%bh%nnz_g = UNSET
      eh%bh%nnz_l_sp = UNSET
      eh%bh%nnz_l_sp1 = UNSET
      eh%bh%pexsi_csc_ready = .false.
      eh%bh%nnz_l_sp2 = UNSET
      eh%bh%siesta_csc_ready = .false.
      eh%ph%n_calls_all = eh%ph%n_calls_all+eh%ph%n_calls
      eh%ph%n_calls = 0
      eh%ph%ill_ovlp = .false.
      eh%ph%n_good = eh%ph%n_basis
      eh%ph%first_blacs_to_ntpoly = .true.
      eh%ph%first_blacs_to_pexsi = .true.
      eh%ph%first_blacs_to_sips = .true.
      eh%ph%first_generic_to_blacs = .true.
      eh%ph%first_generic_to_ntpoly = .true.
      eh%ph%first_generic_to_pexsi = .true.
      eh%ph%first_siesta_to_blacs = .true.
      eh%ph%first_siesta_to_ntpoly = .true.
      eh%ph%first_siesta_to_pexsi = .true.
      eh%ph%first_sips_to_blacs = .true.
      eh%ph%first_sips_to_ntpoly = .true.
      eh%ph%elpa_first = .true.
      eh%ph%omm_first = .true.
      eh%ph%pexsi_first = .true.
      eh%ph%sips_first = .true.
      eh%ph%nt_first = .true.
      eh%ph%exa_first = .true.

      call elsi_cleanup_pexsi(eh%ph)
      call elsi_cleanup_sips(eh%ph)

      if(allocated(eh%evec_real)) then
         call elsi_deallocate(eh%bh,eh%evec_real,"evec_real")
      end if

      if(allocated(eh%evec_cmplx)) then
         call elsi_deallocate(eh%bh,eh%evec_cmplx,"evec_cmplx")
      end if

      if(allocated(eh%ham_real_sp)) then
         call elsi_deallocate(eh%bh,eh%ham_real_sp,"ham_real_sp")
      end if

      if(allocated(eh%ham_cmplx_sp)) then
         call elsi_deallocate(eh%bh,eh%ham_cmplx_sp,"ham_cmplx_sp")
      end if

      if(allocated(eh%ovlp_real_sp)) then
         call elsi_deallocate(eh%bh,eh%ovlp_real_sp,"ovlp_real_sp")
      end if

      if(allocated(eh%ovlp_cmplx_sp)) then
         call elsi_deallocate(eh%bh,eh%ovlp_cmplx_sp,"ovlp_cmplx_sp")
      end if

      if(allocated(eh%dm_real_sp)) then
         call elsi_deallocate(eh%bh,eh%dm_real_sp,"dm_real_sp")
      end if

      if(allocated(eh%dm_cmplx_sp)) then
         call elsi_deallocate(eh%bh,eh%dm_cmplx_sp,"dm_cmplx_sp")
      end if

      if(allocated(eh%row_ind_sp1)) then
         call elsi_deallocate(eh%bh,eh%row_ind_sp1,"row_ind_sp1")
      end if

      if(allocated(eh%col_ptr_sp1)) then
         call elsi_deallocate(eh%bh,eh%col_ptr_sp1,"col_ptr_sp1")
      end if

      if(allocated(eh%row_ind_sp2)) then
         call elsi_deallocate(eh%bh,eh%row_ind_sp2,"row_ind_sp2")
      end if

      if(allocated(eh%col_ptr_sp2)) then
         call elsi_deallocate(eh%bh,eh%col_ptr_sp2,"col_ptr_sp2")
      end if

      if(allocated(eh%row_ind_sp3)) then
         call elsi_deallocate(eh%bh,eh%row_ind_sp3,"row_ind_sp3")
      end if

      if(allocated(eh%col_ind_sp3)) then
         call elsi_deallocate(eh%bh,eh%col_ind_sp3,"col_ind_sp3")
      end if

      if(allocated(eh%map_den)) then
         call elsi_deallocate(eh%bh,eh%map_den,"map_den")
      end if

      if(allocated(eh%map_sp1)) then
         call elsi_deallocate(eh%bh,eh%map_sp1,"map_sp1")
      end if

      if(allocated(eh%perm_sp3)) then
         call elsi_deallocate(eh%bh,eh%perm_sp3,"perm_sp3")
      end if

      if(.not. eh%ph%solver == ELPA_SOLVER) then
         if(allocated(eh%ovlp_real_copy)) then
            call elsi_deallocate(eh%bh,eh%ovlp_real_copy,"ovlp_real_copy")
         end if

         if(allocated(eh%ovlp_cmplx_copy)) then
            call elsi_deallocate(eh%bh,eh%ovlp_cmplx_copy,"ovlp_cmplx_copy")
         end if
      end if
   end if

end subroutine

!>
!! Finalize ELSI.
!!
subroutine elsi_finalize(eh)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle

   character(len=*), parameter :: caller = "elsi_finalize"

   eh%ph%n_calls_all = eh%ph%n_calls_all+eh%ph%n_calls

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_cleanup(eh)

end subroutine

!>
!! Free memory.
!!
subroutine elsi_cleanup(eh)

   implicit none

   type(elsi_handle), intent(inout) :: eh

   character(len=*), parameter :: caller = "elsi_cleanup"

   call elsi_cleanup_elpa(eh%ph)
   call elsi_cleanup_omm(eh%ph)
   call elsi_cleanup_pexsi(eh%ph)
   call elsi_cleanup_eigenexa(eh%ph)
   call elsi_cleanup_ntpoly(eh%ph)
   call elsi_cleanup_sips(eh%ph)
   call elsi_cleanup_magma(eh%ph)

   ! Dense arrays
   if(allocated(eh%ham_real_den)) then
      call elsi_deallocate(eh%bh,eh%ham_real_den,"ham_real_den")
   end if

   if(allocated(eh%ham_cmplx_den)) then
      call elsi_deallocate(eh%bh,eh%ham_cmplx_den,"ham_cmplx_den")
   end if

   if(allocated(eh%ovlp_real_den)) then
      call elsi_deallocate(eh%bh,eh%ovlp_real_den,"ovlp_real_den")
   end if

   if(allocated(eh%ovlp_cmplx_den)) then
      call elsi_deallocate(eh%bh,eh%ovlp_cmplx_den,"ovlp_cmplx_den")
   end if

   if(allocated(eh%eval)) then
      call elsi_deallocate(eh%bh,eh%eval,"eval")
   end if

   if(allocated(eh%evec_real)) then
      call elsi_deallocate(eh%bh,eh%evec_real,"evec_real")
   end if

   if(allocated(eh%evec_cmplx)) then
      call elsi_deallocate(eh%bh,eh%evec_cmplx,"evec_cmplx")
   end if

   if(allocated(eh%dm_real_den)) then
      call elsi_deallocate(eh%bh,eh%dm_real_den,"dm_real_den")
   end if

   if(allocated(eh%dm_cmplx_den)) then
      call elsi_deallocate(eh%bh,eh%dm_cmplx_den,"dm_cmplx_den")
   end if

   ! Sparse arrays
   if(allocated(eh%ham_real_sp)) then
      call elsi_deallocate(eh%bh,eh%ham_real_sp,"ham_real_sp")
   end if

   if(allocated(eh%ham_cmplx_sp)) then
      call elsi_deallocate(eh%bh,eh%ham_cmplx_sp,"ham_cmplx_sp")
   end if

   if(allocated(eh%ovlp_real_sp)) then
      call elsi_deallocate(eh%bh,eh%ovlp_real_sp,"ovlp_real_sp")
   end if

   if(allocated(eh%ovlp_cmplx_sp)) then
      call elsi_deallocate(eh%bh,eh%ovlp_cmplx_sp,"ovlp_cmplx_sp")
   end if

   if(allocated(eh%dm_real_sp)) then
      call elsi_deallocate(eh%bh,eh%dm_real_sp,"dm_real_sp")
   end if

   if(allocated(eh%dm_cmplx_sp)) then
      call elsi_deallocate(eh%bh,eh%dm_cmplx_sp,"dm_cmplx_sp")
   end if

   if(allocated(eh%row_ind_sp1)) then
      call elsi_deallocate(eh%bh,eh%row_ind_sp1,"row_ind_sp1")
   end if

   if(allocated(eh%col_ptr_sp1)) then
      call elsi_deallocate(eh%bh,eh%col_ptr_sp1,"col_ptr_sp1")
   end if

   if(allocated(eh%row_ind_sp2)) then
      call elsi_deallocate(eh%bh,eh%row_ind_sp2,"row_ind_sp2")
   end if

   if(allocated(eh%col_ptr_sp2)) then
      call elsi_deallocate(eh%bh,eh%col_ptr_sp2,"col_ptr_sp2")
   end if

   if(allocated(eh%row_ind_sp3)) then
      call elsi_deallocate(eh%bh,eh%row_ind_sp3,"row_ind_sp3")
   end if

   if(allocated(eh%col_ind_sp3)) then
      call elsi_deallocate(eh%bh,eh%col_ind_sp3,"col_ind_sp3")
   end if

   ! Auxiliary arrays
   if(allocated(eh%ovlp_real_copy)) then
      call elsi_deallocate(eh%bh,eh%ovlp_real_copy,"ovlp_real_copy")
   end if

   if(allocated(eh%ovlp_cmplx_copy)) then
      call elsi_deallocate(eh%bh,eh%ovlp_cmplx_copy,"ovlp_cmplx_copy")
   end if

   if(allocated(eh%occ)) then
      call elsi_deallocate(eh%bh,eh%occ,"occ")
   end if

   if(allocated(eh%omm_c_real)) then
      call elsi_deallocate(eh%bh,eh%omm_c_real,"omm_c_real")
   end if

   if(allocated(eh%omm_c_cmplx)) then
      call elsi_deallocate(eh%bh,eh%omm_c_cmplx,"omm_c_cmplx")
   end if

   if(allocated(eh%pexsi_ne_vec)) then
      call elsi_deallocate(eh%bh,eh%pexsi_ne_vec,"pexsi_ne_vec")
   end if

   if(allocated(eh%map_den)) then
      call elsi_deallocate(eh%bh,eh%map_den,"map_den")
   end if

   if(allocated(eh%map_sp1)) then
      call elsi_deallocate(eh%bh,eh%map_sp1,"map_sp1")
   end if

   if(allocated(eh%perm_sp3)) then
      call elsi_deallocate(eh%bh,eh%perm_sp3,"perm_sp3")
   end if

   if(eh%bh%json_init) then
      call fjson_finish_array(eh%jh)
      call fjson_close_file(eh%jh)
      call fjson_reset_fj_handle(eh%jh)
   end if

   ! Reset handle
   call elsi_reset_param(eh%ph)
   call elsi_reset_basic(eh%bh)

   eh%handle_init = .false.

end subroutine

end module ELSI_SETUP
