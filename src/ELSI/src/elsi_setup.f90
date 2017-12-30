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
!! This module provides routines for setting up an ELSI instance.
!!
module ELSI_SETUP

   use ELSI_CHESS,         only: elsi_set_chess_default
   use ELSI_CONSTANTS,     only: ELPA_SOLVER,OMM_SOLVER,PEXSI_SOLVER,&
                                 CHESS_SOLVER,SIPS_SOLVER,SINGLE_PROC,&
                                 MULTI_PROC,SOLVER_TIMINGS_UNIT_DEFAULT,&
                                 SOLVER_TIMINGS_FILE_DEFAULT
   use ELSI_DATATYPE
   use ELSI_DMP,           only: elsi_set_dmp_default
   use ELSI_ELPA,          only: elsi_set_elpa_default,elsi_get_elpa_comms
   use ELSI_IO,            only: elsi_print_handle_summary,elsi_say,&
                                 elsi_say_setting,elsi_init_file_io,&
                                 elsi_reset_file_io_handle,append_string,&
                                 truncate_string
   use ELSI_MALLOC
   use ELSI_OMM,           only: elsi_set_omm_default
   use ELSI_PEXSI,         only: elsi_set_pexsi_default
   use ELSI_PRECISION,     only: r8,i4
   use ELSI_SIPS,          only: elsi_set_sips_default
   use ELSI_TIMINGS,       only: elsi_init_timer,elsi_init_timings,&
                                 elsi_print_timings,elsi_finalize_timings
   use ELSI_UTILS
   use FOE_BASE,           only: foe_data_deallocate
   use F_PPEXSI_INTERFACE, only: f_ppexsi_plan_finalize
   use MATRIXSWITCH,       only: ms_scalapack_setup,m_deallocate
   use M_QETSC,            only: clean_qetsc
   use SPARSEMATRIX_BASE,  only: deallocate_sparse_matrix,deallocate_matrices

   implicit none

   private

   public  :: elsi_init
   public  :: elsi_finalize
   public  :: elsi_set_mpi
   public  :: elsi_set_mpi_global
   public  :: elsi_set_spin
   public  :: elsi_set_kpoint
   public  :: elsi_set_blacs
   public  :: elsi_set_csc
   public  :: elsi_cleanup

   private :: elsi_final_print

contains

!>
!! This routine initializes ELSI with the solver, parallel mode, matrix storage
!! format, number of basis functions (global size of the Hamiltonian matrix),
!! number of electrons, and number of states.
!!
subroutine elsi_init(e_h,solver,parallel_mode,matrix_format,n_basis,n_electron,&
              n_state)

   implicit none

   type(elsi_handle), intent(out) :: e_h           !< Handle
   integer(kind=i4),  intent(in)  :: solver        !< AUTO,ELPA,LIBOMM,PEXSI,CHESS,SIPS,DMP
   integer(kind=i4),  intent(in)  :: parallel_mode !< SINGLE_PROC,MULTI_PROC
   integer(kind=i4),  intent(in)  :: matrix_format !< BLACS_DENSE,PEXSI_CSC
   integer(kind=i4),  intent(in)  :: n_basis       !< Number of basis functions
   real(kind=r8),     intent(in)  :: n_electron    !< Number of electrons
   integer(kind=i4),  intent(in)  :: n_state       !< Number of states

   character*40, parameter :: caller = "elsi_init"

   ! For safety
   call elsi_cleanup(e_h)

   e_h%handle_init    = .true.
   e_h%n_basis        = n_basis
   e_h%n_nonsing      = n_basis
   e_h%n_electrons    = n_electron
   e_h%n_states       = n_state
   e_h%n_states_solve = n_state
   e_h%n_states_omm   = nint(n_electron/2.0_r8)
   e_h%n_states_dmp   = nint(n_electron/2.0_r8)
   e_h%solver         = solver
   e_h%matrix_format  = matrix_format
   e_h%parallel_mode  = parallel_mode
   e_h%n_elsi_calls   = 0


   if(parallel_mode == SINGLE_PROC) then
      e_h%n_lrow      = n_basis
      e_h%n_lcol      = n_basis
      e_h%blk_row     = n_basis
      e_h%blk_col     = n_basis
      e_h%myid        = 0
      e_h%n_procs     = 1
      e_h%myid_all    = 0
      e_h%n_procs_all = 1
   endif

   ! Set default paramters
   select case(solver)
   case(ELPA_SOLVER)
      call elsi_set_elpa_default(e_h)
   case(OMM_SOLVER)
      call elsi_set_elpa_default(e_h)
      call elsi_set_omm_default(e_h)
   case(PEXSI_SOLVER)
      call elsi_set_pexsi_default(e_h)
   case(CHESS_SOLVER)
      call elsi_set_chess_default(e_h)
   case(SIPS_SOLVER)
      call elsi_set_elpa_default(e_h)
      call elsi_set_sips_default(e_h)
   case(DMP_SOLVER)
      call elsi_set_dmp_default(e_h)
   end select

   ! Initialize stdio handle
   ! By default, ELSI is silent to stdio unless user requests output
   call elsi_init_file_io(e_h%stdio,6,format=HUMAN_READ,print_info=.false.)

   ! Initialize solver timings file handle
   e_h%output_solver_timings = .true.
   ! print_unit and file_name should only have meaning for myid_all.eq.0
   ! as a precautionary measure
   ! However, we don't know which process has myid_all.eq.0 yet,
   ! so we'll unset them once we actually open the file
   call elsi_init_file_io(e_h%solver_timings_file,SOLVER_TIMINGS_UNIT_DEFAULT,&
                          file_name=SOLVER_TIMINGS_FILE_DEFAULT,&
                          format=JSON,print_info=.true.,&
                          comma_json=COMMA_AFTER)

   ! Initialize timer information
   call elsi_init_timer(e_h)
   call elsi_init_timings(e_h%solver_timings, "Solver timings")

end subroutine

!>
!! This routine sets the MPI communicator.
!!
subroutine elsi_set_mpi(e_h,mpi_comm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h      !< Handle
   integer(kind=i4),  intent(in)    :: mpi_comm !< Unit ELSI communicator

   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_set_mpi"

   call elsi_check_handle(e_h,caller)

   if(e_h%parallel_mode == MULTI_PROC) then
      ! Unit ELSI communicator
      ! If there's more than one spin/kpt, each unit communicator
      ! solves one KS problem
      e_h%mpi_comm = mpi_comm

      call MPI_Comm_rank(mpi_comm,e_h%myid,mpierr)
      call MPI_Comm_size(mpi_comm,e_h%n_procs,mpierr)

      e_h%mpi_ready = .true.

      if (e_h%handle_ready) e_h%handle_changed = .true.
   endif

end subroutine

!>
!! This routine sets the global MPI communicator.
!!
subroutine elsi_set_mpi_global(e_h,mpi_comm_all)

   implicit none

   type(elsi_handle), intent(inout) :: e_h          !< Handle
   integer(kind=i4),  intent(in)    :: mpi_comm_all !< Unit ELSI communicator

   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_set_mpi_global"

   call elsi_check_handle(e_h,caller)

   if(e_h%parallel_mode == MULTI_PROC) then
      ! Global ELSI communicator
      e_h%mpi_comm_all = mpi_comm_all

      call MPI_Comm_rank(mpi_comm_all,e_h%myid_all,mpierr)
      call MPI_Comm_size(mpi_comm_all,e_h%n_procs_all,mpierr)

      e_h%global_mpi_ready = .true.

      if (e_h%handle_ready) e_h%handle_changed = .true.
   endif

end subroutine

!>
!! This routine sets the spin information.
!!
subroutine elsi_set_spin(e_h,n_spin,i_spin)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   integer(kind=i4),  intent(in)    :: n_spin !< Number of spin channels
   integer(kind=i4),  intent(in)    :: i_spin !< Spin index

   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%n_spins = n_spin
   e_h%i_spin  = i_spin

end subroutine

!>
!! This routine sets the k-point information.
!!
subroutine elsi_set_kpoint(e_h,n_kpt,i_kpt,weight)

   implicit none

   type(elsi_handle), intent(inout) :: e_h    !< Handle
   integer(kind=i4),  intent(in)    :: n_kpt  !< Number of k-points
   integer(kind=i4),  intent(in)    :: i_kpt  !< K-point index
   real(kind=r8),     intent(in)    :: weight !< Weight

   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%n_kpts   = n_kpt
   e_h%i_kpt    = i_kpt
   e_h%i_weight = weight

end subroutine

!>
!! This routine sets the BLACS context and the block size.
!!
subroutine elsi_set_blacs(e_h,blacs_ctxt,block_size)

   implicit none

   type(elsi_handle), intent(inout) :: e_h        !< Handle
   integer(kind=i4),  intent(in)    :: blacs_ctxt !< BLACS context
   integer(kind=i4),  intent(in)    :: block_size !< Block size

   integer(kind=i4) :: i
   integer(kind=i4) :: i_row
   integer(kind=i4) :: i_col
   integer(kind=i4) :: ierr

   integer(kind=i4), external :: numroc

   character*40, parameter :: caller = "elsi_set_blacs"

   call elsi_check_handle(e_h,caller)

   if(e_h%parallel_mode == MULTI_PROC) then
      e_h%blacs_ctxt = blacs_ctxt
      e_h%blk_row    = block_size
      e_h%blk_col    = block_size

      ! Get processor grid information
      call blacs_gridinfo(e_h%blacs_ctxt,e_h%n_prow,e_h%n_pcol,e_h%my_prow,&
              e_h%my_pcol)

      ! Get local size of matrix
      e_h%n_lrow = numroc(e_h%n_basis,e_h%blk_row,e_h%my_prow,0,e_h%n_prow)
      e_h%n_lcol = numroc(e_h%n_basis,e_h%blk_col,e_h%my_pcol,0,e_h%n_pcol)

      ! Get BLACS descriptor
      call descinit(e_h%sc_desc,e_h%n_basis,e_h%n_basis,e_h%blk_row,&
              e_h%blk_col,0,0,e_h%blacs_ctxt,max(1,e_h%n_lrow),ierr)

      ! Get ELPA communicators
      call elsi_get_elpa_comms(e_h)

      ! Create global-local mapping
      call elsi_allocate(e_h,e_h%loc_row,e_h%n_basis,"loc_row",caller)
      call elsi_allocate(e_h,e_h%loc_col,e_h%n_basis,"loc_col",caller)

      i_row = 0
      i_col = 0

      do i = 1,e_h%n_basis
         if(mod((i-1)/e_h%blk_row,e_h%n_prow) == e_h%my_prow) then
            i_row = i_row+1
            e_h%loc_row(i) = i_row
         endif
         if(mod((i-1)/e_h%blk_col,e_h%n_pcol) == e_h%my_pcol) then
            i_col = i_col+1
            e_h%loc_col(i) = i_col
         endif
      enddo

      ! Set up MatrixSwitch
      if(e_h%solver == OMM_SOLVER) then
         call ms_scalapack_setup(e_h%mpi_comm,e_h%n_prow,'r',e_h%blk_row,&
                 icontxt=e_h%blacs_ctxt)
      endif

      e_h%blacs_ready = .true.
   
      if (e_h%handle_ready) e_h%handle_changed = .true.
   endif

end subroutine

!>
!! This routine sets the sparsity pattern.
!!
subroutine elsi_set_csc(e_h,nnz_g,nnz_l,n_lcol,row_ind,col_ptr)

   implicit none

   type(elsi_handle), intent(inout) :: e_h               !< Handle
   integer(kind=i4),  intent(in)    :: nnz_g             !< Global number of nonzeros
   integer(kind=i4),  intent(in)    :: nnz_l             !< Local number of nonzeros
   integer(kind=i4),  intent(in)    :: n_lcol            !< Local number of columns
   integer(kind=i4),  intent(in)    :: row_ind(nnz_l)    !< Row index
   integer(kind=i4),  intent(in)    :: col_ptr(n_lcol+1) !< Column pointer

   character*40, parameter :: caller = "elsi_set_csc"

   call elsi_check_handle(e_h,caller)
   if (e_h%handle_ready) e_h%handle_changed = .true.

   e_h%nnz_g     = nnz_g
   e_h%nnz_l_sp  = nnz_l
   e_h%n_lcol_sp = n_lcol

   if(e_h%solver == PEXSI_SOLVER) then
      if(allocated(e_h%row_ind_pexsi)) then
         call elsi_deallocate(e_h,e_h%row_ind_pexsi,"row_ind_pexsi")
      endif
      if(allocated(e_h%col_ptr_pexsi)) then
         call elsi_deallocate(e_h,e_h%col_ptr_pexsi,"col_ptr_pexsi")
      endif

      call elsi_allocate(e_h,e_h%row_ind_pexsi,nnz_l,"row_ind_pexsi",caller)
      call elsi_allocate(e_h,e_h%col_ptr_pexsi,n_lcol+1,"col_ptr_pexsi",caller)

      e_h%row_ind_pexsi = row_ind
      e_h%col_ptr_pexsi = col_ptr
   elseif(e_h%solver == CHESS_SOLVER) then
      if(allocated(e_h%row_ind_chess)) then
         call elsi_deallocate(e_h,e_h%row_ind_chess,"row_ind_chess")
      endif
      if(allocated(e_h%col_ptr_chess)) then
         call elsi_deallocate(e_h,e_h%col_ptr_chess,"col_ptr_chess")
      endif

      call elsi_allocate(e_h,e_h%row_ind_chess,nnz_l,"row_ind_chess",caller)
      call elsi_allocate(e_h,e_h%col_ptr_chess,n_lcol+1,"col_ptr_chess",caller)

      e_h%row_ind_chess = row_ind
      e_h%col_ptr_chess = col_ptr
   else
      if(allocated(e_h%row_ind_sips)) then
         call elsi_deallocate(e_h,e_h%row_ind_sips,"row_ind_sips")
      endif
      if(allocated(e_h%col_ptr_sips)) then
         call elsi_deallocate(e_h,e_h%col_ptr_sips,"col_ptr_sips")
      endif

      call elsi_allocate(e_h,e_h%row_ind_sips,nnz_l,"row_ind_sips",caller)
      call elsi_allocate(e_h,e_h%col_ptr_sips,nnz_l,"col_ptr_sips",caller)

      e_h%row_ind_sips = row_ind
      e_h%col_ptr_sips = col_ptr
   endif

   e_h%sparsity_ready = .true.

end subroutine

!>
!! This routine finalizes ELSI.
!!
subroutine elsi_finalize(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character*40, parameter :: caller = "elsi_finalize"

   call elsi_check_handle(e_h,caller)
   call elsi_final_print(e_h)
   call elsi_cleanup(e_h)

end subroutine

!>
!! This routine prints a final output.
!!
subroutine elsi_final_print(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   real(kind=r8) :: sparsity
   character*200 :: info_str

   character*40, parameter :: caller = "elsi_final_print"

   if(e_h%stdio%format.eq.JSON) then
      call elsi_stop("elsi_final_print only supports HUMAN_READ format.",&
                     e_h,caller)
   endif

   call elsi_say(e_h,"  |---------------------------------------------------------------------")
   call elsi_say(e_h,"  | Final ELSI Output                        ")
   call elsi_say(e_h,"  |---------------------------------------------------------------------")

   call append_string(e_h%stdio%prefix,"  | ")
   call elsi_print_handle_summary(e_h)

   call append_string(e_h%stdio%prefix,"  ")
   if(e_h%handle_changed) then
      call elsi_say_setting(e_h,"Was ELSI changed mid-run?","YES")
   else
      call elsi_say_setting(e_h,"Was ELSI changed mid-run?","NO")
   endif
   call elsi_say_setting(e_h,"Number of ELSI calls",e_h%n_elsi_calls)
   call truncate_string(e_h%stdio%prefix,2)

   call elsi_say(e_h,"")
   call elsi_say(e_h,"Timings")
   call append_string(e_h%stdio%prefix,"  ")
   call elsi_print_timings(e_h,e_h%solver_timings)
   call truncate_string(e_h%stdio%prefix,2)

   call truncate_string(e_h%stdio%prefix,4)
   call elsi_say(e_h,"  |--------------------------------------------------------------------")
   call elsi_say(e_h,"  | ELSI Project (c)  elsi-interchange.org   ")
   call elsi_say(e_h,"  |--------------------------------------------------------------------")

end subroutine

!>
!! This routine frees memory.
!!
subroutine elsi_cleanup(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: ierr

   character*40, parameter :: caller = "elsi_cleanup"

   ! ELPA
   if(allocated(e_h%ham_real_elpa)) then
      call elsi_deallocate(e_h,e_h%ham_real_elpa,"ham_real_elpa")
   endif
   if(allocated(e_h%ham_cmplx_elpa)) then
      call elsi_deallocate(e_h,e_h%ham_cmplx_elpa,"ham_cmplx_elpa")
   endif
   if(allocated(e_h%ovlp_real_elpa)) then
      call elsi_deallocate(e_h,e_h%ovlp_real_elpa,"ovlp_real_elpa")
   endif
   if(allocated(e_h%ovlp_cmplx_elpa)) then
      call elsi_deallocate(e_h,e_h%ovlp_cmplx_elpa,"ovlp_cmplx_elpa")
   endif
   if(allocated(e_h%eval_elpa)) then
      call elsi_deallocate(e_h,e_h%eval_elpa,"eval_elpa")
   endif
   if(allocated(e_h%evec_real_elpa)) then
      call elsi_deallocate(e_h,e_h%evec_real_elpa,"evec_real_elpa")
   endif
   if(allocated(e_h%evec_cmplx_elpa)) then
      call elsi_deallocate(e_h,e_h%evec_cmplx_elpa,"evec_cmplx_elpa")
   endif
   if(allocated(e_h%dm_real_elpa)) then
      call elsi_deallocate(e_h,e_h%dm_real_elpa,"dm_real_elpa")
   endif
   if(allocated(e_h%dm_cmplx_elpa)) then
      call elsi_deallocate(e_h,e_h%dm_cmplx_elpa,"dm_cmplx_elpa")
   endif
   if(allocated(e_h%occ_num)) then
      call elsi_deallocate(e_h,e_h%occ_num,"occ_num")
   endif
   if(allocated(e_h%k_weight)) then
      call elsi_deallocate(e_h,e_h%k_weight,"k_weight")
   endif
   if(allocated(e_h%eval_all)) then
      call elsi_deallocate(e_h,e_h%eval_all,"eval_all")
   endif

   ! libOMM
   if(e_h%ham_omm%is_initialized) then
      call m_deallocate(e_h%ham_omm)
   endif
   if(e_h%ovlp_omm%is_initialized) then
      call m_deallocate(e_h%ovlp_omm)
   endif
   if(e_h%dm_omm%is_initialized) then
      call m_deallocate(e_h%dm_omm)
   endif
   if(e_h%coeff%is_initialized) then
      call m_deallocate(e_h%coeff)
   endif
   if(e_h%tdm_omm%is_initialized) then
      call m_deallocate(e_h%tdm_omm)
   endif

   ! PEXSI
   if(allocated(e_h%ham_real_pexsi)) then
      call elsi_deallocate(e_h,e_h%ham_real_pexsi,"ham_real_pexsi")
   endif
   if(allocated(e_h%ham_cmplx_pexsi)) then
      call elsi_deallocate(e_h,e_h%ham_cmplx_pexsi,"ham_cmplx_pexsi")
   endif
   if(allocated(e_h%ovlp_real_pexsi)) then
      call elsi_deallocate(e_h,e_h%ovlp_real_pexsi,"ovlp_real_pexsi")
   endif
   if(allocated(e_h%ovlp_cmplx_pexsi)) then
      call elsi_deallocate(e_h,e_h%ovlp_cmplx_pexsi,"ovlp_cmplx_pexsi")
   endif
   if(allocated(e_h%dm_real_pexsi)) then
      call elsi_deallocate(e_h,e_h%dm_real_pexsi,"dm_real_pexsi")
   endif
   if(allocated(e_h%dm_cmplx_pexsi)) then
      call elsi_deallocate(e_h,e_h%dm_cmplx_pexsi,"dm_cmplx_pexsi")
   endif
   if(allocated(e_h%row_ind_pexsi)) then
      call elsi_deallocate(e_h,e_h%row_ind_pexsi,"row_ind_pexsi")
   endif
   if(allocated(e_h%col_ptr_pexsi)) then
      call elsi_deallocate(e_h,e_h%col_ptr_pexsi,"col_ptr_pexsi")
   endif
   if(allocated(e_h%ne_vec)) then
      call elsi_deallocate(e_h,e_h%ne_vec,"ne_vec")
   endif

   ! CheSS
   if(allocated(e_h%ham_real_chess)) then
      call elsi_deallocate(e_h,e_h%ham_real_chess,"ham_real_chess")
   endif
   if(allocated(e_h%ovlp_real_chess)) then
      call elsi_deallocate(e_h,e_h%ovlp_real_chess,"ovlp_real_chess")
   endif
   if(allocated(e_h%row_ind_chess)) then
      call elsi_deallocate(e_h,e_h%row_ind_chess,"row_ind_chess")
   endif
   if(allocated(e_h%col_ptr_chess)) then
      call elsi_deallocate(e_h,e_h%col_ptr_chess,"col_ptr_chess")
   endif
   if(allocated(e_h%row_ind_buf)) then
      call elsi_deallocate(e_h,e_h%row_ind_buf,"row_ind_buf")
   endif
   if(allocated(e_h%col_ptr_buf)) then
      call elsi_deallocate(e_h,e_h%col_ptr_buf,"col_ptr_buf")
   endif

   ! SIPs
   if(allocated(e_h%ham_real_sips)) then
      call elsi_deallocate(e_h,e_h%ham_real_sips,"ham_real_sips")
   endif
   if(allocated(e_h%ham_cmplx_sips)) then
      call elsi_deallocate(e_h,e_h%ham_cmplx_sips,"ham_cmplx_sips")
   endif
   if(allocated(e_h%ovlp_real_sips)) then
      call elsi_deallocate(e_h,e_h%ovlp_real_sips,"ovlp_real_sips")
   endif
   if(allocated(e_h%ovlp_cmplx_sips)) then
      call elsi_deallocate(e_h,e_h%ovlp_cmplx_sips,"ovlp_cmplx_sips")
   endif
   if(allocated(e_h%dm_real_sips)) then
      call elsi_deallocate(e_h,e_h%dm_real_sips,"dm_real_sips")
   endif
   if(allocated(e_h%dm_cmplx_sips)) then
      call elsi_deallocate(e_h,e_h%dm_cmplx_sips,"dm_cmplx_sips")
   endif
   if(allocated(e_h%row_ind_sips)) then
      call elsi_deallocate(e_h,e_h%row_ind_sips,"row_ind_sips")
   endif
   if(allocated(e_h%col_ptr_sips)) then
      call elsi_deallocate(e_h,e_h%col_ptr_sips,"col_ptr_sips")
   endif
   if(allocated(e_h%slices)) then
      call elsi_deallocate(e_h,e_h%slices,"slices")
   endif

   ! DMP
   if(allocated(e_h%ovlp_real_inv)) then
      call elsi_deallocate(e_h,e_h%ovlp_real_inv,"ovlp_real_inv")
   endif
   if(allocated(e_h%evec1)) then
      call elsi_deallocate(e_h,e_h%evec1,"evec1")
   endif
   if(allocated(e_h%evec2)) then
      call elsi_deallocate(e_h,e_h%evec2,"evec2")
   endif

   ! Auxiliary
   if(allocated(e_h%ham_real_copy)) then
      call elsi_deallocate(e_h,e_h%ham_real_copy,"ham_real_copy")
   endif
   if(allocated(e_h%ovlp_real_copy)) then
      call elsi_deallocate(e_h,e_h%ovlp_real_copy,"ovlp_real_copy")
   endif
   if(allocated(e_h%ovlp_cmplx_copy)) then
      call elsi_deallocate(e_h,e_h%ovlp_cmplx_copy,"ovlp_cmplx_copy")
   endif
   if(allocated(e_h%loc_row)) then
      call elsi_deallocate(e_h,e_h%loc_row,"loc_row")
   endif
   if(allocated(e_h%loc_col)) then
      call elsi_deallocate(e_h,e_h%loc_col,"loc_col")
   endif

   ! Finalize ELPA
   if(e_h%elpa_started) then
      call MPI_Comm_free(e_h%mpi_comm_row,ierr)
      call MPI_Comm_free(e_h%mpi_comm_col,ierr)
   endif

   ! Finalize PEXSI
   if(e_h%pexsi_started) then
      call f_ppexsi_plan_finalize(e_h%pexsi_plan,ierr)
      call MPI_Comm_free(e_h%comm_among_pole,ierr)
      call MPI_Comm_free(e_h%comm_in_pole,ierr)
      call MPI_Comm_free(e_h%comm_among_point,ierr)
      call MPI_Comm_free(e_h%comm_in_point,ierr)
   endif

   ! Finalize CheSS
   if(e_h%chess_started) then
      call deallocate_sparse_matrix(e_h%sparse_mat(1))
      call deallocate_sparse_matrix(e_h%sparse_mat(2))
      call foe_data_deallocate(e_h%ice_obj)
      call foe_data_deallocate(e_h%foe_obj)
      call deallocate_matrices(e_h%ham_chess)
      call deallocate_matrices(e_h%ovlp_chess)
      call deallocate_matrices(e_h%dm_chess)
      call deallocate_matrices(e_h%edm_chess)
      call deallocate_matrices(e_h%ovlp_inv_sqrt(1))
      call f_lib_finalize()
   endif

   ! Finalize SIPs
   if(e_h%sips_started) then
      call clean_qetsc()
   endif

   ! Print final timings
   if(e_h%handle_ready.and.e_h%output_solver_timings) then
      if(e_h%solver_timings_file%format == JSON) then
         ! Closing bracket to signify end of JSON array
         call truncate_string(e_h%solver_timings_file%prefix,2)
         call elsi_say(e_h, "]", e_h%solver_timings_file)
      end if
   end if
   call elsi_finalize_timings(e_h%solver_timings) 

   ! Close open files and finalize file IO handles other than stdio
   if(e_h%handle_ready.and.e_h%output_solver_timings.and.e_h%myid_all.eq.0) then
      close(e_h%solver_timings_file%print_unit)
   end if
   call elsi_reset_file_io_handle(e_h%solver_timings_file)

   ! Close the stdio file handle, then reset e_h
   call elsi_reset_file_io_handle(e_h%stdio)
   call elsi_reset_handle(e_h)

end subroutine

end module ELSI_SETUP
