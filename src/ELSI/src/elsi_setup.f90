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
   use ELSI_CONSTANTS,     only: ELPAA,LIBOMM,PEXSI,CHESS,SIPS,SINGLE_PROC,&
                                 MULTI_PROC
   use ELSI_DATATYPE
   use ELSI_ELPA,          only: elsi_set_elpa_default,elsi_get_elpa_comms
   use ELSI_MALLOC
   use ELSI_MATRICES,      only: elsi_set_row_ind,elsi_set_col_ptr
   use ELSI_OMM,           only: elsi_set_omm_default
   use ELSI_PEXSI,         only: elsi_set_pexsi_default
   use ELSI_PRECISION,     only: r8,i4
   use ELSI_SIPS,          only: elsi_set_sips_default
   use ELSI_UTILS
   use FOE_BASE,           only: foe_data_deallocate
   use F_PPEXSI_INTERFACE, only: f_ppexsi_plan_finalize
   use MATRIXSWITCH,       only: ms_scalapack_setup,m_deallocate
   use M_QETSC,            only: clean_qetsc
   use SPARSEMATRIX_BASE,  only: deallocate_sparse_matrix,deallocate_matrices

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
subroutine elsi_init(e_h,solver,parallel_mode,matrix_format,n_basis,n_electron,&
              n_state)

   implicit none

   type(elsi_handle), intent(out) :: e_h           !< Handle
   integer(kind=i4),  intent(in)  :: solver        !< AUTO,ELPA,LIBOMM,PEXSI,CHESS,SIPS
   integer(kind=i4),  intent(in)  :: parallel_mode !< SINGLE_PROC,MULTI_PROC
   integer(kind=i4),  intent(in)  :: matrix_format !< BLACS_DENSE,PEXSI_CSC
   integer(kind=i4),  intent(in)  :: n_basis       !< Number of basis functions
   real(kind=r8),     intent(in)  :: n_electron    !< Number of electrons
   integer(kind=i4),  intent(in)  :: n_state       !< Number of states

   character*40, parameter :: caller = "elsi_init"

   ! For safety
   call elsi_cleanup(e_h)

   e_h%handle_ready   = .true.
   e_h%n_basis        = n_basis
   e_h%n_nonsing      = n_basis
   e_h%n_electrons    = n_electron
   e_h%n_states       = n_state
   e_h%n_states_solve = n_state
   e_h%n_states_omm   = nint(n_electron/2.0_r8)
   e_h%solver         = solver
   e_h%matrix_format  = matrix_format
   e_h%parallel_mode  = parallel_mode
   e_h%n_elsi_calls   = 0

   if(parallel_mode == SINGLE_PROC) then
      e_h%n_l_rows    = n_basis
      e_h%n_l_cols    = n_basis
      e_h%n_b_rows    = n_basis
      e_h%n_b_cols    = n_basis
      e_h%myid        = 0
      e_h%n_procs     = 1
      e_h%myid_all    = 0
      e_h%n_procs_all = 1
   endif

   ! Set ELPA default
   if(solver == ELPAA) then
      call elsi_set_elpa_default(e_h)
   endif

   ! Set libOMM default
   if(solver == LIBOMM) then
      call elsi_set_elpa_default(e_h)
      call elsi_set_omm_default(e_h)
   endif

   ! Set PEXSI default
   if(solver == PEXSI) then
      call elsi_set_pexsi_default(e_h)
   endif

   ! Set CheSS default
   if(solver == CHESS) then
      call elsi_set_chess_default(e_h)
   endif

   ! Set SIPs default
   if(solver == SIPS) then
      call elsi_set_sips_default(e_h)
   endif

   call elsi_init_timer(e_h)

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

   integer(kind=i4) :: i,i_row,i_col
   integer(kind=i4) :: blacs_info

   integer(kind=i4), external :: numroc

   character*40, parameter :: caller = "elsi_set_blacs"

   call elsi_check_handle(e_h,caller)

   if(e_h%parallel_mode == MULTI_PROC) then
      e_h%blacs_ctxt = blacs_ctxt
      e_h%n_b_rows = block_size
      e_h%n_b_cols = block_size

      ! Get processor grid information
      call blacs_gridinfo(e_h%blacs_ctxt,e_h%n_p_rows,e_h%n_p_cols,&
              e_h%my_p_row,e_h%my_p_col)

      ! Get local size of matrix
      e_h%n_l_rows = numroc(e_h%n_basis,e_h%n_b_rows,e_h%my_p_row,0,&
                        e_h%n_p_rows)
      e_h%n_l_cols = numroc(e_h%n_basis,e_h%n_b_cols,e_h%my_p_col,0,&
                        e_h%n_p_cols)

      ! Get BLACS descriptor
      call descinit(e_h%sc_desc,e_h%n_basis,e_h%n_basis,e_h%n_b_rows,&
              e_h%n_b_cols,0,0,e_h%blacs_ctxt,max(1,e_h%n_l_rows),blacs_info)

      ! Get ELPA communicators
      call elsi_get_elpa_comms(e_h)

      ! Create global-local mapping
      call elsi_allocate(e_h,e_h%loc_row,e_h%n_basis,"loc_row",caller)
      call elsi_allocate(e_h,e_h%loc_col,e_h%n_basis,"loc_col",caller)

      i_row = 0
      i_col = 0

      do i = 1,e_h%n_basis
         if(mod((i-1)/e_h%n_b_rows,e_h%n_p_rows) == e_h%my_p_row) then
            i_row = i_row+1
            e_h%loc_row(i) = i_row
         endif
         if(mod((i-1)/e_h%n_b_cols,e_h%n_p_cols) == e_h%my_p_col) then
            i_col = i_col+1
            e_h%loc_col(i) = i_col
         endif
      enddo

      ! Set up MatrixSwitch
      if(e_h%solver == LIBOMM) then
         call ms_scalapack_setup(e_h%mpi_comm,e_h%n_p_rows,'r',e_h%n_b_rows,&
                 icontxt=e_h%blacs_ctxt)
      endif

      e_h%blacs_ready = .true.
   endif

end subroutine

!>
!! This routine sets the sparsity pattern.
!!
subroutine elsi_set_csc(e_h,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                 !< Handle
   integer(kind=i4),  intent(in)    :: nnz_g               !< Global number of nonzeros
   integer(kind=i4),  intent(in)    :: nnz_l               !< Local number of nonzeros
   integer(kind=i4),  intent(in)    :: n_l_cols            !< Local number of columns
   integer(kind=i4),  intent(in)    :: row_ind(nnz_l)      !< Row index
   integer(kind=i4),  intent(in)    :: col_ptr(n_l_cols+1) !< Column pointer

   character*40, parameter :: caller = "elsi_set_csc"

   call elsi_check_handle(e_h,caller)

   e_h%nnz_g       = nnz_g
   e_h%nnz_l_sp    = nnz_l
   e_h%n_l_cols_sp = n_l_cols

   call elsi_set_row_ind(e_h,row_ind)
   call elsi_set_col_ptr(e_h,col_ptr)

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

   type(elsi_handle), intent(in) :: e_h !< Handle

   real(kind=r8) :: sparsity
   character*200 :: info_str

   character*40, parameter :: caller = "elsi_final_print"

   call elsi_say("  |------------------------------------------",e_h)
   call elsi_say("  | Final ELSI Output                        ",e_h)
   call elsi_say("  |------------------------------------------",e_h)

   write(info_str,"(A,I13)") "  | Number of basis functions :",e_h%n_basis
   call elsi_say(info_str,e_h)

   if(e_h%parallel_mode == MULTI_PROC) then
      write(info_str,"(A,I13)") "  | Number of nonzeros        :",e_h%nnz_g
      call elsi_say(info_str,e_h)

      sparsity = 1.0_r8-(1.0_r8*e_h%nnz_g/e_h%n_basis/e_h%n_basis)
      write(info_str,"(A,F13.3)") "  | Sparsity                  :",sparsity
      call elsi_say(info_str,e_h)
   endif

   write(info_str,"(A,F13.1)") "  | Number of electrons       :",e_h%n_electrons
   call elsi_say(info_str,e_h)
   write(info_str,"(A,I13)") "  | Number of spins           :",e_h%n_spins
   call elsi_say(info_str,e_h)
   write(info_str,"(A,I13)") "  | Number of k-points        :",e_h%n_kpts
   call elsi_say(info_str,e_h)

   if(e_h%solver == ELPAA .or. e_h%solver == SIPS) then
      write(info_str,"(A,I13)") "  | Number of states          :",e_h%n_states
      call elsi_say(info_str,e_h)
   endif

   if(e_h%solver == ELPAA) then
      call elsi_say("  | Solver                    :         ELPA ",e_h)
   elseif(e_h%solver == LIBOMM) then
      call elsi_say("  | Solver                    :       libOMM ",e_h)
   elseif(e_h%solver == PEXSI) then
      call elsi_say("  | Solver                    :        PEXSI ",e_h)
   elseif(e_h%solver == CHESS) then
      call elsi_say("  | Solver                    :        CheSS ",e_h)
   elseif(e_h%solver == SIPS) then
      call elsi_say("  | Solver                    :         SIPs ",e_h)
   endif

   if(e_h%parallel_mode == MULTI_PROC) then
      call elsi_say("  | Parallel mode             :   MULTI_PROC ",e_h)
   elseif(e_h%parallel_mode == SINGLE_PROC) then
      call elsi_say("  | Parallel mode             :  SINGLE_PROC ",e_h)
   endif

   if(e_h%matrix_format == BLACS_DENSE) then
      call elsi_say("  | Matrix format             :  BLACS_DENSE ",e_h)
   elseif(e_h%matrix_format == PEXSI_CSC) then
      call elsi_say("  | Matrix format             :    PEXSI_CSC ",e_h)
   endif

   call elsi_say("  |------------------------------------------",e_h)
   call elsi_say("  | ELSI Project (c)  elsi-interchange.org   ",e_h)
   call elsi_say("  |------------------------------------------",e_h)

end subroutine

!>
!! This routine frees memory.
!!
subroutine elsi_cleanup(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   integer(kind=i4) :: ierr

   character*40, parameter :: caller = "elsi_cleanup"

   ! Nullify pointers
   if(associated(e_h%ham_real)) then
      nullify(e_h%ham_real)
   endif
   if(associated(e_h%ham_cmplx)) then
      nullify(e_h%ham_cmplx)
   endif
   if(associated(e_h%ovlp_real)) then
      nullify(e_h%ovlp_real)
   endif
   if(associated(e_h%ovlp_cmplx)) then
      nullify(e_h%ovlp_cmplx)
   endif
   if(associated(e_h%eval)) then
      nullify(e_h%eval)
   endif
   if(associated(e_h%evec_real)) then
      nullify(e_h%evec_real)
   endif
   if(associated(e_h%evec_cmplx)) then
      nullify(e_h%evec_cmplx)
   endif
   if(associated(e_h%dm_real)) then
      nullify(e_h%dm_real)
   endif
   if(associated(e_h%dm_cmplx)) then
      nullify(e_h%dm_cmplx)
   endif
   if(associated(e_h%ham_real_ccs)) then
      nullify(e_h%ham_real_ccs)
   endif
   if(associated(e_h%ham_cmplx_ccs)) then
      nullify(e_h%ham_cmplx_ccs)
   endif
   if(associated(e_h%ovlp_real_ccs)) then
      nullify(e_h%ovlp_real_ccs)
   endif
   if(associated(e_h%ovlp_cmplx_ccs)) then
      nullify(e_h%ovlp_cmplx_ccs)
   endif
   if(associated(e_h%dm_real_ccs)) then
      nullify(e_h%dm_real_ccs)
   endif
   if(associated(e_h%dm_cmplx_ccs)) then
      nullify(e_h%dm_cmplx_ccs)
   endif
   if(associated(e_h%row_ind_ccs)) then
      nullify(e_h%row_ind_ccs)
   endif
   if(associated(e_h%col_ptr_ccs)) then
      nullify(e_h%col_ptr_ccs)
   endif

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
   if(allocated(e_h%ovlp_real_omm)) then
      call elsi_deallocate(e_h,e_h%ovlp_real_omm,"ovlp_real_omm")
   endif
   if(allocated(e_h%ovlp_cmplx_omm)) then
      call elsi_deallocate(e_h,e_h%ovlp_cmplx_omm,"ovlp_cmplx_omm")
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

   if(allocated(e_h%loc_row)) then
      call elsi_deallocate(e_h,e_h%loc_row,"loc_row")
   endif
   if(allocated(e_h%loc_col)) then
      call elsi_deallocate(e_h,e_h%loc_col,"loc_col")
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

   ! Reset e_h
   call elsi_reset_handle(e_h)

end subroutine

end module ELSI_SETUP
