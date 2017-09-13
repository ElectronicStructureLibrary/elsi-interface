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

   use ELSI_CHESS, only: elsi_set_chess_default
   use ELSI_CONSTANTS, only: ELPA,LIBOMM,PEXSI,CHESS,SIPS,SINGLE_PROC,MULTI_PROC
   use ELSI_DATATYPE
   use ELSI_ELPA, only: elsi_set_elpa_default,elsi_get_elpa_comms
   use ELSI_OMM, only: elsi_set_omm_default
   use ELSI_PEXSI, only: elsi_set_pexsi_default
   use ELSI_PRECISION, only: r8,i4
   use ELSI_SIPS, only: elsi_set_sips_default
   use ELSI_UTILS
   use MATRIXSWITCH, only: ms_scalapack_setup

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
   if(solver == ELPA) then
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

   call elsi_init_timers(e_h)

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

   write(info_str,"(A)") "  |------------------------------------------"
   call elsi_statement_print(info_str,e_h)
   write(info_str,"(A)") "  | Final ELSI Output"
   call elsi_statement_print(info_str,e_h)
   write(info_str,"(A)") "  |------------------------------------------"
   call elsi_statement_print(info_str,e_h)

   write(info_str,"(A,I13)") "  | Number of basis functions :",e_h%n_basis
   call elsi_statement_print(info_str,e_h)

   if(e_h%solver == PEXSI .or. e_h%solver == SIPS) then
      write(info_str,"(A,I13)") "  | Number of nonzeros        :",e_h%nnz_g
      call elsi_statement_print(info_str,e_h)

      sparsity = 1.0_r8-(1.0_r8*e_h%nnz_g/e_h%n_basis/e_h%n_basis)
      write(info_str,"(A,F13.3)") "  | Sparsity                  :",sparsity
      call elsi_statement_print(info_str,e_h)
   endif

   write(info_str,"(A,F13.1)") "  | Number of electrons       :",e_h%n_electrons
   call elsi_statement_print(info_str,e_h)
   write(info_str,"(A,I13)") "  | Number of spins           :",e_h%n_spins
   call elsi_statement_print(info_str,e_h)
   write(info_str,"(A,I13)") "  | Number of k-points        :",e_h%n_kpts
   call elsi_statement_print(info_str,e_h)

   if(e_h%solver == ELPA .or. e_h%solver == SIPS) then
      write(info_str,"(A,I13)") "  | Number of states          :",e_h%n_states
      call elsi_statement_print(info_str,e_h)
   endif

   if(e_h%solver == ELPA) then
      write(info_str,"(A,A13)") "  | Solver                    :","ELPA"
      call elsi_statement_print(info_str,e_h)
   elseif(e_h%solver == LIBOMM) then
      write(info_str,"(A,A13)") "  | Solver                    :","libOMM"
      call elsi_statement_print(info_str,e_h)
   elseif(e_h%solver == PEXSI) then
      write(info_str,"(A,A13)") "  | Solver                    :","PEXSI"
      call elsi_statement_print(info_str,e_h)
   elseif(e_h%solver == CHESS) then
      write(info_str,"(A,A13)") "  | Solver                    :","CheSS"
      call elsi_statement_print(info_str,e_h)
   elseif(e_h%solver == SIPS) then
      write(info_str,"(A,A13)") "  | Solver                    :","SIPs"
      call elsi_statement_print(info_str,e_h)
   endif

   if(e_h%parallel_mode == MULTI_PROC) then
      write(info_str,"(A,A13)") "  | Parallel mode             :","MULTI_PROC"
      call elsi_statement_print(info_str,e_h)
   elseif(e_h%parallel_mode == SINGLE_PROC) then
      write(info_str,"(A,A13)") "  | Parallel mode             :","SINGLE_PROC"
      call elsi_statement_print(info_str,e_h)
   endif

   if(e_h%matrix_format == BLACS_DENSE) then
      write(info_str,"(A,A13)") "  | Matrix format             :","BLACS_DENSE"
      call elsi_statement_print(info_str,e_h)
   elseif(e_h%matrix_format == PEXSI_CSC) then
      write(info_str,"(A,A13)") "  | Matrix format             :","PEXSI_CSC"
      call elsi_statement_print(info_str,e_h)
   endif

   write(info_str,"(A)") "  |------------------------------------------"
   call elsi_statement_print(info_str,e_h)
   write(info_str,"(A)") "  | ELSI Project (c)  elsi-interchange.org"
   call elsi_statement_print(info_str,e_h)
   write(info_str,"(A)") "  |------------------------------------------"
   call elsi_statement_print(info_str,e_h)

end subroutine

end module ELSI_SETUP
