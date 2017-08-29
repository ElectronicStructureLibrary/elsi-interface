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
   use ELSI_DATATYPE, only: elsi_handle
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
subroutine elsi_init(elsi_h,solver,parallel_mode,matrix_format,n_basis,&
              n_electron,n_state)

   implicit none

   type(elsi_handle), intent(out) :: elsi_h        !< Handle
   integer(kind=i4),  intent(in)  :: solver        !< AUTO,ELPA,LIBOMM,PEXSI,CHESS,SIPS
   integer(kind=i4),  intent(in)  :: parallel_mode !< SINGLE_PROC,MULTI_PROC
   integer(kind=i4),  intent(in)  :: matrix_format !< BLACS_DENSE,PEXSI_CSC
   integer(kind=i4),  intent(in)  :: n_basis       !< Number of basis functions
   real(kind=r8),     intent(in)  :: n_electron    !< Number of electrons
   integer(kind=i4),  intent(in)  :: n_state       !< Number of states

   character*40, parameter :: caller = "elsi_init"

   ! For safety
   call elsi_cleanup(elsi_h)

   elsi_h%handle_ready   = .true.
   elsi_h%n_basis        = n_basis
   elsi_h%n_nonsing      = n_basis
   elsi_h%n_electrons    = n_electron
   elsi_h%n_states       = n_state
   elsi_h%n_states_solve = n_state
   elsi_h%n_states_omm   = nint(n_electron/2.0_r8)
   elsi_h%solver         = solver
   elsi_h%matrix_format  = matrix_format
   elsi_h%parallel_mode  = parallel_mode
   elsi_h%n_elsi_calls   = 0

   if(parallel_mode == SINGLE_PROC) then
      elsi_h%n_l_rows    = n_basis
      elsi_h%n_l_cols    = n_basis
      elsi_h%n_b_rows    = n_basis
      elsi_h%n_b_cols    = n_basis
      elsi_h%myid        = 0
      elsi_h%n_procs     = 1
      elsi_h%myid_all    = 0
      elsi_h%n_procs_all = 1
   endif

   ! Set ELPA default
   if(solver == ELPA) then
      call elsi_set_elpa_default(elsi_h)
   endif

   ! Set libOMM default
   if(solver == LIBOMM) then
      call elsi_set_elpa_default(elsi_h)
      call elsi_set_omm_default(elsi_h)
   endif

   ! Set PEXSI default
   if(solver == PEXSI) then
      call elsi_set_pexsi_default(elsi_h)
   endif

   ! Set CheSS default
   if(solver == CHESS) then
      call elsi_set_chess_default(elsi_h)
   endif

   ! Set SIPs default
   if(solver == SIPS) then
      call elsi_set_sips_default(elsi_h)
   endif

   call elsi_init_timers(elsi_h)

end subroutine

!>
!! This routine sets the MPI communicator.
!!
subroutine elsi_set_mpi(elsi_h,mpi_comm)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h   !< Handle
   integer(kind=i4),  intent(in)    :: mpi_comm !< Unit ELSI communicator

   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_set_mpi"

   call elsi_check_handle(elsi_h,caller)

   if(elsi_h%parallel_mode == MULTI_PROC) then
      ! Unit ELSI communicator
      ! If there's more than one spin/kpt, each unit communicator
      ! solves one KS problem
      elsi_h%mpi_comm = mpi_comm

      call MPI_Comm_rank(mpi_comm,elsi_h%myid,mpierr)
      call MPI_Comm_size(mpi_comm,elsi_h%n_procs,mpierr)

      elsi_h%mpi_ready = .true.
   endif

end subroutine

!>
!! This routine sets the global MPI communicator.
!!
subroutine elsi_set_mpi_global(elsi_h,mpi_comm_all)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h       !< Handle
   integer(kind=i4),  intent(in)    :: mpi_comm_all !< Unit ELSI communicator

   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_set_mpi_global"

   call elsi_check_handle(elsi_h,caller)

   if(elsi_h%parallel_mode == MULTI_PROC) then
      ! Global ELSI communicator
      elsi_h%mpi_comm_all = mpi_comm_all

      call MPI_Comm_rank(mpi_comm_all,elsi_h%myid_all,mpierr)
      call MPI_Comm_size(mpi_comm_all,elsi_h%n_procs_all,mpierr)

      elsi_h%global_mpi_ready = .true.
   endif

end subroutine

!>
!! This routine sets the spin information.
!!
subroutine elsi_set_spin(elsi_h,n_spin,i_spin)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   integer(kind=i4),  intent(in)    :: n_spin !< Number of spin channels
   integer(kind=i4),  intent(in)    :: i_spin !< Spin index

   elsi_h%n_spins = n_spin
   elsi_h%i_spin  = i_spin

end subroutine

!>
!! This routine sets the k-point information.
!!
subroutine elsi_set_kpoint(elsi_h,n_kpt,i_kpt,weight)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle
   integer(kind=i4),  intent(in)    :: n_kpt  !< Number of k-points
   integer(kind=i4),  intent(in)    :: i_kpt  !< K-point index
   real(kind=r8),     intent(in)    :: weight !< Weight

   elsi_h%n_kpts   = n_kpt
   elsi_h%i_kpt    = i_kpt
   elsi_h%i_weight = weight

end subroutine

!>
!! This routine sets the BLACS context and the block size.
!!
subroutine elsi_set_blacs(elsi_h,blacs_ctxt,block_size)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h     !< Handle
   integer(kind=i4),  intent(in)    :: blacs_ctxt !< BLACS context
   integer(kind=i4),  intent(in)    :: block_size !< Block size

   integer(kind=i4) :: i,i_row,i_col
   integer(kind=i4) :: blacs_info

   integer(kind=i4), external :: numroc

   character*40, parameter :: caller = "elsi_set_blacs"

   call elsi_check_handle(elsi_h,caller)

   if(elsi_h%parallel_mode == MULTI_PROC) then
      elsi_h%blacs_ctxt = blacs_ctxt
      elsi_h%n_b_rows = block_size
      elsi_h%n_b_cols = block_size

      ! Get processor grid information
      call blacs_gridinfo(elsi_h%blacs_ctxt,elsi_h%n_p_rows,elsi_h%n_p_cols,&
              elsi_h%my_p_row,elsi_h%my_p_col)

      ! Get local size of matrix
      elsi_h%n_l_rows = numroc(elsi_h%n_basis,elsi_h%n_b_rows,&
         elsi_h%my_p_row,0,elsi_h%n_p_rows)
      elsi_h%n_l_cols = numroc(elsi_h%n_basis,elsi_h%n_b_cols,&
         elsi_h%my_p_col,0,elsi_h%n_p_cols)

      ! Get BLACS descriptor
      call descinit(elsi_h%sc_desc,elsi_h%n_basis,elsi_h%n_basis,elsi_h%n_b_rows,&
              elsi_h%n_b_cols,0,0,elsi_h%blacs_ctxt,max(1,elsi_h%n_l_rows),blacs_info)

      ! Get ELPA communicators
      call elsi_get_elpa_comms(elsi_h)

      ! Compute global-local mapping
      call elsi_allocate(elsi_h,elsi_h%local_row,elsi_h%n_basis,"local_row",caller)
      call elsi_allocate(elsi_h,elsi_h%local_col,elsi_h%n_basis,"local_col",caller)

      i_row = 0
      i_col = 0

      do i = 1,elsi_h%n_basis
         if(mod((i-1)/elsi_h%n_b_rows,elsi_h%n_p_rows) == elsi_h%my_p_row) then
            i_row = i_row+1
            elsi_h%local_row(i) = i_row
         endif
         if(mod((i-1)/elsi_h%n_b_cols,elsi_h%n_p_cols) == elsi_h%my_p_col) then
            i_col = i_col+1
            elsi_h%local_col(i) = i_col
         endif
      enddo

      ! Set up MatrixSwitch
      if(elsi_h%solver == LIBOMM) then
         call ms_scalapack_setup(elsi_h%mpi_comm,elsi_h%n_p_rows,'r',&
                 elsi_h%n_b_rows,icontxt=elsi_h%blacs_ctxt)
      endif

      elsi_h%blacs_ready = .true.
   endif

end subroutine

!>
!! This routine sets the sparsity pattern.
!!
subroutine elsi_set_csc(elsi_h,nnz_g,nnz_l,n_l_cols,row_ind,col_ptr)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h              !< Handle
   integer(kind=i4),  intent(in)    :: nnz_g               !< Global number of nonzeros
   integer(kind=i4),  intent(in)    :: nnz_l               !< Local number of nonzeros
   integer(kind=i4),  intent(in)    :: n_l_cols            !< Local number of columns
   integer(kind=i4)                 :: row_ind(nnz_l)      !< Row index
   integer(kind=i4)                 :: col_ptr(n_l_cols+1) !< Column pointer

   character*40, parameter :: caller = "elsi_set_csc"

   call elsi_check_handle(elsi_h,caller)

   elsi_h%nnz_g       = nnz_g
   elsi_h%nnz_l_sp    = nnz_l
   elsi_h%n_l_cols_sp = n_l_cols

   call elsi_set_row_ind(elsi_h,row_ind)
   call elsi_set_col_ptr(elsi_h,col_ptr)

   elsi_h%sparsity_ready = .true.

end subroutine

!>
!! This routine finalizes ELSI.
!!
subroutine elsi_finalize(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   character*40, parameter :: caller = "elsi_finalize"

   call elsi_check_handle(elsi_h,caller)
   call elsi_final_print(elsi_h)
   call elsi_cleanup(elsi_h)

end subroutine

end module ELSI_SETUP
