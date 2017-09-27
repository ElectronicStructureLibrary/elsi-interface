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
!! This module performs parallel matrix IO.
!!
module ELSI_IO

   use, intrinsic :: ISO_C_BINDING
   use ELSI_CONSTANTS, only: HEADER_SIZE,BLACS_DENSE,DENSE_FILE,CSC_FILE,&
                             READ_FILE,WRITE_FILE,REAL_VALUES,COMPLEX_VALUES,&
                             FILE_VERSION
   use ELSI_DATATYPE
   use ELSI_MALLOC
   use ELSI_MATCONV,   only: elsi_pexsi_to_blacs_dm,elsi_blacs_to_sips_hs
   use ELSI_MATRICES,  only: elsi_set_sparse_dm
   use ELSI_PRECISION, only: r8,i4,i8
   use ELSI_SETUP,     only: elsi_init,elsi_set_mpi,elsi_set_blacs,&
                             elsi_set_csc,elsi_cleanup
   use ELSI_UTILS

   implicit none

   private

   public :: elsi_init_rw
   public :: elsi_finalize_rw
   public :: elsi_set_rw_mpi
   public :: elsi_set_rw_blacs
   public :: elsi_set_rw_csc
   public :: elsi_set_rw_output
   public :: elsi_set_rw_write_unit
   public :: elsi_set_rw_zero_def
   public :: elsi_set_rw_header
   public :: elsi_get_rw_header
   public :: elsi_read_mat_dim
   public :: elsi_read_mat_dim_sparse
   public :: elsi_read_mat_real
   public :: elsi_read_mat_real_sparse
   public :: elsi_read_mat_complex
   public :: elsi_read_mat_complex_sparse
   public :: elsi_write_mat_real
   public :: elsi_write_mat_real_sparse
   public :: elsi_write_mat_complex
   public :: elsi_write_mat_complex_sparse

contains

!>
!! This routine initializes a handle for reading and writing matrices.
!!
subroutine elsi_init_rw(rw_h,task,parallel_mode,file_format,n_basis,n_electron)

   implicit none

   type(elsi_rw_handle), intent(out) :: rw_h          !< Handle
   integer(kind=i4),     intent(in)  :: task          !< READ_MAT,WRITE_MAT
   integer(kind=i4),     intent(in)  :: parallel_mode !< SINGLE_PROC,MULTI_PROC
   integer(kind=i4),     intent(in)  :: file_format   !< DENSE_FILE,CSC_FILE
   integer(kind=i4),     intent(in)  :: n_basis       !< Number of basis functions
   real(kind=r8),        intent(in)  :: n_electron    !< Number of electrons

   character*40, parameter :: caller = "elsi_init_io"

   ! For safety
   call elsi_reset_rw_handle(rw_h)

   rw_h%handle_ready  = .true.
   rw_h%rw_task       = task
   rw_h%parallel_mode = parallel_mode
   rw_h%file_format   = file_format
   rw_h%n_basis       = n_basis
   rw_h%n_electrons   = n_electron

   if(parallel_mode == SINGLE_PROC) then
      rw_h%n_lrow  = n_basis
      rw_h%n_lcol  = n_basis
      rw_h%myid    = 0
      rw_h%n_procs = 1
   endif

end subroutine

!>
!! This routine finalizes a handle.
!!
subroutine elsi_finalize_rw(rw_h)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h !< Handle

   character*40, parameter :: caller = "elsi_finalize_rw"

   call elsi_check_rw_handle(rw_h,caller)
   call elsi_reset_rw_handle(rw_h)

end subroutine

!>
!! This routine sets the MPI communicator.
!!
subroutine elsi_set_rw_mpi(rw_h,mpi_comm)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h     !< Handle
   integer(kind=i4),     intent(in)    :: mpi_comm !< Communicator

   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_set_rw_mpi"

   call elsi_check_rw_handle(rw_h,caller)

   if(rw_h%parallel_mode == MULTI_PROC) then
      rw_h%mpi_comm = mpi_comm

      call MPI_Comm_rank(mpi_comm,rw_h%myid,mpierr)
      call MPI_Comm_size(mpi_comm,rw_h%n_procs,mpierr)

      rw_h%mpi_ready = .true.
   endif

end subroutine

!>
!! This routine sets the BLACS context and the block size.
!!
subroutine elsi_set_rw_blacs(rw_h,blacs_ctxt,block_size)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h       !< Handle
   integer(kind=i4),     intent(in)    :: blacs_ctxt !< BLACS context
   integer(kind=i4),     intent(in)    :: block_size !< Block size

   integer(kind=i4) :: n_prow
   integer(kind=i4) :: n_pcol
   integer(kind=i4) :: my_prow
   integer(kind=i4) :: my_pcol

   integer(kind=i4), external :: numroc

   character*40, parameter :: caller = "elsi_set_rw_blacs"

   call elsi_check_rw_handle(rw_h,caller)

   if(rw_h%parallel_mode == MULTI_PROC) then
      rw_h%blacs_ctxt = blacs_ctxt
      rw_h%blk        = block_size

      if(rw_h%rw_task == WRITE_FILE) then
         ! Get processor grid information
         call blacs_gridinfo(rw_h%blacs_ctxt,n_prow,n_pcol,my_prow,my_pcol)

         ! Get local size of matrix
         rw_h%n_lrow = numroc(rw_h%n_basis,rw_h%blk,my_prow,0,n_prow)
         rw_h%n_lcol = numroc(rw_h%n_basis,rw_h%blk,my_pcol,0,n_pcol)
      endif

      rw_h%blacs_ready = .true.
   endif

end subroutine

!>
!! This routine sets the sparsity pattern.
!!
subroutine elsi_set_rw_csc(rw_h,nnz_g,nnz_l_sp,n_lcol_sp)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h      !< Handle
   integer(kind=i4),     intent(in)    :: nnz_g     !< Global number of nonzeros
   integer(kind=i4),     intent(in)    :: nnz_l_sp  !< Local number of nonzeros
   integer(kind=i4),     intent(in)    :: n_lcol_sp !< Local number of columns

   character*40, parameter :: caller = "elsi_set_rw_csc"

   call elsi_check_rw_handle(rw_h,caller)

   rw_h%nnz_g     = nnz_g
   rw_h%nnz_l_sp  = nnz_l_sp
   rw_h%n_lcol_sp = n_lcol_sp

   rw_h%sparsity_ready = .true.

end subroutine

!>
!! This routine sets ELSI output level.
!!
subroutine elsi_set_rw_output(rw_h,out_level)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h      !< Handle
   integer(kind=i4),     intent(in)    :: out_level !< Output level

   character*40, parameter :: caller = "elsi_set_rw_output"

   call elsi_check_rw_handle(rw_h,caller)

   if(out_level <= 0) then
      rw_h%print_info = .false.
      rw_h%print_mem  = .false.
   elseif(out_level == 1) then
      rw_h%print_info = .true.
      rw_h%print_mem  = .false.
   elseif(out_level == 2) then
      rw_h%print_info = .true.
      rw_h%print_mem  = .false.
   else
      rw_h%print_info = .true.
      rw_h%print_mem  = .true.
   endif

end subroutine

!>
!! This routine sets the unit to be used by ELSI output.
!!
subroutine elsi_set_rw_write_unit(rw_h,write_unit)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h       !< Handle
   integer(kind=i4),     intent(in)    :: write_unit !< Unit

   character*40, parameter :: caller = "elsi_set_rw_write_unit"

   call elsi_check_rw_handle(rw_h,caller)

   rw_h%print_unit = write_unit

end subroutine

!>
!! This routine sets the threshold to define "zero".
!!
subroutine elsi_set_rw_zero_def(rw_h,zero_def)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h     !< Handle
   real(kind=r8),        intent(in)    :: zero_def !< Zero tolerance

   character*40, parameter :: caller = "elsi_set_rw_zero_def"

   call elsi_check_rw_handle(rw_h,caller)

   rw_h%zero_def = zero_def

end subroutine

!>
!! This routine sets a matrix file header.
!!
subroutine elsi_set_rw_header(rw_h,header_user)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h           !< Handle
   integer(kind=i4),     intent(in)    :: header_user(8) !< User's header

   character*40, parameter :: caller = "elsi_set_rw_header"

   call elsi_check_rw_handle(rw_h,caller)

   rw_h%header_user = header_user

end subroutine

!>
!! This routine gets a matrix file header.
!!
subroutine elsi_get_rw_header(rw_h,header_user)

   implicit none

   type(elsi_rw_handle), intent(in)  :: rw_h           !< Handle
   integer(kind=i4),     intent(out) :: header_user(8) !< User's header

   character*40, parameter :: caller = "elsi_get_rw_header"

   call elsi_check_rw_handle(rw_h,caller)

   header_user = rw_h%header_user

end subroutine

!>
!! This routine checks whether a handle has been properly initialized for
!! reading and writing matrices.
!!
subroutine elsi_check_rw_handle(rw_h,caller)

   implicit none

   type(elsi_rw_handle), intent(in) :: rw_h   !< Handle
   character(len=*),     intent(in) :: caller !< Caller

   if(.not. rw_h%handle_ready) then
      call elsi_rw_stop(" Invalid handle! Not initialized.",rw_h,caller)
   endif

end subroutine

!>
!! This routine resets a handle.
!!
subroutine elsi_reset_rw_handle(rw_h)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h !< Handle

   character*40, parameter :: caller = "elsi_reset_rw_handle"

   rw_h%handle_ready   = .false.
   rw_h%rw_task        = UNSET
   rw_h%parallel_mode  = UNSET
   rw_h%file_format    = UNSET
   rw_h%print_info     = .false.
   rw_h%print_mem      = .false.
   rw_h%print_unit     = 6
   rw_h%myid           = UNSET
   rw_h%n_procs        = UNSET
   rw_h%mpi_comm       = UNSET
   rw_h%mpi_ready      = .false.
   rw_h%blacs_ctxt     = UNSET
   rw_h%blk            = UNSET
   rw_h%n_lrow         = UNSET
   rw_h%n_lcol         = UNSET
   rw_h%blacs_ready    = .false.
   rw_h%nnz_g          = UNSET
   rw_h%nnz_l_sp       = UNSET
   rw_h%n_lcol_sp      = UNSET
   rw_h%zero_def       = 1.0e-15_r8
   rw_h%sparsity_ready = .false.
   rw_h%n_electrons    = 0.0_r8
   rw_h%n_basis        = UNSET
   rw_h%header_user    = UNSET

end subroutine

!>
!! Clean shutdown in case of errors.
!!
subroutine elsi_rw_stop(info,rw_h,caller)

   implicit none

   character(len=*),     intent(in) :: info   !< Error message
   type(elsi_rw_handle), intent(in) :: rw_h   !< Handle
   character(len=*),     intent(in) :: caller !< Caller

   character*800 :: info_str
   integer       :: mpierr

   if(rw_h%mpi_ready) then
      write(info_str,"(A,I7,5A)") "**Error! MPI task ",rw_h%myid," in ",&
         trim(caller),": ",trim(info)," Exiting..."
      write(rw_h%print_unit,"(A)") trim(info_str)

      if(rw_h%n_procs > 1) then
         call MPI_Abort(rw_h%mpi_comm,0,mpierr)
      endif
   else
      write(info_str,"(5A)") "**Error! ",trim(caller),": ",trim(info),&
         " Exiting..."
      write(rw_h%print_unit,"(A)") trim(info_str)
   endif

   stop

end subroutine

!>
!! This routine reads the dimensions of a 2D block-cyclic dense matrix from
!! file.
!!
subroutine elsi_read_mat_dim(rw_h,f_name,n_electron,n_basis,n_lrow,n_lcol)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h       !< Handle
   character(*),         intent(in)    :: f_name     !< File name
   real(kind=r8),        intent(out)   :: n_electron !< Number of electrons
   integer(kind=i4),     intent(out)   :: n_basis    !< Matrix size
   integer(kind=i4),     intent(out)   :: n_lrow     !< Local number of rows
   integer(kind=i4),     intent(out)   :: n_lcol     !< Local number of columns

   character*40, parameter :: caller = "elsi_read_mat_dim"

   if(rw_h%parallel_mode == MULTI_PROC) then
      call elsi_read_mat_dim_mp(rw_h,f_name,n_electron,n_basis,n_lrow,n_lcol)
   else
      call elsi_read_mat_dim_sp(rw_h,f_name,n_electron,n_basis,n_lrow,n_lcol)
   endif

end subroutine

!>
!! This routine reads a 2D block-cyclic dense matrix from file.
!!
subroutine elsi_read_mat_real(rw_h,f_name,mat)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h                         !< Handle
   character(*),         intent(in)    :: f_name                       !< File name
   real(kind=r8),        intent(out)   :: mat(rw_h%n_lrow,rw_h%n_lcol) !< Matrix

   character*40, parameter :: caller = "elsi_read_mat_real"

   if(rw_h%parallel_mode == MULTI_PROC) then
      call elsi_read_mat_real_mp(rw_h,f_name,mat)
   else
      call elsi_read_mat_real_sp(rw_h,f_name,mat)
   endif

end subroutine

!>
!! This routine writes a 2D block-cyclic dense matrix to file.
!!
subroutine elsi_write_mat_real(rw_h,f_name,mat)

   implicit none

   type(elsi_rw_handle), intent(in) :: rw_h                         !< Handle
   character(*),         intent(in) :: f_name                       !< File name
   real(kind=r8),        intent(in) :: mat(rw_h%n_lrow,rw_h%n_lcol) !< Matrix

   character*40, parameter :: caller = "elsi_write_mat_real"

   if(rw_h%parallel_mode == MULTI_PROC) then
      call elsi_write_mat_real_mp(rw_h,f_name,mat)
   else
      call elsi_write_mat_real_sp(rw_h,f_name,mat)
   endif

end subroutine

!>
!! This routine reads a 2D block-cyclic dense matrix from file.
!!
subroutine elsi_read_mat_complex(rw_h,f_name,mat)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h                         !< Handle
   character(*),         intent(in)    :: f_name                       !< File name
   complex(kind=r8),     intent(out)   :: mat(rw_h%n_lrow,rw_h%n_lcol) !< Matrix

   character*40, parameter :: caller = "elsi_read_mat_complex"

   if(rw_h%parallel_mode == MULTI_PROC) then
      call elsi_read_mat_complex_mp(rw_h,f_name,mat)
   else
      call elsi_read_mat_complex_sp(rw_h,f_name,mat)
   endif

end subroutine

!>
!! This routine writes a 2D block-cyclic dense matrix to file.
!!
subroutine elsi_write_mat_complex(rw_h,f_name,mat)

   implicit none

   type(elsi_rw_handle), intent(in) :: rw_h                         !< Handle
   character(*),         intent(in) :: f_name                       !< File name
   complex(kind=r8),     intent(in) :: mat(rw_h%n_lrow,rw_h%n_lcol) !< Matrix

   character*40, parameter :: caller = "elsi_write_mat_complex"

   if(rw_h%parallel_mode == MULTI_PROC) then
      call elsi_write_mat_complex_mp(rw_h,f_name,mat)
   else
      call elsi_write_mat_complex_sp(rw_h,f_name,mat)
   endif

end subroutine

!>
!! This routine reads the dimensions of a 2D block-cyclic dense matrix from
!! file.
!!
subroutine elsi_read_mat_dim_mp(rw_h,f_name,n_electron,n_basis,n_lrow,n_lcol)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h       !< Handle
   character(*),         intent(in)    :: f_name     !< File name
   real(kind=r8),        intent(out)   :: n_electron !< Number of electrons
   integer(kind=i4),     intent(out)   :: n_basis    !< Matrix size
   integer(kind=i4),     intent(out)   :: n_lrow     !< Local number of rows
   integer(kind=i4),     intent(out)   :: n_lcol     !< Local number of columns

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_prow
   integer(kind=i4) :: n_pcol
   integer(kind=i4) :: my_prow
   integer(kind=i4) :: my_pcol
   integer(kind=i8) :: offset
   logical          :: file_ok

   integer(kind=i4), external :: numroc

   character*40, parameter :: caller = "elsi_read_mat_dim_mp"

   inquire(file=f_name,exist=file_ok)

   if(.not. file_ok) then
      call elsi_rw_stop(" File does not exist.",rw_h,caller)
   endif

   ! Open file
   f_mode = mpi_mode_rdonly

   call MPI_File_open(rw_h%mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Read header
   if(rw_h%myid == 0) then
      offset = 0

      call MPI_File_read_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Broadcast header
   call MPI_Bcast(header,HEADER_SIZE,mpi_integer4,0,rw_h%mpi_comm,mpierr)

   n_basis    = header(4)
   n_electron = real(header(5),kind=r8)

   ! Get processor grid information
   call blacs_gridinfo(rw_h%blacs_ctxt,n_prow,n_pcol,my_prow,my_pcol)

   ! Get local size of matrix
   n_lrow = numroc(n_basis,rw_h%blk,my_prow,0,n_prow)
   n_lcol = numroc(n_basis,rw_h%blk,my_pcol,0,n_pcol)

   rw_h%n_basis     = n_basis
   rw_h%n_electrons = n_electron
   rw_h%n_lrow      = n_lrow
   rw_h%n_lcol      = n_lcol
   rw_h%nnz_g       = header(6)

end subroutine

!>
!! This routine reads the dimensions of a 1D block CSC matrix from file.
!!
subroutine elsi_read_mat_dim_sparse(rw_h,f_name,n_electron,n_basis,nnz_g,&
              nnz_l_sp,n_lcol_sp)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h       !< Handle
   character(*),         intent(in)    :: f_name     !< File name
   real(kind=r8),        intent(out)   :: n_electron !< Number of electrons
   integer(kind=i4),     intent(out)   :: n_basis    !< Matrix size
   integer(kind=i4),     intent(out)   :: nnz_g      !< Global number of nonzeros
   integer(kind=i4),     intent(out)   :: nnz_l_sp   !< Local number of nonzeros
   integer(kind=i4),     intent(out)   :: n_lcol_sp  !< Local number of columns

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_lcol0
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   logical          :: file_ok

   integer(kind=i4), allocatable :: col_ptr(:)

   character*40, parameter :: caller = "elsi_read_mat_dim_sparse"

   inquire(file=f_name,exist=file_ok)

   if(.not. file_ok) then
      call elsi_rw_stop(" File does not exist.",rw_h,caller)
   endif

   ! Open file
   f_mode = mpi_mode_rdonly

   call MPI_File_open(rw_h%mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Read header
   if(rw_h%myid == 0) then
      offset = 0

      call MPI_File_read_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Broadcast header
   call MPI_Bcast(header,HEADER_SIZE,mpi_integer4,0,rw_h%mpi_comm,mpierr)

   n_basis    = header(4)
   n_electron = real(header(5),kind=r8)
   nnz_g      = header(6)

   ! Compute n_lcol
   n_lcol_sp = n_basis/rw_h%n_procs
   n_lcol0   = n_lcol_sp
   if(rw_h%myid == rw_h%n_procs-1) then
      n_lcol_sp = n_basis-(rw_h%n_procs-1)*n_lcol0
   endif

   allocate(col_ptr(n_lcol_sp+1))

   ! Read column pointer
   offset = HEADER_SIZE*4+rw_h%myid*n_lcol0*4

   call MPI_File_read_at_all(f_handle,offset,col_ptr,n_lcol_sp+1,mpi_integer4,&
           mpi_status_ignore,mpierr)

   if(rw_h%myid == rw_h%n_procs-1) then
      col_ptr(n_lcol_sp+1) = nnz_g+1
   endif

   ! Shift column pointer
   prev_nnz = col_ptr(1)-1
   col_ptr  = col_ptr-prev_nnz

   ! Compute nnz_l_sp
   nnz_l_sp = col_ptr(n_lcol_sp+1)-col_ptr(1)

   deallocate(col_ptr)

   rw_h%n_basis     = n_basis
   rw_h%n_electrons = n_electron
   rw_h%nnz_g       = nnz_g
   rw_h%nnz_l_sp    = nnz_l_sp
   rw_h%n_lcol_sp   = n_lcol_sp

end subroutine

!>
!! This routine reads a 2D block-cyclic dense matrix from file.
!!
subroutine elsi_read_mat_real_mp(rw_h,f_name,mat)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h                         !< Handle
   character(*),         intent(in)    :: f_name                       !< File name
   real(kind=r8),        intent(out)   :: mat(rw_h%n_lrow,rw_h%n_lcol) !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_lcol0
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   logical          :: file_ok
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)
   real(kind=r8),    allocatable :: nnz_val(:)

   type(elsi_handle) :: aux_h

   character*40, parameter :: caller = "elsi_read_mat_real_mp"

   inquire(file=f_name,exist=file_ok)

   if(.not. file_ok) then
      call elsi_rw_stop(" File does not exist.",rw_h,caller)
   endif

   call elsi_init(aux_h,ELPA_SOLVER,MULTI_PROC,BLACS_DENSE,rw_h%n_basis,&
           rw_h%n_electrons,0)
   call elsi_set_mpi(aux_h,rw_h%mpi_comm)
   call elsi_set_blacs(aux_h,rw_h%blacs_ctxt,rw_h%blk)

   ! Output
   aux_h%myid_all   = rw_h%myid
   aux_h%print_info = rw_h%print_info
   aux_h%print_mem  = rw_h%print_mem
   aux_h%print_unit = rw_h%print_unit

   call elsi_get_time(aux_h,t0)

   ! Open file
   f_mode = mpi_mode_rdonly

   call MPI_File_open(rw_h%mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Compute n_lcol_sp
   rw_h%n_lcol_sp  = rw_h%n_basis/rw_h%n_procs
   n_lcol0         = rw_h%n_lcol_sp
   if(rw_h%myid == rw_h%n_procs-1) then
      rw_h%n_lcol_sp = rw_h%n_basis-(rw_h%n_procs-1)*n_lcol0
   endif

   call elsi_allocate(aux_h,col_ptr,rw_h%n_lcol_sp+1,"col_ptr",caller)

   ! Read column pointer
   offset = HEADER_SIZE*4+rw_h%myid*n_lcol0*4

   call MPI_File_read_at_all(f_handle,offset,col_ptr,rw_h%n_lcol_sp+1,&
           mpi_integer4,mpi_status_ignore,mpierr)

   if(rw_h%myid == rw_h%n_procs-1) then
      col_ptr(rw_h%n_lcol_sp+1) = rw_h%nnz_g+1
   endif

   ! Shift column pointer
   prev_nnz = col_ptr(1)-1
   col_ptr  = col_ptr-prev_nnz

   ! Compute nnz_l_sp
   rw_h%nnz_l_sp = col_ptr(rw_h%n_lcol_sp+1)-col_ptr(1)

   call elsi_allocate(aux_h,row_ind,rw_h%nnz_l_sp,"row_ind",caller)

   ! Read row index
   offset = HEADER_SIZE*4+rw_h%n_basis*4+prev_nnz*4

   call MPI_File_read_at_all(f_handle,offset,row_ind,rw_h%nnz_l_sp,&
           mpi_integer4,mpi_status_ignore,mpierr)

   ! Read non-zero value
   offset = HEADER_SIZE*4+rw_h%n_basis*4+rw_h%nnz_g*4+prev_nnz*8

   call elsi_allocate(aux_h,nnz_val,rw_h%nnz_l_sp,"nnz_val",caller)

   call MPI_File_read_at_all(f_handle,offset,nnz_val,rw_h%nnz_l_sp,mpi_real8,&
           mpi_status_ignore,mpierr)

   ! Close file
   call MPI_File_close(f_handle,mpierr)

   ! Redistribute matrix
   call elsi_set_csc(aux_h,rw_h%nnz_g,rw_h%nnz_l_sp,rw_h%n_lcol_sp,row_ind,&
           col_ptr)
   call elsi_set_sparse_dm(aux_h,nnz_val)

   aux_h%n_elsi_calls  = 1
   aux_h%my_prow_pexsi = 0
   aux_h%np_per_pole   = rw_h%n_procs

   call elsi_pexsi_to_blacs_dm(aux_h,mat)

   call elsi_deallocate(aux_h,col_ptr,"col_ptr")
   call elsi_deallocate(aux_h,row_ind,"row_ind")
   call elsi_deallocate(aux_h,nnz_val,"nnz_val")

   call elsi_get_time(aux_h,t1)

   write(info_str,"('  Finished reading matrix')")
   call elsi_say(info_str,aux_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,aux_h)

   call elsi_cleanup(aux_h)

end subroutine

!>
!! This routine reads a 1D block CSC matrix from file.
!!
subroutine elsi_read_mat_real_sparse(rw_h,f_name,row_ind,col_ptr,mat)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h                      !< Handle
   character(*),         intent(in)    :: f_name                    !< File name
   integer(kind=i4),     intent(out)   :: row_ind(rw_h%nnz_l_sp)    !< Row index
   integer(kind=i4),     intent(out)   :: col_ptr(rw_h%n_lcol_sp+1) !< Column pointer
   real(kind=r8),        intent(out)   :: mat(rw_h%nnz_l_sp)        !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: n_lcol0
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   logical          :: file_ok
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   type(elsi_handle) :: aux_h

   character*40, parameter :: caller = "elsi_read_mat_real_sparse"

   inquire(file=f_name,exist=file_ok)

   if(.not. file_ok) then
      call elsi_rw_stop(" File does not exist.",rw_h,caller)
   endif

   ! Output
   aux_h%myid_all   = rw_h%myid
   aux_h%print_info = rw_h%print_info
   aux_h%print_unit = rw_h%print_unit

   call elsi_init_timer(aux_h)
   call elsi_get_time(aux_h,t0)

   ! Open file
   f_mode = mpi_mode_rdonly

   call MPI_File_open(rw_h%mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Compute n_lcol0
   n_lcol0 = rw_h%n_basis/rw_h%n_procs

   ! Read column pointer
   offset = HEADER_SIZE*4+rw_h%myid*n_lcol0*4

   call MPI_File_read_at_all(f_handle,offset,col_ptr,rw_h%n_lcol_sp+1,&
           mpi_integer4,mpi_status_ignore,mpierr)

   if(rw_h%myid == rw_h%n_procs-1) then
      col_ptr(rw_h%n_lcol_sp+1) = rw_h%nnz_g+1
   endif

   ! Shift column pointer
   prev_nnz = col_ptr(1)-1
   col_ptr  = col_ptr-prev_nnz

   ! Read row index
   offset = HEADER_SIZE*4+rw_h%n_basis*4+prev_nnz*4

   call MPI_File_read_at_all(f_handle,offset,row_ind,rw_h%nnz_l_sp,&
           mpi_integer4,mpi_status_ignore,mpierr)

   ! Read non-zero value
   offset = HEADER_SIZE*4+rw_h%n_basis*4+rw_h%nnz_g*4+prev_nnz*8

   call MPI_File_read_at_all(f_handle,offset,mat,rw_h%nnz_l_sp,mpi_real8,&
           mpi_status_ignore,mpierr)

   call MPI_File_close(f_handle,mpierr)

   call elsi_get_time(aux_h,t1)

   write(info_str,"('  Finished reading matrix')")
   call elsi_say(info_str,aux_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,aux_h)

end subroutine

!>
!! This routine reads a 2D block-cyclic dense matrix from file.
!!
subroutine elsi_read_mat_complex_mp(rw_h,f_name,mat)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h                         !< Handle
   character(*),         intent(in)    :: f_name                       !< File name
   complex(kind=r8),     intent(out)   :: mat(rw_h%n_lrow,rw_h%n_lcol) !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_lcol0
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   logical          :: file_ok
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)
   complex(kind=r8), allocatable :: nnz_val(:)

   type(elsi_handle) :: aux_h

   character*40, parameter :: caller = "elsi_read_mat_complex_mp"

   inquire(file=f_name,exist=file_ok)

   if(.not. file_ok) then
      call elsi_rw_stop(" File does not exist.",rw_h,caller)
   endif

   call elsi_init(aux_h,ELPA_SOLVER,MULTI_PROC,BLACS_DENSE,rw_h%n_basis,&
           rw_h%n_electrons,0)
   call elsi_set_mpi(aux_h,rw_h%mpi_comm)
   call elsi_set_blacs(aux_h,rw_h%blacs_ctxt,rw_h%blk)

   ! Output
   aux_h%myid_all   = rw_h%myid
   aux_h%print_info = rw_h%print_info
   aux_h%print_mem  = rw_h%print_mem
   aux_h%print_unit = rw_h%print_unit

   call elsi_get_time(aux_h,t0)

   ! Open file
   f_mode = mpi_mode_rdonly

   call MPI_File_open(rw_h%mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Compute n_lcol_sp
   rw_h%n_lcol_sp = rw_h%n_basis/rw_h%n_procs
   n_lcol0        = rw_h%n_lcol_sp
   if(rw_h%myid == rw_h%n_procs-1) then
      rw_h%n_lcol_sp = rw_h%n_basis-(rw_h%n_procs-1)*n_lcol0
   endif

   call elsi_allocate(aux_h,col_ptr,rw_h%n_lcol_sp+1,"col_ptr",caller)

   ! Read column pointer
   offset = HEADER_SIZE*4+rw_h%myid*n_lcol0*4

   call MPI_File_read_at_all(f_handle,offset,col_ptr,rw_h%n_lcol_sp+1,&
           mpi_integer4,mpi_status_ignore,mpierr)

   if(rw_h%myid == rw_h%n_procs-1) then
      col_ptr(rw_h%n_lcol_sp+1) = rw_h%nnz_g+1
   endif

   ! Shift column pointer
   prev_nnz = col_ptr(1)-1
   col_ptr  = col_ptr-prev_nnz

   ! Compute nnz_l_sp
   rw_h%nnz_l_sp = col_ptr(rw_h%n_lcol_sp+1)-col_ptr(1)

   call elsi_allocate(aux_h,row_ind,rw_h%nnz_l_sp,"row_ind",caller)

   ! Read row index
   offset = HEADER_SIZE*4+rw_h%n_basis*4+prev_nnz*4

   call MPI_File_read_at_all(f_handle,offset,row_ind,rw_h%nnz_l_sp,&
           mpi_integer4,mpi_status_ignore,mpierr)

   ! Read non-zero value
   offset = HEADER_SIZE*4+rw_h%n_basis*4+rw_h%nnz_g*4+prev_nnz*16

   call elsi_allocate(aux_h,nnz_val,rw_h%nnz_l_sp,"nnz_val",caller)

   call MPI_File_read_at_all(f_handle,offset,nnz_val,rw_h%nnz_l_sp,&
           mpi_complex16,mpi_status_ignore,mpierr)

   ! Close file
   call MPI_File_close(f_handle,mpierr)

   ! Redistribute matrix
   call elsi_set_csc(aux_h,rw_h%nnz_g,rw_h%nnz_l_sp,rw_h%n_lcol_sp,row_ind,&
           col_ptr)
   call elsi_set_sparse_dm(aux_h,nnz_val)

   aux_h%n_elsi_calls  = 1
   aux_h%my_prow_pexsi = 0
   aux_h%np_per_pole   = rw_h%n_procs

   call elsi_pexsi_to_blacs_dm(aux_h,mat)

   call elsi_deallocate(aux_h,col_ptr,"col_ptr")
   call elsi_deallocate(aux_h,row_ind,"row_ind")
   call elsi_deallocate(aux_h,nnz_val,"nnz_val")

   call elsi_get_time(aux_h,t1)

   write(info_str,"('  Finished reading matrix')")
   call elsi_say(info_str,aux_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,aux_h)

   call elsi_cleanup(aux_h)

end subroutine

!>
!! This routine reads a 1D block CSC matrix from file.
!!
subroutine elsi_read_mat_complex_sparse(rw_h,f_name,row_ind,col_ptr,mat)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h                      !< Handle
   character(*),         intent(in)    :: f_name                    !< File name
   integer(kind=i4),     intent(out)   :: row_ind(rw_h%nnz_l_sp)    !< Row index
   integer(kind=i4),     intent(out)   :: col_ptr(rw_h%n_lcol_sp+1) !< Column pointer
   complex(kind=r8),     intent(out)   :: mat(rw_h%nnz_l_sp)        !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: n_lcol0
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   logical          :: file_ok
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   type(elsi_handle) :: aux_h

   character*40, parameter :: caller = "elsi_read_mat_complex_sparse"

   inquire(file=f_name,exist=file_ok)

   if(.not. file_ok) then
      call elsi_rw_stop(" File does not exist.",rw_h,caller)
   endif

   ! Output
   aux_h%myid_all   = rw_h%myid
   aux_h%print_info = rw_h%print_info
   aux_h%print_unit = rw_h%print_unit

   call elsi_init_timer(aux_h)
   call elsi_get_time(aux_h,t0)

   ! Open file
   f_mode = mpi_mode_rdonly

   call MPI_File_open(rw_h%mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Compute n_lcol0
   n_lcol0 = rw_h%n_basis/rw_h%n_procs

   ! Read column pointer
   offset = HEADER_SIZE*4+rw_h%myid*n_lcol0*4

   call MPI_File_read_at_all(f_handle,offset,col_ptr,rw_h%n_lcol_sp+1,&
           mpi_integer4,mpi_status_ignore,mpierr)

   if(rw_h%myid == rw_h%n_procs-1) then
      col_ptr(rw_h%n_lcol_sp+1) = rw_h%nnz_g+1
   endif

   ! Shift column pointer
   prev_nnz = col_ptr(1)-1
   col_ptr  = col_ptr-prev_nnz

   ! Read row index
   offset = HEADER_SIZE*4+rw_h%n_basis*4+prev_nnz*4

   call MPI_File_read_at_all(f_handle,offset,row_ind,rw_h%nnz_l_sp,&
           mpi_integer4,mpi_status_ignore,mpierr)

   ! Read non-zero value
   offset = HEADER_SIZE*4+rw_h%n_basis*4+rw_h%nnz_g*4+prev_nnz*16

   call MPI_File_read_at_all(f_handle,offset,mat,rw_h%nnz_l_sp,mpi_complex16,&
           mpi_status_ignore,mpierr)

   call MPI_File_close(f_handle,mpierr)

   call elsi_get_time(aux_h,t1)

   write(info_str,"('  Finished reading matrix')")
   call elsi_say(info_str,aux_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,aux_h)

   call elsi_cleanup(aux_h)

end subroutine

!>
!! This routine writes a 2D block-cyclic dense matrix to file.
!!
subroutine elsi_write_mat_real_mp(rw_h,f_name,mat)

   implicit none

   type(elsi_rw_handle), intent(in) :: rw_h                         !< Handle
   character(*),         intent(in) :: f_name                       !< File name
   real(kind=r8),        intent(in) :: mat(rw_h%n_lrow,rw_h%n_lcol) !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_lcol0
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   real(kind=r8), allocatable :: dummy(:,:)

   type(elsi_handle) :: aux_h

   character*40, parameter :: caller = "elsi_write_mat_real_mp"

   call elsi_init(aux_h,ELPA_SOLVER,MULTI_PROC,BLACS_DENSE,rw_h%n_basis,&
           rw_h%n_electrons,0)
   call elsi_set_mpi(aux_h,rw_h%mpi_comm)
   call elsi_set_blacs(aux_h,rw_h%blacs_ctxt,rw_h%blk)
   call elsi_get_time(aux_h,t0)

   ! Output
   aux_h%myid_all   = rw_h%myid
   aux_h%print_info = rw_h%print_info
   aux_h%print_mem  = rw_h%print_mem
   aux_h%print_unit = rw_h%print_unit

   aux_h%zero_def     = rw_h%zero_def
   aux_h%ovlp_is_unit = .true.
   aux_h%n_elsi_calls = 1
   aux_h%n_lcol_sp    = rw_h%n_basis/rw_h%n_procs
   if(rw_h%myid == rw_h%n_procs-1) then
      aux_h%n_lcol_sp = rw_h%n_basis-(rw_h%n_procs-1)*aux_h%n_lcol_sp
   endif

   call elsi_allocate(aux_h,dummy,1,1,"dummy",caller)
   call elsi_blacs_to_sips_hs(aux_h,mat,dummy)
   call elsi_deallocate(aux_h,dummy,"dummy")

   ! Open file
   f_mode = mpi_mode_wronly+mpi_mode_create

   call MPI_File_open(rw_h%mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Write header
   header(1)    = FILE_VERSION
   header(2)    = rw_h%file_format
   header(3)    = REAL_VALUES
   header(4)    = rw_h%n_basis
   header(5)    = int(rw_h%n_electrons,kind=i4)
   header(6)    = aux_h%nnz_g
   header(7:8)  = 0
   header(9:16) = rw_h%header_user

   if(rw_h%myid == 0) then
      offset = 0

      call MPI_File_write_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Compute shift of column pointers
   prev_nnz = 0

   call MPI_Exscan(aux_h%nnz_l_sp,prev_nnz,1,mpi_integer4,mpi_sum,&
           rw_h%mpi_comm,mpierr)

   ! Shift column pointer
   aux_h%col_ptr_sips = aux_h%col_ptr_sips+prev_nnz

   ! Write column pointer
   n_lcol0 = rw_h%n_basis/rw_h%n_procs
   offset  = HEADER_SIZE*4+rw_h%myid*n_lcol0*4

   call MPI_File_write_at_all(f_handle,offset,aux_h%col_ptr_sips,&
           aux_h%n_lcol_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Write row index
   offset = HEADER_SIZE*4+rw_h%n_basis*4+prev_nnz*4

   call MPI_File_write_at_all(f_handle,offset,aux_h%row_ind_sips,&
           aux_h%nnz_l_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Write non-zero value
   offset = HEADER_SIZE*4+rw_h%n_basis*4+aux_h%nnz_g*4+prev_nnz*8

   call MPI_File_write_at_all(f_handle,offset,aux_h%ham_real_sips,&
           aux_h%nnz_l_sp,mpi_real8,mpi_status_ignore,mpierr)

   call MPI_File_close(f_handle,mpierr)

   call elsi_get_time(aux_h,t1)

   write(info_str,"('  Finished writing matrix')")
   call elsi_say(info_str,aux_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,aux_h)

   call elsi_cleanup(aux_h)

end subroutine

!>
!! This routine writes a 2D block-cyclic dense matrix to file.
!!
subroutine elsi_write_mat_complex_mp(rw_h,f_name,mat)

   implicit none

   type(elsi_rw_handle), intent(in) :: rw_h                         !< Handle
   character(*),         intent(in) :: f_name                       !< File name
   complex(kind=r8),     intent(in) :: mat(rw_h%n_lrow,rw_h%n_lcol) !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: n_lcol0
   integer(kind=i4) :: prev_nnz
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   complex(kind=r8), allocatable :: dummy(:,:)

   type(elsi_handle) :: aux_h

   character*40, parameter :: caller = "elsi_write_mat_complex_mp"

   call elsi_init(aux_h,ELPA_SOLVER,MULTI_PROC,BLACS_DENSE,rw_h%n_basis,&
           rw_h%n_electrons,0)
   call elsi_set_mpi(aux_h,rw_h%mpi_comm)
   call elsi_set_blacs(aux_h,rw_h%blacs_ctxt,rw_h%blk)
   call elsi_get_time(aux_h,t0)

   ! Output
   aux_h%myid_all   = rw_h%myid
   aux_h%print_info = rw_h%print_info
   aux_h%print_mem  = rw_h%print_mem
   aux_h%print_unit = rw_h%print_unit

   aux_h%zero_def     = rw_h%zero_def
   aux_h%ovlp_is_unit = .true.
   aux_h%n_elsi_calls = 1
   aux_h%n_lcol_sp    = rw_h%n_basis/rw_h%n_procs
   if(rw_h%myid == rw_h%n_procs-1) then
      aux_h%n_lcol_sp = rw_h%n_basis-(rw_h%n_procs-1)*aux_h%n_lcol_sp
   endif

   call elsi_allocate(aux_h,dummy,1,1,"dummy",caller)
   call elsi_blacs_to_sips_hs(aux_h,mat,dummy)
   call elsi_deallocate(aux_h,dummy,"dummy")

   ! Open file
   f_mode = mpi_mode_wronly+mpi_mode_create

   call MPI_File_open(rw_h%mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Write header
   header(1)    = FILE_VERSION
   header(2)    = rw_h%file_format
   header(3)    = COMPLEX_VALUES
   header(4)    = rw_h%n_basis
   header(5)    = int(rw_h%n_electrons,kind=i4)
   header(6)    = aux_h%nnz_g
   header(7:8)  = 0
   header(9:16) = rw_h%header_user

   if(rw_h%myid == 0) then
      offset = 0

      call MPI_File_write_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Compute shift of column pointers
   prev_nnz = 0

   call MPI_Exscan(aux_h%nnz_l_sp,prev_nnz,1,mpi_integer4,mpi_sum,&
           rw_h%mpi_comm,mpierr)

   ! Shift column pointer
   aux_h%col_ptr_sips = aux_h%col_ptr_sips+prev_nnz

   ! Write column pointer
   n_lcol0 = rw_h%n_basis/rw_h%n_procs
   offset  = HEADER_SIZE*4+rw_h%myid*n_lcol0*4

   call MPI_File_write_at_all(f_handle,offset,aux_h%col_ptr_sips,&
           aux_h%n_lcol_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Write row index
   offset = HEADER_SIZE*4+rw_h%n_basis*4+prev_nnz*4

   call MPI_File_write_at_all(f_handle,offset,aux_h%row_ind_sips,&
           aux_h%nnz_l_sp,mpi_integer4,mpi_status_ignore,mpierr)

   ! Write non-zero value
   offset = HEADER_SIZE*4+rw_h%n_basis*4+aux_h%nnz_g*4+prev_nnz*16

   call MPI_File_write_at_all(f_handle,offset,aux_h%ham_cmplx_sips,&
           aux_h%nnz_l_sp,mpi_complex16,mpi_status_ignore,mpierr)

   call MPI_File_close(f_handle,mpierr)

   call elsi_get_time(aux_h,t1)

   write(info_str,"('  Finished writing matrix')")
   call elsi_say(info_str,aux_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,aux_h)

   call elsi_cleanup(aux_h)

end subroutine

!>
!! This routine writes a 1D block CSC matrix to file.
!!
subroutine elsi_write_mat_real_sparse(rw_h,f_name,row_ind,col_ptr,mat)

   implicit none

   type(elsi_rw_handle), intent(in) :: rw_h                      !< Handle
   character(*),         intent(in) :: f_name                    !< File name
   integer(kind=i4),     intent(in) :: row_ind(rw_h%nnz_l_sp)    !< Row index
   integer(kind=i4),     intent(in) :: col_ptr(rw_h%n_lcol_sp+1) !< Column pointer
   real(kind=r8),        intent(in) :: mat(rw_h%nnz_l_sp)        !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: prev_nnz
   integer(kind=i4) :: n_lcol0
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), allocatable :: col_ptr_shift(:)

   type(elsi_handle) :: aux_h

   character*40, parameter :: caller = "elsi_write_mat_real_sparse"

   ! Output
   aux_h%myid_all   = rw_h%myid
   aux_h%print_info = rw_h%print_info
   aux_h%print_unit = rw_h%print_unit

   call elsi_init_timer(aux_h)
   call elsi_get_time(aux_h,t0)

   ! Open file
   f_mode = mpi_mode_wronly+mpi_mode_create

   call MPI_File_open(rw_h%mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Write header
   header(1)    = FILE_VERSION
   header(2)    = rw_h%file_format
   header(3)    = REAL_VALUES
   header(4)    = rw_h%n_basis
   header(5)    = int(rw_h%n_electrons,kind=i4)
   header(6)    = rw_h%nnz_g
   header(7:8)  = 0
   header(9:16) = rw_h%header_user

   if(rw_h%myid == 0) then
      offset = 0

      call MPI_File_write_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Compute shift of column pointers
   prev_nnz = 0

   call MPI_Exscan(rw_h%nnz_l_sp,prev_nnz,1,mpi_integer4,mpi_sum,rw_h%mpi_comm,&
           mpierr)

   ! Shift column pointer
   allocate(col_ptr_shift(rw_h%n_lcol_sp+1))

   col_ptr_shift = col_ptr+prev_nnz

   ! Write column pointer
   n_lcol0 = rw_h%n_basis/rw_h%n_procs
   offset  = HEADER_SIZE*4+rw_h%myid*n_lcol0*4

   call MPI_File_write_at_all(f_handle,offset,col_ptr_shift,rw_h%n_lcol_sp,&
           mpi_integer4,mpi_status_ignore,mpierr)

   deallocate(col_ptr_shift)

   ! Write row index
   offset = HEADER_SIZE*4+rw_h%n_basis*4+prev_nnz*4

   call MPI_File_write_at_all(f_handle,offset,row_ind,rw_h%nnz_l_sp,&
           mpi_integer4,mpi_status_ignore,mpierr)

   ! Write non-zero value
   offset = HEADER_SIZE*4+rw_h%n_basis*4+rw_h%nnz_g*4+prev_nnz*8

   call MPI_File_write_at_all(f_handle,offset,mat,rw_h%nnz_l_sp,mpi_real8,&
           mpi_status_ignore,mpierr)

   call MPI_File_close(f_handle,mpierr)

   call elsi_get_time(aux_h,t1)

   write(info_str,"('  Finished writing matrix')")
   call elsi_say(info_str,aux_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,aux_h)

end subroutine

!>
!! This routine writes a 1D block CSC matrix to file.
!!
subroutine elsi_write_mat_complex_sparse(rw_h,f_name,row_ind,col_ptr,mat)

   implicit none

   type(elsi_rw_handle), intent(in) :: rw_h                      !< Handle
   character(*),         intent(in) :: f_name                    !< File name
   integer(kind=i4),     intent(in) :: row_ind(rw_h%nnz_l_sp)    !< Row index
   integer(kind=i4),     intent(in) :: col_ptr(rw_h%n_lcol_sp+1) !< Column pointer
   complex(kind=r8),     intent(in) :: mat(rw_h%nnz_l_sp)        !< Matrix

   integer(kind=i4) :: mpierr
   integer(kind=i4) :: f_handle
   integer(kind=i4) :: f_mode
   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: prev_nnz
   integer(kind=i4) :: n_lcol0
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), allocatable :: col_ptr_shift(:)

   type(elsi_handle) :: aux_h

   character*40, parameter :: caller = "elsi_write_mat_complex_sparse"

   ! Output
   aux_h%myid_all   = rw_h%myid
   aux_h%print_info = rw_h%print_info
   aux_h%print_unit = rw_h%print_unit

   call elsi_init_timer(aux_h)
   call elsi_get_time(aux_h,t0)

   ! Open file
   f_mode = mpi_mode_wronly+mpi_mode_create

   call MPI_File_open(rw_h%mpi_comm,f_name,f_mode,mpi_info_null,f_handle,mpierr)

   ! Write header
   header(1)    = FILE_VERSION
   header(2)    = rw_h%file_format
   header(3)    = COMPLEX_VALUES
   header(4)    = rw_h%n_basis
   header(5)    = int(rw_h%n_electrons,kind=i4)
   header(6)    = rw_h%nnz_g
   header(7:8)  = 0
   header(9:16) = rw_h%header_user

   if(rw_h%myid == 0) then
      offset = 0

      call MPI_File_write_at(f_handle,offset,header,HEADER_SIZE,mpi_integer4,&
              mpi_status_ignore,mpierr)
   endif

   ! Compute shift of column pointers
   prev_nnz = 0

   call MPI_Exscan(rw_h%nnz_l_sp,prev_nnz,1,mpi_integer4,mpi_sum,rw_h%mpi_comm,&
           mpierr)

   ! Shift column pointer
   allocate(col_ptr_shift(rw_h%n_lcol_sp+1))

   col_ptr_shift = col_ptr+prev_nnz

   ! Write column pointer
   n_lcol0 = rw_h%n_basis/rw_h%n_procs
   offset  = HEADER_SIZE*4+rw_h%myid*n_lcol0*4

   call MPI_File_write_at_all(f_handle,offset,col_ptr_shift,rw_h%n_lcol_sp,&
           mpi_integer4,mpi_status_ignore,mpierr)

   deallocate(col_ptr_shift)

   ! Write row index
   offset = HEADER_SIZE*4+rw_h%n_basis*4+prev_nnz*4

   call MPI_File_write_at_all(f_handle,offset,row_ind,rw_h%nnz_l_sp,&
           mpi_integer4,mpi_status_ignore,mpierr)

   ! Write non-zero value
   offset = HEADER_SIZE*4+rw_h%n_basis*4+rw_h%nnz_g*4+prev_nnz*16

   call MPI_File_write_at_all(f_handle,offset,mat,rw_h%nnz_l_sp,mpi_complex16,&
           mpi_status_ignore,mpierr)

   call MPI_File_close(f_handle,mpierr)

   call elsi_get_time(aux_h,t1)

   write(info_str,"('  Finished writing matrix')")
   call elsi_say(info_str,aux_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,aux_h)

end subroutine

!>
!! This routine reads the dimensions of a 2D block-cyclic dense matrix from
!! file.
!!
subroutine elsi_read_mat_dim_sp(rw_h,f_name,n_electron,n_basis,n_lrow,n_lcol)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h       !< Handle
   character(*),         intent(in)    :: f_name     !< File name
   real(kind=r8),        intent(out)   :: n_electron !< Number of electrons
   integer(kind=i4),     intent(out)   :: n_basis    !< Matrix size
   integer(kind=i4),     intent(out)   :: n_lrow     !< Local number of rows
   integer(kind=i4),     intent(out)   :: n_lcol     !< Local number of columns

   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i8) :: offset
   logical          :: file_ok

   character*40, parameter :: caller = "elsi_read_mat_dim_mp"

   inquire(file=f_name,exist=file_ok)

   if(.not. file_ok) then
      call elsi_rw_stop(" File does not exist.",rw_h,caller)
   endif

   ! Open file
   open(file=f_name,unit=99,access="stream",form="unformatted")

   ! Read header
   offset = 1

   read(unit=99,pos=offset) header

   close(unit=99)

   n_basis    = header(4)
   n_electron = real(header(5),kind=r8)
   n_lrow     = n_basis
   n_lcol     = n_basis

   rw_h%n_basis     = n_basis
   rw_h%n_electrons = n_electron
   rw_h%n_lrow      = n_lrow
   rw_h%n_lcol      = n_lcol

end subroutine

!>
!! This routine reads a 2D block-cyclic dense matrix from file.
!!
subroutine elsi_read_mat_real_sp(rw_h,f_name,mat)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h                         !< Handle
   character(*),         intent(in)    :: f_name                       !< File name
   real(kind=r8),        intent(out)   :: mat(rw_h%n_lrow,rw_h%n_lcol) !< Matrix

   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: this_nnz
   integer(kind=i8) :: offset
   logical          :: file_ok
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)
   real(kind=r8),    allocatable :: nnz_val(:)

   type(elsi_handle) :: aux_h

   character*40, parameter :: caller = "elsi_read_mat_real_sp"

   inquire(file=f_name,exist=file_ok)

   if(.not. file_ok) then
      call elsi_rw_stop(" File does not exist.",rw_h,caller)
   endif

   ! Output
   aux_h%myid_all   = rw_h%myid
   aux_h%print_info = rw_h%print_info
   aux_h%print_mem  = rw_h%print_mem
   aux_h%print_unit = rw_h%print_unit

   call elsi_init_timer(aux_h)
   call elsi_get_time(aux_h,t0)

   ! Open file
   open(file=f_name,unit=99,access="stream",form="unformatted")

   ! Read header
   offset = 1

   read(unit=99,pos=offset) header

   rw_h%nnz_g     = header(6)
   rw_h%nnz_l_sp  = header(6)
   rw_h%n_lcol_sp = rw_h%n_basis

   call elsi_allocate(aux_h,col_ptr,rw_h%n_lcol_sp+1,"col_ptr",caller)

   ! Read column pointer
   offset = 1+HEADER_SIZE*4

   read(unit=99,pos=offset) col_ptr(1:rw_h%n_lcol_sp)

   col_ptr(rw_h%n_lcol_sp+1) = rw_h%nnz_g+1

   call elsi_allocate(aux_h,row_ind,rw_h%nnz_l_sp,"row_ind",caller)

   ! Read row index
   offset = 1+HEADER_SIZE*4+rw_h%n_basis*4

   read(unit=99,pos=offset) row_ind

   ! Read non-zero value
   offset = 1+HEADER_SIZE*4+rw_h%n_basis*4+rw_h%nnz_g*4

   call elsi_allocate(aux_h,nnz_val,rw_h%nnz_l_sp,"nnz_val",caller)

   read(unit=99,pos=offset) nnz_val

   ! Close file
   close(unit=99)

   ! Convert to dense
   mat = 0.0_r8

   i_val = 0
   do i = 1,rw_h%n_basis
      this_nnz = col_ptr(i+1)-col_ptr(i)

      do j = i_val+1,i_val+this_nnz
         mat(row_ind(j),i) = nnz_val(j)
      enddo

      i_val = i_val+this_nnz
   enddo

   call elsi_deallocate(aux_h,col_ptr,"col_ptr")
   call elsi_deallocate(aux_h,row_ind,"row_ind")
   call elsi_deallocate(aux_h,nnz_val,"nnz_val")

   call elsi_get_time(aux_h,t1)

   write(info_str,"('  Finished reading matrix')")
   call elsi_say(info_str,aux_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,aux_h)

end subroutine

!>
!! This routine reads a 2D block-cyclic dense matrix from file.
!!
subroutine elsi_read_mat_complex_sp(rw_h,f_name,mat)

   implicit none

   type(elsi_rw_handle), intent(inout) :: rw_h                         !< Handle
   character(*),         intent(in)    :: f_name                       !< File name
   complex(kind=r8),     intent(out)   :: mat(rw_h%n_lrow,rw_h%n_lcol) !< Matrix

   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: this_nnz
   integer(kind=i8) :: offset
   logical          :: file_ok
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)
   complex(kind=r8), allocatable :: nnz_val(:)

   type(elsi_handle) :: aux_h

   character*40, parameter :: caller = "elsi_read_mat_complex_sp"

   inquire(file=f_name,exist=file_ok)

   if(.not. file_ok) then
      call elsi_rw_stop(" File does not exist.",rw_h,caller)
   endif

   ! Output
   aux_h%myid_all   = rw_h%myid
   aux_h%print_info = rw_h%print_info
   aux_h%print_unit = rw_h%print_unit

   call elsi_init_timer(aux_h)
   call elsi_get_time(aux_h,t0)

   ! Open file
   open(file=f_name,unit=99,access="stream",form="unformatted")

   ! Read header
   offset = 1

   read(unit=99,pos=offset) header

   rw_h%nnz_g     = header(6)
   rw_h%nnz_l_sp  = header(6)
   rw_h%n_lcol_sp = rw_h%n_basis

   call elsi_allocate(aux_h,col_ptr,rw_h%n_lcol_sp+1,"col_ptr",caller)

   ! Read column pointer
   offset = 1+HEADER_SIZE*4

   read(unit=99,pos=offset) col_ptr(1:rw_h%n_lcol_sp)

   col_ptr(rw_h%n_lcol_sp+1) = rw_h%nnz_g+1

   call elsi_allocate(aux_h,row_ind,rw_h%nnz_l_sp,"row_ind",caller)

   ! Read row index
   offset = 1+HEADER_SIZE*4+rw_h%n_basis*4

   read(unit=99,pos=offset) row_ind

   ! Read non-zero value
   offset = 1+HEADER_SIZE*4+rw_h%n_basis*4+rw_h%nnz_g*4

   call elsi_allocate(aux_h,nnz_val,rw_h%nnz_l_sp,"nnz_val",caller)

   read(unit=99,pos=offset) nnz_val

   ! Close file
   close(unit=99)

   ! Convert to dense
   mat = (0.0_r8,0.0_r8)

   i_val = 0
   do i = 1,rw_h%n_basis
      this_nnz = col_ptr(i+1)-col_ptr(i)

      do j = i_val+1,i_val+this_nnz
         mat(row_ind(j),i) = nnz_val(j)
      enddo

      i_val = i_val+this_nnz
   enddo

   call elsi_deallocate(aux_h,col_ptr,"col_ptr")
   call elsi_deallocate(aux_h,row_ind,"row_ind")
   call elsi_deallocate(aux_h,nnz_val,"nnz_val")

   call elsi_get_time(aux_h,t1)

   write(info_str,"('  Finished reading matrix')")
   call elsi_say(info_str,aux_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,aux_h)

end subroutine

!>
!! This routine writes a 2D block-cyclic dense matrix to file.
!!
subroutine elsi_write_mat_real_sp(rw_h,f_name,mat)

   implicit none

   type(elsi_rw_handle), intent(in) :: rw_h                         !< Handle
   character(*),         intent(in) :: f_name                       !< File name
   real(kind=r8),        intent(in) :: mat(rw_h%n_lrow,rw_h%n_lcol) !< Matrix

   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: this_nnz
   integer(kind=i4) :: nnz_g
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)
   real(kind=r8),    allocatable :: nnz_val(:)

   type(elsi_handle) :: aux_h

   character*40, parameter :: caller = "elsi_write_mat_real_sp"

   ! Output
   aux_h%myid_all   = rw_h%myid
   aux_h%print_info = rw_h%print_info
   aux_h%print_unit = rw_h%print_unit
   aux_h%zero_def   = rw_h%zero_def

   call elsi_init_timer(aux_h)
   call elsi_get_time(aux_h,t0)

   ! Compute nnz
   call elsi_get_local_nnz(aux_h,mat,rw_h%n_lrow,rw_h%n_lcol,nnz_g)

   ! Convert to CSC
   call elsi_allocate(aux_h,col_ptr,rw_h%n_basis+1,"col_ptr",caller)
   call elsi_allocate(aux_h,row_ind,nnz_g,"row_ind",caller)
   call elsi_allocate(aux_h,nnz_val,nnz_g,"nnz_val",caller)

   i_val   = 0
   col_ptr = 1

   do i = 1,rw_h%n_lcol
      this_nnz = 0

      do j = 1,rw_h%n_lrow
         if(abs(mat(j,i)) > rw_h%zero_def) then
            this_nnz       = this_nnz+1
            i_val          = i_val+1
            nnz_val(i_val) = mat(j,i)
            row_ind(i_val) = j
         endif
      enddo

      col_ptr(i+1) = col_ptr(i)+this_nnz
   enddo

   ! Open file
   open(file=f_name,unit=99,access="stream",form="unformatted")

   ! Write header
   header(1)    = FILE_VERSION
   header(2)    = rw_h%file_format
   header(3)    = REAL_VALUES
   header(4)    = rw_h%n_basis
   header(5)    = int(rw_h%n_electrons,kind=i4)
   header(6)    = nnz_g
   header(7:8)  = 0
   header(9:16) = rw_h%header_user

   offset = 1
   write(unit=99,pos=offset) header

   ! Write column pointer
   offset = 1+HEADER_SIZE*4

   write(unit=99,pos=offset) col_ptr(1:rw_h%n_basis)

   ! Write row index
   offset = 1+HEADER_SIZE*4+rw_h%n_basis*4

   write(unit=99,pos=offset) row_ind

   ! Write non-zero value
   offset = 1+HEADER_SIZE*4+rw_h%n_basis*4+nnz_g*4

   write(unit=99,pos=offset) nnz_val

   ! Close file
   close(unit=99)

   call elsi_deallocate(aux_h,col_ptr,"col_ptr")
   call elsi_deallocate(aux_h,row_ind,"row_ind")
   call elsi_deallocate(aux_h,nnz_val,"nnz_val")

   call elsi_get_time(aux_h,t1)

   write(info_str,"('  Finished writing matrix')")
   call elsi_say(info_str,aux_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,aux_h)

end subroutine

!>
!! This routine writes a 2D block-cyclic dense matrix to file.
!!
subroutine elsi_write_mat_complex_sp(rw_h,f_name,mat)

   implicit none

   type(elsi_rw_handle), intent(in) :: rw_h                         !< Handle
   character(*),         intent(in) :: f_name                       !< File name
   complex(kind=r8),     intent(in) :: mat(rw_h%n_lrow,rw_h%n_lcol) !< Matrix

   integer(kind=i4) :: header(HEADER_SIZE)
   integer(kind=i4) :: i_val
   integer(kind=i4) :: i
   integer(kind=i4) :: j
   integer(kind=i4) :: this_nnz
   integer(kind=i4) :: nnz_g
   integer(kind=i8) :: offset
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   character*200    :: info_str

   integer(kind=i4), allocatable :: row_ind(:)
   integer(kind=i4), allocatable :: col_ptr(:)
   complex(kind=r8), allocatable :: nnz_val(:)

   type(elsi_handle) :: aux_h

   character*40, parameter :: caller = "elsi_write_mat_complex_sp"

   ! Output
   aux_h%myid_all   = rw_h%myid
   aux_h%print_info = rw_h%print_info
   aux_h%print_unit = rw_h%print_unit
   aux_h%zero_def   = rw_h%zero_def

   call elsi_init_timer(aux_h)
   call elsi_get_time(aux_h,t0)

   ! Compute nnz
   call elsi_get_local_nnz(aux_h,mat,rw_h%n_lrow,rw_h%n_lcol,nnz_g)

   ! Convert to CSC
   call elsi_allocate(aux_h,col_ptr,rw_h%n_basis+1,"col_ptr",caller)
   call elsi_allocate(aux_h,row_ind,nnz_g,"row_ind",caller)
   call elsi_allocate(aux_h,nnz_val,nnz_g,"nnz_val",caller)

   i_val   = 0
   col_ptr = 1

   do i = 1,rw_h%n_lcol
      this_nnz = 0

      do j = 1,rw_h%n_lrow
         if(abs(mat(j,i)) > rw_h%zero_def) then
            this_nnz       = this_nnz+1
            i_val          = i_val+1
            nnz_val(i_val) = mat(j,i)
            row_ind(i_val) = j
         endif
      enddo

      col_ptr(i+1) = col_ptr(i)+this_nnz
   enddo

   ! Open file
   open(file=f_name,unit=99,access="stream",form="unformatted")

   ! Write header
   header(1)    = FILE_VERSION
   header(2)    = rw_h%file_format
   header(3)    = COMPLEX_VALUES
   header(4)    = rw_h%n_basis
   header(5)    = int(rw_h%n_electrons,kind=i4)
   header(6)    = nnz_g
   header(7:8)  = 0
   header(9:16) = rw_h%header_user

   offset = 1
   write(unit=99,pos=offset) header

   ! Write column pointer
   offset = 1+HEADER_SIZE*4

   write(unit=99,pos=offset) col_ptr(1:rw_h%n_basis)

   ! Write row index
   offset = 1+HEADER_SIZE*4+rw_h%n_basis*4

   write(unit=99,pos=offset) row_ind

   ! Write non-zero value
   offset = 1+HEADER_SIZE*4+rw_h%n_basis*4+nnz_g*4

   write(unit=99,pos=offset) nnz_val

   ! Close file
   close(unit=99)

   call elsi_deallocate(aux_h,col_ptr,"col_ptr")
   call elsi_deallocate(aux_h,row_ind,"row_ind")
   call elsi_deallocate(aux_h,nnz_val,"nnz_val")

   call elsi_get_time(aux_h,t1)

   write(info_str,"('  Finished writing matrix')")
   call elsi_say(info_str,aux_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(info_str,aux_h)

end subroutine

end module ELSI_IO
