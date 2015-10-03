!Copyright (c) 2015, ELSI consortium
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
! * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer in the
!   documentation and/or other materials provided with the distribution.
! * Neither the name of the ELSI project nor the
!   names of its contributors may be used to endorse or promote products
!   derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!> 
!! This module contains wrapper functions to the MPI library 
!!

module ELSI_MPI_TOOLS

  use iso_c_binding

  use ELSI_DIMENSIONS
  use matrixswitch

  implicit none
  private

  !> Calculates the local matrix dimensions based on the BLACS grid
  integer, external :: numroc

  public :: elsi_initialize_mpi
  public :: elsi_initialize_blacs 
  public :: elsi_set_mpi
  public :: elsi_set_blacs 
  public :: elsi_get_global_row
  public :: elsi_get_global_col
  public :: elsi_get_global_dimensions
  public :: elsi_get_local_dimensions
  public :: elsi_get_processor_grid
  public :: elsi_get_comm_grid
  public :: elsi_get_myid
  public :: elsi_finalize_mpi
  public :: elsi_finalize_blacs
  public :: elsi_stop

  contains

!> 
!! Initialize MPI
!!
subroutine elsi_initialize_mpi()

   implicit none
   include 'mpif.h'

   ! MPI Initialization

   external_mpi = .False.

   call mpi_init(mpierr)
   if (mpierr /= 0) then
      write(*,'(a)') "MPI: Failed to initialize MPI."
      stop
   end if
   call mpi_comm_rank(mpi_comm_world, myid, mpierr)
   if (mpierr /= 0) then
      write(*,'(a)') "MPI: Failed to determine MPI rank."
      stop
   end if

   mpi_comm_global = mpi_comm_world

   call mpi_comm_size(mpi_comm_global, n_procs, mpierr)
   if (mpierr /= 0) then
      write(*,'(a)') "MPI: Failed to determine MPI size."
      stop
   end if

end subroutine

!> 
!! Set MPI from external
!!
subroutine elsi_set_mpi(global_comm, n_procs_in, myid_in)

   implicit none
   include 'mpif.h'

   integer,intent(in) :: myid_in            !< local process id
   integer,intent(in) :: n_procs_in         !< number of mpi processes
   integer,intent(in) :: global_comm !< global mpi communicator


   ! MPI Initialization

   external_mpi = .True.
   mpi_comm_global = global_comm
   n_procs         = n_procs_in
   myid            = myid_in
end subroutine


!> 
!! Finalize MPI
!!
subroutine elsi_finalize_mpi()

   implicit none
   include 'mpif.h'

   ! MPI Finalization

   call mpi_finalize(mpierr)
   if (mpierr /= 0) then
      write(*,'(a)') "MPI: Failed to finalize MPI."
      stop
   end if

end subroutine


!> 
!! Initialize BLACS Grid
!!
subroutine elsi_initialize_blacs()

   implicit none

   include "mpif.h"

   ! Exception how to deviate from blocksize for OMM
   integer :: exception(2) = (/0,0/)

   external_blacs = .False.
  ! Define blockcyclic setup
  do n_p_cols = NINT(SQRT(REAL(n_procs))),2,-1
     if(mod(n_procs,n_p_cols) == 0 ) exit
  enddo

  n_p_rows = n_procs / n_p_cols

  ! Set up BLACS and MPI communicators

  blacs_ctxt = mpi_comm_global
  call BLACS_Gridinit( blacs_ctxt, 'C', n_p_rows, n_p_cols )
  call BLACS_Gridinfo( blacs_ctxt, n_p_rows, n_p_cols, my_p_row, my_p_col )

  call mpi_comm_split(mpi_comm_global,my_p_col,my_p_row,mpi_comm_row,mpierr)

   if (mpierr /= 0) then
      write(*,'(a)') "ELPA: Failed to get ELPA row communicators."
      stop
   end if
 
  call mpi_comm_split(mpi_comm_global,my_p_row,my_p_col,mpi_comm_col,mpierr)

  if (mpierr /= 0) then
      write(*,'(a)') "ELPA: Failed to get ELPA column communicators."
      stop
   end if

  n_l_rows = numroc(n_g_rank, n_b_rows, my_p_row, 0, n_p_rows)
  n_l_cols = numroc(n_g_rank, n_b_cols, my_p_col, 0, n_p_cols)

  if(myid == 0) print *, 'Global matrixsize: ',n_g_rank, ' x ', n_g_rank
  if(myid == 0) print *, 'Blocksize: ',n_b_rows, ' x ', n_b_cols
  if(myid == 0) print *, 'Processor grid: ',n_p_rows, ' x ', n_p_cols
  if(myid == 0) print *, 'Local Matrixsize: ',n_l_rows, ' x ', n_l_cols

  call descinit( sc_desc, n_g_rank, n_g_rank, n_b_rows, n_b_cols, 0, 0, &
                 blacs_ctxt, MAX(1,n_l_rows), blacs_info )

  if (method == OMM_DENSE) then
     call ms_scalapack_setup (n_procs, n_p_rows, "C", n_b_rows, exception,&
           blacs_ctxt)
  end if


end subroutine

!> 
!! Set BLACS Grid from external
!!
subroutine elsi_set_blacs( blacs_ctxt_in, n_p_rows_in, n_p_cols_in, &
      my_p_row_in, my_p_col_in, mpi_comm_row_in, mpi_comm_col_in,&
      n_l_rows_in, n_l_cols_in, sc_desc_in)

   implicit none

   include "mpif.h"
   ! Exception how to deviate from blocksize for OMM
   integer :: exception(2) = (/0,0/)

  integer,intent(in) :: blacs_ctxt_in     !< local blacs context
  integer,intent(in) :: sc_desc_in(9)     !< local blacs context
  integer,intent(in) :: mpi_comm_row_in   !< row    communicatior
  integer,intent(in) :: mpi_comm_col_in   !< column communicatior
  integer,intent(in) :: my_p_row_in       !< process row    position
  integer,intent(in) :: my_p_col_in       !< process column position
  integer,intent(in) :: n_p_rows_in       !< Number of processes in row
  integer,intent(in) :: n_p_cols_in       !< Number of processes in column
  integer,intent(in) :: n_l_rows_in       !< Number of local rows
  integer,intent(in) :: n_l_cols_in       !< Number of local columns
   
  external_blacs = .True.
  blacs_ctxt = blacs_ctxt_in 
  n_p_rows = n_p_rows_in 
  n_p_cols = n_p_cols_in
  my_p_row = my_p_row_in
  my_p_col = my_p_col_in
  mpi_comm_row = mpi_comm_row_in
  mpi_comm_col = mpi_comm_col_in
  n_l_rows = n_l_rows_in
  n_l_cols = n_l_cols_in
  sc_desc = sc_desc_in

  if(myid == 0) print *, 'Global matrixsize: ',n_g_rank, ' x ', n_g_rank
  if(myid == 0) print *, 'Blocksize: ',n_b_rows, ' x ', n_b_cols
  if(myid == 0) print *, 'Processor grid: ',n_p_rows, ' x ', n_p_cols
  if(myid == 0) print *, 'Local Matrixsize: ',n_l_rows, ' x ', n_l_cols
  call MPI_BARRIER(mpi_comm_global,mpierr)

  print *, myid," Id mpi_comm_global: ", mpi_comm_global
  print *, myid," Id ip_row/col: ", my_p_row, "/", my_p_col
  print *, myid," Id MPI_COMM_rows/cols: ", mpi_comm_row, "/", mpi_comm_col
  call MPI_BARRIER(mpi_comm_global,mpierr)
  
  if (method == OMM_DENSE) then
     call ms_scalapack_setup (n_procs, n_p_rows, "C", n_b_rows, exception,&
           blacs_ctxt)
  end if

end subroutine


!> 
!! Finalize BLACS grid
!!
subroutine elsi_finalize_blacs()

   implicit none

   call blacs_gridexit(blacs_ctxt)

end subroutine

!> 
!! Computes global row index based on local row index
!!
subroutine elsi_get_global_row (global_idx, local_idx)

   implicit none

   integer, intent(out) :: global_idx  !< Global index 
   integer, intent(in)  :: local_idx   !< Local index

   integer :: block !< local block 
   integer :: idx   !< local index in block


   block = FLOOR( 1d0 * (local_idx - 1) / n_b_rows) 
   idx = local_idx - block * n_b_rows

   global_idx = my_p_row * n_b_rows + block * n_b_rows * n_p_rows + idx
  
end subroutine


!> 
!! Computes global column index based on local column index
!!
subroutine elsi_get_global_col (global_idx, local_idx)

   implicit none

   integer, intent(out) :: global_idx   !< global index
   integer, intent(in)  :: local_idx    !< local index

   integer :: block !< local block
   integer :: idx   !< local index in block


   block = FLOOR( 1d0 * (local_idx - 1) / n_b_cols) 
   idx = local_idx - block * n_b_cols

   global_idx = my_p_col * n_b_cols + block * n_p_cols * n_b_cols + idx
  
end subroutine

!> 
!! Gets global matrix specifications such as rank and block dimensions
!!
subroutine elsi_get_global_dimensions (g_dim, g_block_rows, g_block_cols)

   implicit none

   integer, intent(out) :: g_dim   !< global rank
   integer, intent(out) :: g_block_rows !< global blocksize
   integer, intent(out) :: g_block_cols !< global blocksize

   g_dim   = n_g_rank 
   g_block_rows = n_b_rows
   g_block_cols = n_b_cols

end subroutine

!> 
!! Gets the local matrix dimensions
!!
subroutine elsi_get_local_dimensions (rows, cols)

   implicit none

   integer, intent(out) :: rows   !< output local rows
   integer, intent(out) :: cols   !< output local cols

   rows  = n_l_rows 
   cols  = n_l_cols

end subroutine

!> 
!! Get the processor grid as specified by BLACS
!!
subroutine elsi_get_processor_grid (rows, cols)

   implicit none

   integer, intent(out) :: rows   !< output processor rows
   integer, intent(out) :: cols   !< output processor cols

   rows  = n_p_rows 
   cols  = n_p_cols

end subroutine


!> 
!! Gets the row and column communicators for ELPA
!!
subroutine elsi_get_comm_grid (rows, cols)

   implicit none

   integer, intent(out) :: rows   !< output row communicator
   integer, intent(out) :: cols   !< output col communicator

   rows  = mpi_comm_row 
   cols  = mpi_comm_col

end subroutine

!> 
!! Gets the process id 
!!

subroutine elsi_get_myid (id)

   implicit none

   integer, intent(out) :: id   !< output process id

   id = myid 

end subroutine

!>
!! Clean shutdown in case of error
!!
subroutine elsi_stop(message, caller)

      implicit none

      character(len=*), intent(in) :: message
      character(len=*), intent(in) :: caller

      character(LEN=4096) :: string_message

      write(string_message, "(1X,'*** Proc',I5,' in ',A,': ',A)") &
           & myid, trim(caller), trim(message)

      write(*,'(A)') trim(string_message)
        
      if (n_procs > 1) then
         call MPI_Abort(mpi_comm_global, 0, mpierr)
      end if
        
      stop

end subroutine elsi_stop

end module ELSI_MPI_TOOLS
