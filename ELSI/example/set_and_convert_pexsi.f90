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
!! This program illustrates how to construct and write an eigenvalue problem
!! using the ELSI Interface software
!!

program set_and_convert

  use iso_c_binding
  use ELSI

  implicit none
  include 'mpif.h'

  ! This is the ELSI test suite
  ! First we will test the writing and reading of a matrix to a file

  ! Ok, let us first create a matrix for the ELSI interface
  ! Something like

  integer :: matrixsize = 10, blocksize = 2
 
  ! Local variables
  real*8  :: element
  real*8, allocatable  :: H_matrix(:,:)
  real*8, allocatable  :: S_matrix(:,:)
  real*8, allocatable  :: R_matrix(:,:)

  ! Local MPI vars
  integer :: mpierr
  integer :: myid
  integer :: root
  integer :: i_proc
  integer :: n_procs

  ! Local BLACS vars
  integer,external :: numroc
  integer :: blacs_ctxt
  integer :: blacs_info
  integer :: sc_desc(9)
  integer :: mpi_comm_row
  integer :: mpi_comm_col
  integer :: n_process_rows
  integer :: n_process_cols
  integer :: my_process_row
  integer :: my_process_col
  ! Grid 2
  integer :: sc_desc_2(9)
  integer :: block_rows_2
  integer :: block_cols_2


  ! Local Matrix dimensions
  integer :: n_rows
  integer :: n_cols
  integer :: n_rows_2
  integer :: n_cols_2
  integer :: n_rows_Pexsi
  integer :: n_cols_Pexsi
  integer :: n_cols_add
  integer :: my_col_offset
  integer :: cols_sent

  ! Positions
  integer :: global_row, local_row
  integer :: global_col, local_col
  integer :: block, offset

  ! Helpers
  integer :: offset_sent
  integer :: offset_receive
  integer :: my_counts
  integer, allocatable :: strides(:)
  integer, allocatable :: rcounts(:)
  integer :: share 


  !  Pharse command line argumnents, if given
   INTEGER*4 :: iargc
   character*16 arg1
   character*16 arg2

   if (iargc() == 2) then
      call getarg(1, arg1)
      call getarg(2, arg2)
      read(arg1, *) matrixsize
      read(arg2, *) blocksize
   endif

  ! Simulate external MPI environment
  call mpi_init(mpierr)
  if (mpierr /= 0) then
      write(*,'(a)') "MPI: Failed to initialize MPI."
      stop
  end if
  call mpi_comm_rank(mpi_comm_world, myid, mpierr)
  if (mpierr /= 0) then
      write(*,'(a)') "MPI: Failed to initialize MPI."
      stop
   end if
  call mpi_comm_size(mpi_comm_world, n_procs, mpierr)
  if (mpierr /= 0) then
      write(*,'(a)') "MPI: Failed to initialize MPI."
      stop
   end if

  ! Simulate external blacs environment
  do n_process_cols = NINT(SQRT(REAL(n_procs))),n_procs,1
     if(mod(n_procs,n_process_cols) == 0 ) exit
  enddo

  n_process_rows = n_procs / n_process_cols

  blacs_ctxt = mpi_comm_world
  call BLACS_Gridinit( blacs_ctxt, 'R', n_process_rows, n_process_cols )
  call BLACS_Gridinfo( blacs_ctxt, n_process_rows, n_process_cols, &
        my_process_row, my_process_col )
  call mpi_comm_split(mpi_comm_world,my_process_col,my_process_row,&
        mpi_comm_row, mpierr)
  call mpi_comm_split(mpi_comm_world,my_process_row,my_process_col,&
        mpi_comm_col, mpierr)
  n_rows = numroc(matrixsize, blocksize, my_process_row, 0, n_process_rows)
  n_cols = numroc(matrixsize, blocksize, my_process_col, 0, n_process_cols)
  call descinit( sc_desc, matrixsize, matrixsize, blocksize, blocksize, 0, 0, &
                 blacs_ctxt, MAX(1,n_rows), blacs_info )


  ! First set the parallel treatment
  call elsi_set_mpi(mpi_comm_world,n_procs,myid)
  
  ! Second ELSI Specifications
  call elsi_set_method(PEXSI)
  call elsi_set_mode(REAL_VALUES)
  
  ! Third Define the problem
  call elsi_initialize_problem(matrixsize, blocksize, blocksize)
  
  ! Forth Set the parallel distribution
  call elsi_initialize_blacs() 

  ! Simulate external matrix setup
  allocate(H_matrix(n_rows,n_cols))
  allocate(S_matrix(n_rows,n_cols))
  H_matrix = 0d0
  S_matrix = 0d0

  do local_row = 1, n_rows
    block = FLOOR( 1d0 * (local_row - 1) / blocksize)
    offset = local_row - block * blocksize
    global_row = my_process_row * blocksize + block * blocksize * n_process_rows + offset
    do local_col = 1, n_cols
      block = FLOOR( 1d0 * (local_col - 1) / blocksize)
      offset = local_col - block * blocksize
      global_col = my_process_col * blocksize + block * blocksize * n_process_cols + offset

      element = 1d3 * global_row + 1d0 * global_col
      H_matrix(local_row, local_col) = element
      S_matrix(local_row, local_col) = element + 0.5d0
    end do
  end do

  call scalapack_dense_to_pexsi_sparse( H_matrix, S_matrix, n_rows, n_cols,&
        mpi_comm_world, blacs_ctxt, sc_desc)

  ! Write eigenvalue problem to another file
  call elsi_write_ev_problem("elsi_eigenvalue_problem_out.hdf5")
  
  ! elsi shutdown
  call elsi_finalize()

  deallocate(H_matrix)
  deallocate(S_matrix)

  call blacs_gridexit(blacs_ctxt)
  call mpi_finalize(mpierr)

end program
