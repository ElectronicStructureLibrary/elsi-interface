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
  real*8, allocatable  :: buffer(:,:)
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
  do n_process_cols = NINT(SQRT(REAL(n_procs))),2,-1
     if(mod(n_procs,n_process_cols) == 0 ) exit
  enddo

  n_process_rows = n_procs / n_process_cols

  blacs_ctxt = mpi_comm_world
  call BLACS_Gridinit( blacs_ctxt, 'C', n_process_rows, n_process_cols )
  call BLACS_Gridinfo( blacs_ctxt, n_process_rows, n_process_cols, &
        my_process_row, my_process_col )
  call mpi_comm_split(mpi_comm_world,my_process_col,my_process_row,&
        mpi_comm_row,mpierr)
  call mpi_comm_split(mpi_comm_world,my_process_row,my_process_col,&
        mpi_comm_col,mpierr)
  n_rows = numroc(matrixsize, blocksize, my_process_row, 0, n_process_rows)
  n_cols = numroc(matrixsize, blocksize, my_process_col, 0, n_process_cols)
  call descinit( sc_desc, matrixsize, matrixsize, blocksize, blocksize, 0, 0, &
                 blacs_ctxt, MAX(1,n_rows), blacs_info )

  if ((n_process_rows == 1 .or. n_process_cols == 1) .and. n_procs /= 1) then
     if (myid == 0) print *, "We encountered an scalapack bug when working "//&
     "with prime process numbers and using pdtran to transform a matrix "  //&
     "to a different blacs transcriptor. We stop here and wait for a "     //&
     "scalapack patch. While this setup is an inefficient corner case, "   //&
     "restart your calculation choosing a process number which in a "      //&
     "square number in best case for optimal performance and refrain from "//&
     "using prime numbers."
     call MPI_ABORT(mpi_comm_world)
     stop
  end if

  ! First set the parallel treatment
  call elsi_set_mpi(mpi_comm_world,n_procs,myid)
  
  ! Second ELSI Specifications
  call elsi_set_method(ELPA)
  call elsi_set_mode(REAL_VALUES)
  
  ! Third Define the problem
  call elsi_initialize_problem(matrixsize, blocksize, blocksize)
  
  ! Forth Set the parallel distribution
  call elsi_set_blacs(blacs_ctxt, n_process_rows, n_process_cols,&
        my_process_row, my_process_col, mpi_comm_row, mpi_comm_col, &
        n_rows, n_cols, sc_desc)
  call MPI_BARRIER(mpi_comm_world,mpierr)

  ! Simulate external matrix setup
  allocate(H_matrix(n_rows,n_cols))
  H_matrix = 0d0

  do local_row = 1, n_rows
    call elsi_get_global_row(global_row, local_row)
    do local_col = 1, n_cols
      call elsi_get_global_col(global_col, local_col)
      element = 1d3 * global_row + 1d0 * global_col
      H_matrix(local_row, local_col) = element
    end do
  end do

  if (myid == 0) print *, "H_matrix"
  call MPI_BARRIER(mpi_comm_world,mpierr)
  do i_proc = 0, n_procs - 1
    if (myid == i_proc) print *, myid, ": ", H_matrix
    call MPI_BARRIER(mpi_comm_world,mpierr)
  end do

  ! Now set up the other grid
  ! Simulate external blacs environment
  block_rows_2 = matrixsize
  block_cols_2 = FLOOR(1d0 * matrixsize / n_process_cols)

  n_rows_2 = matrixsize
  n_cols_2 = numroc(matrixsize, block_cols_2, my_process_col, 0, n_process_cols)
  call descinit( sc_desc_2, matrixsize, matrixsize, block_rows_2, block_cols_2,&
        0, 0, blacs_ctxt, MAX(1,n_rows_2), blacs_info )

  allocate(buffer(n_rows_2,n_cols_2))
  buffer = 0d0

  ! Hamiltonian is symmetric, this is a way rearrange the memory kayou
  call pdtran(matrixsize, matrixsize, &
        1d0, H_matrix, 1, 1, sc_desc, &
        0d0, buffer, 1, 1, sc_desc_2)

  ! However, no inter blacs is possible so we have to do some by hand communication
  call MPI_Bcast(buffer, n_rows_2 * n_cols_2, MPI_DOUBLE, 0, mpi_comm_row,&
        mpierr) 

  ! PEXSI setup
  n_rows_Pexsi = matrixsize
  n_cols_Pexsi = Floor(1d0 * matrixsize / n_procs)
  n_cols_add = matrixsize - n_cols_Pexsi * n_procs

  ! The last process gathers all remaining columns
  if (myid == n_procs - 1) then
     allocate(R_matrix(n_rows_Pexsi,n_cols_Pexsi+n_cols_add))
  else
     allocate(R_matrix(n_rows_Pexsi,n_cols_Pexsi))
  end if
  R_matrix = 0d0
     
  ! Get my id within the the current process row
  call mpi_comm_rank(mpi_comm_row, my_col_offset, mpierr)
  my_col_offset = 1 + my_col_offset * Floor(1d0 * matrixsize / n_procs)

  ! Fill elements from the buffer
  R_matrix(:,1:n_cols_Pexsi) = &
     buffer(:,my_col_offset:my_col_offset+n_cols_Pexsi-1)

  ! Only the last process in the last process row gathers the remaining columns.
  if (n_cols_add > 0) then 
     
     ! How many columns do I have to share?
     cols_sent = n_cols_2 - n_process_rows * n_cols_Pexsi

     ! If I receive elements where do I start to put them?
     if (myid .eq. n_procs - 1) then
        offset_receive = n_cols_Pexsi + 1
     else 
        offset_receive = 1
     end if

     ! If I sent elements where do I start to take them?
     if (cols_sent > 0) then
        offset_sent = n_cols_2 - cols_sent + 1
     else 
        offset_sent = 1
     end if
    
     allocate(strides(n_process_cols))
     allocate(rcounts(n_process_cols))

     strides = 0
     rcounts = 0

     ! How many elements do I sent
     my_counts = cols_sent * n_rows_2
     
     ! Tell Host process how many elements I sent  
     call MPI_Gather (my_counts, 1, MPI_INT, &
           rcounts, 1, MPI_INT, &
           n_process_cols - 1, mpi_comm_col, mpierr)

     ! Where will I put the recieved elements from process i
     do i_proc = 2, n_process_cols
        strides(i_proc) = SUM(rcounts(1:i_proc-1))
     end do

     if (my_process_row .eq. n_process_rows - 1) then     
       call MPI_GatherV (buffer(1,offset_sent), &
            n_rows_2 * cols_sent, MPI_DOUBLE, &
            R_matrix(1,offset_receive), &
            rcounts, strides, &
            MPI_DOUBLE, n_process_cols - 1, mpi_comm_col, mpierr)
     end if

     deallocate(strides)
     deallocate(rcounts)
  
  end if

  if (myid == 0) print *, "R_matrix"
  call MPI_BARRIER(mpi_comm_world,mpierr)
  do i_proc = 0, n_procs - 1
    if (myid == i_proc) print *, myid, ": ", R_matrix
    call MPI_BARRIER(mpi_comm_world,mpierr)
  end do

  deallocate(H_matrix)
  deallocate(buffer)
  deallocate(R_matrix) 

  call blacs_gridexit(blacs_ctxt)
  call mpi_finalize(mpierr)

end program
