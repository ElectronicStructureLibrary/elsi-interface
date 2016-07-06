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

program set_and_write

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

  ! Local MPI vars
  integer :: mpierr
  integer :: myid
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

  ! Local Matrix dimensions
  integer :: n_rows
  integer :: n_cols

  ! Positions
  integer :: global_row, local_row
  integer :: global_col, local_col


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
  call mpi_comm_rank(mpi_comm_world, myid, mpierr)
  call mpi_comm_size(mpi_comm_world, n_procs, mpierr)

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
        mpi_comm_row,mpierr)
  call mpi_comm_split(mpi_comm_world,my_process_row,my_process_col,&
        mpi_comm_col,mpierr)
  n_rows = numroc(matrixsize, blocksize, my_process_row, 0, n_process_rows)
  n_cols = numroc(matrixsize, blocksize, my_process_col, 0, n_process_cols)
  call descinit( sc_desc, matrixsize, matrixsize, blocksize, blocksize, 0, 0, &
                 blacs_ctxt, MAX(1,n_rows), blacs_info )

  ! First set the parallel treatment
  call elsi_set_mpi(mpi_comm_world,n_procs,myid)
  
  ! Second ELSI Specifications
  call elsi_set_method(LIBOMM)
  call elsi_set_mode(REAL_VALUES)
  
  ! Third define the problem
  call elsi_init_problem(matrixsize, blocksize, blocksize)
  
  ! Forth: Set the parallel distribution
  call elsi_set_blacs(blacs_ctxt, 'R', blocksize, blocksize, n_process_rows, &
        n_process_cols, my_process_row, my_process_col, n_rows, n_cols, &
        sc_desc, mpi_comm_row, mpi_comm_col)

  ! Simulate external matrix setup
  allocate(H_matrix(n_rows,n_cols))
  allocate(S_matrix(n_rows,n_cols))
  H_matrix = 0d0
  S_matrix = 0d0 
  do local_row = 1, n_rows
    call elsi_get_global_row(global_row, local_row)
    do local_col = 1, n_cols
      call elsi_get_global_col(global_col, local_col)
      if (global_row >= global_col) then
        element = 1d3 * global_row + 1d0 * global_col
        H_matrix(local_row,local_col) = element
        if (global_row == global_col) then
           S_matrix(local_row,local_col) = 1d0
        else 
           S_matrix(local_row,local_col) = 0d0
        end if
      end if
    end do
  end do

  ! Construct H and S
  call elsi_set_hamiltonian(H_matrix)
  call elsi_set_overlap(S_matrix)
  call elsi_symmetrize_hamiltonian()
  call elsi_symmetrize_overlap()

  ! Write eigenvalue problem
  call elsi_write_evp("elsi_eigenvalue_problem.hdf5")

  ! elsi shutdown
  call elsi_finalize()

  deallocate(H_matrix)
  deallocate(S_matrix)

  call blacs_gridexit(blacs_ctxt)
  call mpi_finalize(mpierr)

end program
