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

program construct_and_write

  use iso_c_binding
  use ELSI
  use ELSI_MPI_TOOLS

  implicit none
  include 'mpif.h'

  ! This is the ELSI test suite
  ! First we will test the writing and reading of a matrix to a file

  ! Ok, let us first create a matrix for the ELSI interface
  ! Something like

  integer :: matrixsize = 10, blocksize = 2
 
  ! Random number generator
  integer :: iseed(4096) !< Random seed

  ! Local variables
  real*8  :: element
  integer :: myid
  integer :: n_rows, n_cols, my_p_row, my_p_col
  integer :: l_row, l_col, i_row, i_col


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

  ! Now set some ELSI specifications
  call elsi_initialize_mpi()
  call elsi_initialize_problem(matrixsize, blocksize, blocksize)
  call elsi_initialize_blacs()
  call elsi_set_method(ELPA)
  call elsi_set_mode(REAL_VALUES)
  
   ! Initialize the data space
  call elsi_allocate_matrices()

  call elsi_get_local_dimensions(n_rows,n_cols)

  call elsi_get_myid(myid)
  iseed(:) = myid + 1 
  call RANDOM_SEED(put=iseed)

  ! Construct H and S
  do l_row = 1, n_rows
    call elsi_get_global_row(i_row, l_row)
    do l_col = l_row, n_cols
      call elsi_get_global_col(i_col, l_col)
      call RANDOM_NUMBER(element)
      if (i_row == i_col) then
         call elsi_set_hamiltonian_element(element, l_row, l_col)
         call elsi_set_overlap_element    (1.0d0,   l_row, l_col)
      else 
         call elsi_set_hamiltonian_element(element, l_row, l_col)
         call elsi_set_overlap_element    (0.0d0,   l_row, l_col)
      end if
    end do
  end do

  call elsi_symmetrize_hamiltonian()
  call elsi_symmetrize_overlap()

  ! Write eigenvalue problem
  call elsi_write_ev_problem("elsi_eigenvalue_problem.hdf5")

  ! elsi shutdown
  call elsi_finalize()

end program
