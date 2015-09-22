!> ELSI Interchange
!! 
!! Copyright of the original code rests with the authors inside the ELSI
!! consortium. The copyright of any additional modifications shall rest
!! with their original authors, but shall adhere to the licensing terms
!! distributed along with the original code in the file "COPYING".

program readwrite

  use iso_c_binding
  use ELSI
  use MPI_TOOLS

  implicit none
  include 'mpif.h'

  ! This is the ELSI test suite
  ! First we will test the writing and reading of a matrix to a file

  ! Ok, let us first create a matrix for the ELSI interface
  ! Something like

  integer :: matrixsize = 10, blocksize = 2
 
  ! Random number generator
  integer :: iseed(4096) !< Random seed

  ! Test matrices
  real*8, allocatable :: matrix_h(:,:)
  real*8, allocatable :: matrix_s(:,:)

  ! Local variables
  real*8 :: errmax
  real*8 :: colerr
  integer :: i
  integer :: myid
  integer :: n_procs
  integer :: error
  integer :: i_row, i_col, p_row, p_col, l_row, l_col, g_row, g_col


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
  call elsi_initialize_mpi(myid)
  call elsi_initialize_problem(matrixsize, blocksize)
  call elsi_set_method(ELPA)
  call elsi_set_mode(REAL_VALUES)

  !do i_row = 1, matrixsize
  !  my_task = elsi_get_local_row(i_row, p_row, l_row)
  !  call elsi_get_global_row(g_row, p_row, l_row)
  !  do i_col = 1, matrixsize
  !    my_task = elsi_get_local_col(i_col, p_col, l_col)
  !    call elsi_get_global_col(g_col, p_col, l_col)
  !    if(myid==0) write(*,'(8(a,I3),a)') &
  !      '(',i_row,',',i_col,') mapped to process (', &
  !          p_row,',',p_col,') and local (',l_row,',',l_col,') &
  !      mapped back to (', g_row,',',g_col,')'
  !   end do
  !end do

  ! Generate some data
  allocate(matrix_h(n_rows,n_cols))

  !call RANDOM_NUMBER(buffer)
  !h = buffer

  ! Processor matrix
  do l_row = 1, n_rows
    do l_col = 1, n_cols
      matrix_h(l_row,l_col) = 1.d0 * myid
    end do
  end do

  allocate(matrix_s(n_rows,n_cols))
  do l_row = 1, n_rows
    call elsi_get_global_row(i_row, p_row, l_row)
    do l_col = 1, n_cols
      call elsi_get_global_col(i_col, p_col, l_col)
      if (i_row == i_col) then
         matrix_s(l_row,l_col) = 1.d0
      else 
         matrix_s(l_row,l_col) = 0.d0
      end if
    end do
  end do


  ! Initialize the data space
  call elsi_allocate_matrices()

  ! Set matrices
  call elsi_set_hamiltonian(matrix_h, n_rows, n_cols)
  call elsi_set_overlap(matrix_s, n_rows, n_cols)

  ! Write eigenvalue problem
  call elsi_write_ev_problem("elsi_eigenvalue_problem.hdf5")

  ! Deallocations
  deallocate(matrix_h) 
  deallocate(matrix_s) 

  ! elsi shutdown
  call elsi_finalize()

end program
