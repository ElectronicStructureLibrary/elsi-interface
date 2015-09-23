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
  call elsi_initialize_mpi(myid)
  call elsi_initialize_problem(matrixsize, blocksize, blocksize)
  call elsi_set_method(ELPA)
  call elsi_set_mode(REAL_VALUES)
  
   ! Initialize the data space
  call elsi_allocate_matrices()

  call elsi_get_local_dimensions(n_rows,n_cols)

  iseed(:) = myid + 1 
  call RANDOM_SEED(put=iseed)

  ! Construct H
  do l_row = 1, n_rows
    call elsi_get_global_row(i_row, my_p_row, l_row)
    do l_col = l_row, n_cols
      call elsi_get_global_col(i_col, my_p_col, l_col)
      call RANDOM_NUMBER(element)
      if (i_row == i_col) then
         call elsi_set_overlap_element(0.5d0 * element, l_row, l_col)
      else 
         call elsi_set_overlap_element(element, l_row, l_col)
      end if
    end do
  end do

  call elsi_symmetrize_hamiltonian()

  ! Construct S = 1
  do l_row = 1, n_rows
    call elsi_get_global_row(i_row, my_p_row, l_row)
    do l_col = l_row, n_cols
      call elsi_get_global_col(i_col, my_p_col, l_col)
      if (i_row == i_col) then
         call elsi_set_overlap_element(0.5d0, l_row, l_col)
      else 
         call elsi_set_overlap_element(0.0d0, l_row, l_col)
      end if
    end do
  end do
  
  call elsi_symmetrize_overlap()

  ! Write eigenvalue problem
  call elsi_write_ev_problem("elsi_eigenvalue_problem.hdf5")

  ! elsi shutdown
  call elsi_finalize()

end program
