!> ELSI Interchange
!! 
!! Copyright of the original code rests with the authors inside the ELSI
!! consortium. The copyright of any additional modifications shall rest
!! with their original authors, but shall adhere to the licensing terms
!! distributed along with the original code in the file "COPYING".

program readwrite

  use iso_c_binding
  use ELSI
  use ELSI_MPI_TOOLS

  implicit none
  include 'mpif.h'

  ! This is the ELSI test suite
  ! First we will test the writing and reading of a matrix to a file

  integer :: n_eigenvectors = 10

  integer :: myid
  real*8, allocatable :: eigenvals(:)
  integer :: matrixsize, block_rows, block_cols


  !  Pharse command line argumnents, if given
   INTEGER*4 :: iargc
   character*16 arg1

   if (iargc() == 2) then
      call getarg(1, arg1)
      read(arg1, *) n_eigenvectors
   endif


  ! Now set some ELSI specifications
  call elsi_initialize_mpi()
  call elsi_initialize_problem_from_file("elsi_eigenvalue_problem.hdf5")
  call elsi_initialize_blacs()
  call elsi_set_method(ELPA)
  call elsi_set_mode(REAL_VALUES)

  ! Read eigenvalue problem
  call elsi_allocate_matrices()
  call elsi_read_ev_problem("elsi_eigenvalue_problem.hdf5")

  call elsi_get_global_dimensions(matrixsize,block_rows,block_cols)
  allocate (eigenvals(matrixsize))

  ! Solve the eigenvalue problem
  call elsi_solve_ev_problem(n_eigenvectors)
  
  call elsi_get_eigenvalues(eigenvals, matrixsize)

  if (myid == 0) print *, "Eigenvalues : ", eigenvals

  deallocate(eigenvals)

  ! elsi shutdown
  call elsi_finalize()

end program
