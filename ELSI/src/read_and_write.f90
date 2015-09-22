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

  ! This is the ELSI test suite
  ! First we will test the writing and reading of a matrix to a file

  ! Local variables
  integer :: myid

  ! Now set some ELSI specifications
  call elsi_initialize_mpi(myid)
  call elsi_initialize_problem_from_file("elsi_eigenvalue_problem.hdf5")
  call elsi_set_method(ELPA)
  call elsi_set_mode(REAL_VALUES)

  ! Read eigenvalue problem
  call elsi_allocate_matrices()
  call elsi_read_ev_problem("elsi_eigenvalue_problem.hdf5")

  ! Write eigenvalue problem to another file
  call elsi_write_ev_problem("elsi_eigenvalue_problem_out.hdf5")

  ! elsi shutdown
  call elsi_finalize()

end program
