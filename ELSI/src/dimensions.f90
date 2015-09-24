!> ELSI Interchange
!! 
!! Copyright of the original code rests with the authors inside the ELSI
!! consortium. The copyright of any additional modifications shall rest
!! with their original authors, but shall adhere to the licensing terms
!! distributed along with the original code in the file "COPYING".

module DIMENSIONS

  use iso_c_binding

  implicit none

  !> Global Matrix Dimensions
  integer :: n_g_rank       !< Rank of eigenvalue problem 
  
  !> Global Matrix Blocks
  integer :: n_b_rows       !< Number of rows 
  integer :: n_b_cols       !< Number of columns 

  !> Processor grid
  integer :: n_p_rows       !< Number of processor rows 
  integer :: n_p_cols       !< Number of processor columns
  
  ! Local Matrix Dimensions
  integer :: n_l_rows        !< Number of rows for this process
  integer :: n_l_cols        !< Number of columns for this process

  ! MPI variables
  integer :: myid            !< local process id
  integer :: mpierr          !< mpi error handler
  integer :: n_procs         !< number of mpi processes
  integer :: mpi_comm_global !< global mpi communicator

  ! BLACS variables
  integer :: blacs_ctxt     !< local blacs context
  integer :: sc_desc(9)     !< local blacs context
  integer :: mpi_comm_row   !< row    communicatior
  integer :: mpi_comm_col   !< column communicatior
  integer :: my_p_row       !< process row    position
  integer :: my_p_col       !< process column position
  integer :: blacs_info

  ! HDF5 variables
  integer :: h5err

  contains

  !EMPTY

end module DIMENSIONS
