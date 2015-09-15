!> ELSI Interchange
!! 
!! Copyright of the original code rests with the authors inside the ELSI
!! consortium. The copyright of any additional modifications shall rest
!! with their original authors, but shall adhere to the licensing terms
!! distributed along with the original code in the file "COPYING".

module MPI_TOOLS

  use iso_c_binding

  implicit none
  include 'mpif.h'

  ! Coarse grid
  integer :: np_rows        !< Number of rows (coarse grid)
  integer :: np_cols        !< Number of cols (coarse grid)
  
  ! Sub grid
  integer :: n_rows        !< Total number of rows for this process
  integer :: n_cols        !< Total number of columns for this process

  ! MPI variables
  integer :: myid           !< local process id
  integer :: mpierr         !< mpi error handler
  integer :: n_procs        !< number of mpi processes

  ! BLACS variable
  integer :: blacs_ctxt     !< local blacs context
  integer :: sc_desc(9)     !< local blacs context
  integer :: mpi_comm_rows  !< row communicatior
  integer :: mpi_comm_cols  !< row communicatior
  integer :: ip_row         !< local row    position
  integer :: ip_col         !< local column position
  integer :: info

  contains

subroutine elsi_initialize_mpi()

   implicit none

   ! MPI Initialization

   call mpi_init(mpierr)
   if (mpierr) then
      write(*,'(a)') "MPI: Failed to initialize MPI."
      stop
   end if
   call mpi_comm_rank(mpi_comm_world, myid, mpierr)
   if (mpierr) then
      write(*,'(a)') "MPI: Failed to determine MPI rank."
      stop
   end if

   call mpi_comm_size(mpi_comm_world, n_procs, mpierr)
   if (mpierr) then
      write(*,'(a)') "MPI: Failed to determine MPI size."
      stop
   end if


end subroutine

subroutine elsi_initialize_blacs(n_dim,n_block)

   implicit none

   integer, intent(in) :: n_dim !< dimension of matrix
   integer, intent(in) :: n_block !< dimension of sub blocks

  ! Define blockcyclic setup
  do np_cols = NINT(SQRT(REAL(n_procs))),2,-1
     if(mod(n_procs,np_cols) == 0 ) exit
  enddo

  np_rows = n_procs / np_cols


  ! Set up BLACS and MPI communicators

  blacs_ctxt = mpi_comm_world
  call BLACS_Gridinit( blacs_ctxt, 'C', np_rows, np_cols )
  call BLACS_Gridinfo( blacs_ctxt, np_rows, np_cols, ip_row, ip_col )

  call descinit( sc_desc, n_dim, n_dim, n_block, n_block, 0, 0, &
                 blacs_ctxt, n_rows, info )


end subroutine


end module MPI_TOOLS
