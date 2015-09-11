!> ELSI Interchange
!! 
!! Copyright of the original code rests with the authors inside the ELSI
!! consortium. The copyright of any additional modifications shall rest
!! with their original authors, but shall adhere to the licensing terms
!! distributed along with the original code in the file "COPYING".

subroutine readwrite

  use iso_c_binding
  use ELSI

  implicit none
  include 'mpif.h'

  ! This is the ELSI test suite
  ! First we will test the writing and reading of a matrix to a file

  ! Ok, let us first create a matrix for the ELSI interface
  ! Something like

  integer, parameter :: n_dim = 100

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

  ! Random number generator
  integer :: iseed(4096) !< Random seed

  ! Test matrices
  real*8, allocatable :: h(:,:)
  real*8, allocatable :: s(:,:)
  real*8, allocatable :: buffer(:,:) 

  ! Local variables
  real*8 :: errmax
  real*8 :: colerr
  integer :: i

  ! MPI Initialization

  call mpi_init(mpierr)
  call mpi_comm_rank(mpi_comm_world, myid, mpierr)
  call mpi_comm_size(mpi_comm_world, n_procs, mpierr)
 
  ! Define blockcyclic setup
  do np_cols = NINT(SQRT(REAL(n_procs))),2,-1
     if(mod(n_procs,np_cols) == 0 ) exit
  enddo

  np_rows = n_procs / np_cols

  ! Set up BLACS and MPI communicators

  blacs_ctxt = mpi_comm_world
  call BLACS_Gridinit( blacs_ctxt, 'C', np_rows, np_cols )
  call BLACS_Gridinfo( blacs_ctxt, np_rows, np_cols, ip_row, ip_col )

  call descinit( sc_desc, n_dim, n_dim, 16, 16, 0, 0, &
                 blacs_ctxt, n_rows, info )

  call get_elpa_row_col_comms( mpi_comm_world, ip_row, ip_col, &
                               mpi_comm_rows, mpi_comm_cols)

  ! Generate some data
  allocate(h(n_dim,n_dim))
  allocate(s(n_dim,n_dim))

  iseed(:) = myid
  call RANDOM_SEED(put=iseed)
  call RANDOM_NUMBER(buffer)

  ! Symmetrize hamiltonian
  h = buffer
  ! h = h + buffer**T
  call pdtran( n_dim, n_dim, 1.d0, buffer, 1, 1, sc_desc, &
               1.d0, h, 1, 1, sc_desc)

  call RANDOM_NUMBER(buffer)
  s = buffer
  ! s = s + buffer**T
  call pdtran( n_dim, n_dim, 1.d0, buffer, 1, 1, sc_desc, &
               1.d0, s, 1, 1, sc_desc)


  ! Now set some ELSI specifications
  call elsi_set_method(ELPA)
  call elsi_set_mode(REAL_VALUES)

  ! Initialize the data space
  call elsi_allocate_matrices(n_rows,n_cols)

  ! Set matrices
  call elsi_set_hamiltonian(h)
  call elsi_set_overlap(s)

  ! Write eigenvalue problem
  call elsi_write_ev_problem("elsi_eigenvalue_problem.hdf5")

  ! Read eigenvalue problem 
  call elsi_read_ev_problem("elsi_eigenvalue_problem.hdf5")

  call elsi_get_hamiltonian(buffer)
  buffer = buffer - h

  errmax = 0
  do i=1,n_cols
      colerr = 0
      call pdnrm2( n_dim, colerr, buffer, 1, i, sc_desc, 1)
      if (colerr > errmax) errmax = colerr
  enddo

  ! Get maximum error norm over all processors
  colerr = errmax
  call mpi_allreduce(colerr,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
  if(myid==0) print *
  if(myid==0) print *,'Error Residual for H    :',errmax

  call elsi_get_overlap(buffer)
  buffer = buffer - s

  errmax = 0
  do i=1,n_cols
      colerr = 0
      call pdnrm2( n_dim, colerr, buffer, 1, i, sc_desc, 1)
      if (colerr > errmax) errmax = colerr
  enddo

  ! Get maximum error norm over all processors
  colerr = errmax
  call mpi_allreduce(colerr,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
  if(myid==0) print *
  if(myid==0) print *,'Error Residual for S    :',errmax


  ! Deallocations
  deallocate(h) 
  deallocate(s) 


  ! exit BLACS grid
  call blacs_gridexit(blacs_ctxt)

  ! finalize mpi
  call mpi_finalize(mpierr)

end subroutine readwrite
