!> ELSI Interchange
!! 
!! Copyright of the original code rests with the authors inside the ELSI
!! consortium. The copyright of any additional modifications shall rest
!! with their original authors, but shall adhere to the licensing terms
!! distributed along with the original code in the file "COPYING".

program readwrite

  use iso_c_binding
  use ELSI

  implicit none

  ! This is the ELSI test suite
  ! First we will test the writing and reading of a matrix to a file

  ! Ok, let us first create a matrix for the ELSI interface
  ! Something like

  integer, parameter :: matrixsize = 100, blocksize = 16
 
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
  integer :: error
  integer :: i_row, i_col, p_row, p_col, l_row, l_col, g_row, g_col
  logical :: my_task

  if(myid==0) print *, 'Initialize ELSI MPI'
  call elsi_initialize_mpi()
  if(myid==0) print *, 'DONE'

  if(myid==0) print *, 'Initialize ELSI BLACS'
  call elsi_initialize_blacs(matrixsize, blocksize)
  if(myid==0) print *, 'DONE'

  do i_row = 1, matrixsize
    my_task = elsi_get_local_row(i_row, p_row, l_row)
    call elsi_get_global_row(g_row, p_row, l_row)
    do i_col = 1, matrixsize
      my_task = elsi_get_local_col(i_col, p_col, l_col)
      call elsi_get_global_col(g_col, p_col, l_col)
      if(myid==0) write(*,'(8(a,I3),a)') &
        '(',i_row,',',i_col,') mapped to process (', &
            p_row,',',p_col,') and local (',l_row,',',l_col,') &
        mapped back to (', g_row,',',g_col,')'
     end do
  end do

  call MPI_BARRIER(mpi_comm_world,error)

  if(myid==0) print *, 'Initialize H'
  ! Generate some data
  allocate(h(matrixsize,matrixsize))

  iseed(:) = myid
  call RANDOM_SEED(put=iseed)
  !call RANDOM_NUMBER(buffer)
  !h = buffer

  ! Processor matrix
  do i_row = 1, matrixsize
    if (elsi_get_local_row(i_row, p_row, l_row)) then
      do i_col = 1, matrixsize
        if (elsi_get_local_col(i_col, p_col, l_col)) then
          h(l_row,l_col) = 1.0 * myid
        end if
      end do
    end if
  end do

  ! h = h + buffer**T
  !call pdtran( matrixsize, matrixsize, 1.d0, buffer, 1, 1, sc_desc, &
  !             1.d0, h, 1, 1, sc_desc)

  if(myid==0) print *, 'DONE'


  if(myid==0) print *, 'Initialize S'
  allocate(s(matrixsize,matrixsize))
  allocate(buffer(matrixsize,matrixsize))
  call RANDOM_NUMBER(buffer)
  s = buffer
  ! s = s + buffer**T
  call pdtran( matrixsize, matrixsize, 1.d0, buffer, 1, 1, sc_desc, &
               1.d0, s, 1, 1, sc_desc)
  if(myid==0) print *, 'DONE'

  if(myid==0) print *, 'Set ELSI Modes'
  ! Now set some ELSI specifications
  call elsi_set_method(ELPA)
  call elsi_set_mode(REAL_VALUES)
  if(myid==0) print *, 'DONE'

  if(myid==0) print *, 'ELSI allocate Matrices'
  ! Initialize the data space
  call elsi_allocate_matrices(n_rows,n_cols)
  if(myid==0) print *, 'DONE'

  if(myid==0) print *, 'ELSI set Matrices'
  ! Set matrices
  call elsi_set_hamiltonian(h, n_rows, n_cols)
  call elsi_set_overlap(s, n_rows, n_cols)
  if(myid==0) print *, 'DONE'

  if(myid==0) print *, 'ELSI write Matrices'
  ! Write eigenvalue problem
  call elsi_write_ev_problem("elsi_eigenvalue_problem.hdf5")
  if(myid==0) print *, 'DONE'

  if(myid==0) print *, 'ELSI read Matrices'
  ! Read eigenvalue problem 
  call elsi_read_ev_problem("elsi_eigenvalue_problem.hdf5")
  if(myid==0) print *, 'DONE'

  call elsi_get_hamiltonian(buffer, n_rows, n_cols)
  buffer = buffer - h

  errmax = 0
  do i=1,n_cols
      colerr = 0
      call pdnrm2( matrixsize, colerr, buffer, 1, i, sc_desc, 1)
      if (colerr > errmax) errmax = colerr
  enddo

  ! Get maximum error norm over all processors
  colerr = errmax
  call mpi_allreduce(colerr,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
  if(myid==0) print *
  if(myid==0) print *,'Error Residual for H    :',errmax

  call elsi_get_real_overlap(buffer, n_rows, n_cols)
  buffer = buffer - s

  errmax = 0
  do i=1,n_cols
      colerr = 0
      call pdnrm2( matrixsize, colerr, buffer, 1, i, sc_desc, 1)
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

end program
