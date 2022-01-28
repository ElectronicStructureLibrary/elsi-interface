











module global_product
  use precision
  implicit none
  private

  public :: global_product_double
  public :: global_product_single

  contains

! real double precision first
















subroutine global_product_&
&double&
&(obj, z, n, mpi_comm_rows, mpi_comm_cols, npc_0, npc_n)
  ! This routine calculates the global product of z.
  use precision
  use elpa_abstract_impl
  use elpa_mpi
  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)               :: mpi_comm_cols, mpi_comm_rows
  integer(kind=ik), intent(in)               :: npc_0, npc_n
  integer(kind=MPI_KIND)                     :: mpierr,my_pcolMPI, np_colsMPI, np_rowsMPI
  integer(kind=ik)                           :: n, my_pcol, np_cols, np_rows
  real(kind=rk8)                   :: z(n)

  real(kind=rk8)                   :: tmp(n)
  integer(kind=ik)                           :: np
  integer(kind=MPI_KIND)                     :: allreduce_request1, allreduce_request2
  logical                                    :: useNonBlockingCollectivesCols
  logical                                    :: useNonBlockingCollectivesRows
  integer(kind=c_int)                        :: non_blocking_collectives_rows, error, &
                                                non_blocking_collectives_cols

 call obj%get("nbc_row_global_product", non_blocking_collectives_rows, error)
 if (error .ne. ELPA_OK) then
   print *,"Problem setting option for non blocking collectives for rows in global_product. Aborting..."
   stop
 endif

 call obj%get("nbc_col_global_product", non_blocking_collectives_cols, error)
 if (error .ne. ELPA_OK) then
   print *,"Problem setting option for non blocking collectives for cols in global_product. Aborting..."
   stop
 endif

 if (non_blocking_collectives_rows .eq. 1) then
   useNonBlockingCollectivesRows = .true.
 else
   useNonBlockingCollectivesRows = .false.
 endif

 if (non_blocking_collectives_cols .eq. 1) then
   useNonBlockingCollectivesCols = .true.
 else
   useNonBlockingCollectivesCols = .false.
 endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
  np_rows = int(np_rowsMPI,kind=c_int)

  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
  my_pcol = int(my_pcolMPI,kind=c_int)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)
  np_cols = int(np_colsMPI,kind=c_int)

  !!my_pcol = int(my_pcolMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  if (npc_n==1 .and. np_rows==1) return ! nothing to do

  ! Do an mpi_allreduce over processor rows

  if (useNonBlockingCollectivesRows) then
    call obj%timer%start("mpi_nbc_communication")
    call mpi_iallreduce(z, tmp, int(n,kind=MPI_KIND), MPI_REAL8, MPI_PROD, int(mpi_comm_rows,kind=MPI_KIND), &
           allreduce_request1, mpierr)
    call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
    call obj%timer%stop("mpi_nbc_communication")
  else
    call obj%timer%start("mpi_communication")
    call mpi_allreduce(z, tmp, int(n,kind=MPI_KIND), MPI_REAL8, MPI_PROD, int(mpi_comm_rows,kind=MPI_KIND), &
           mpierr)
    call obj%timer%stop("mpi_communication")
  endif

  ! If only 1 processor column, we are done
  if (npc_n==1) then
    z(:) = tmp(:)
    return
  endif

  ! If all processor columns are involved, we can use mpi_allreduce
  if (npc_n==np_cols) then
    if (useNonBlockingCollectivesCols) then
      call obj%timer%start("mpi_nbc_communication")
      call mpi_iallreduce(tmp, z, int(n,kind=MPI_KIND), MPI_REAL8, MPI_PROD, int(mpi_comm_cols,kind=MPI_KIND), &
            allreduce_request2, mpierr)
      call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
      call obj%timer%stop("mpi_nbc_communication")
    else
      call obj%timer%start("mpi_communication")
      call mpi_allreduce(tmp, z, int(n,kind=MPI_KIND), MPI_REAL8, MPI_PROD, int(mpi_comm_cols,kind=MPI_KIND), &
             mpierr)
      call obj%timer%stop("mpi_communication")
    endif
    return
  endif

  ! We send all vectors to the first proc, do the product there
  ! and redistribute the result.

  if (my_pcol == npc_0) then
    z(1:n) = tmp(1:n)
    do np = npc_0+1, npc_0+npc_n-1
      call obj%timer%start("mpi_communication")
      call mpi_recv(tmp, int(n,kind=MPI_KIND), MPI_REAL8, int(np,kind=MPI_KIND), 1117_MPI_KIND, &
                    int(mpi_comm_cols,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
      call obj%timer%stop("mpi_communication")
      z(1:n) = z(1:n)*tmp(1:n)
    enddo
    do np = npc_0+1, npc_0+npc_n-1
      call obj%timer%start("mpi_communication")
      call mpi_send(z, int(n,kind=MPI_KIND), MPI_REAL8, int(np,kind=MPI_KIND), 1118_MPI_KIND, &
                    int(mpi_comm_cols,kind=MPI_KIND), mpierr)
      call obj%timer%stop("mpi_communication")
    enddo
  else
    call obj%timer%start("mpi_communication")
    call mpi_send(tmp, int(n,kind=MPI_KIND), MPI_REAL8, int(npc_0,kind=MPI_KIND), 1117_MPI_KIND, &
                  int(mpi_comm_cols,kind=MPI_KIND), mpierr)
    call mpi_recv(z, int(n,kind=MPI_KIND), MPI_REAL8, int(npc_0,kind=MPI_KIND), 1118_MPI_KIND, &
                  int(mpi_comm_cols,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
    call obj%timer%stop("mpi_communication")

  endif
end subroutine global_product_&
        &double

! real single precision first


















subroutine global_product_&
&single&
&(obj, z, n, mpi_comm_rows, mpi_comm_cols, npc_0, npc_n)
  ! This routine calculates the global product of z.
  use precision
  use elpa_abstract_impl
  use elpa_mpi
  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)               :: mpi_comm_cols, mpi_comm_rows
  integer(kind=ik), intent(in)               :: npc_0, npc_n
  integer(kind=MPI_KIND)                     :: mpierr,my_pcolMPI, np_colsMPI, np_rowsMPI
  integer(kind=ik)                           :: n, my_pcol, np_cols, np_rows
  real(kind=rk4)                   :: z(n)

  real(kind=rk4)                   :: tmp(n)
  integer(kind=ik)                           :: np
  integer(kind=MPI_KIND)                     :: allreduce_request1, allreduce_request2
  logical                                    :: useNonBlockingCollectivesCols
  logical                                    :: useNonBlockingCollectivesRows
  integer(kind=c_int)                        :: non_blocking_collectives_rows, error, &
                                                non_blocking_collectives_cols

 call obj%get("nbc_row_global_product", non_blocking_collectives_rows, error)
 if (error .ne. ELPA_OK) then
   print *,"Problem setting option for non blocking collectives for rows in global_product. Aborting..."
   stop
 endif

 call obj%get("nbc_col_global_product", non_blocking_collectives_cols, error)
 if (error .ne. ELPA_OK) then
   print *,"Problem setting option for non blocking collectives for cols in global_product. Aborting..."
   stop
 endif

 if (non_blocking_collectives_rows .eq. 1) then
   useNonBlockingCollectivesRows = .true.
 else
   useNonBlockingCollectivesRows = .false.
 endif

 if (non_blocking_collectives_cols .eq. 1) then
   useNonBlockingCollectivesCols = .true.
 else
   useNonBlockingCollectivesCols = .false.
 endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
  np_rows = int(np_rowsMPI,kind=c_int)

  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
  my_pcol = int(my_pcolMPI,kind=c_int)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)
  np_cols = int(np_colsMPI,kind=c_int)

  !!my_pcol = int(my_pcolMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  if (npc_n==1 .and. np_rows==1) return ! nothing to do

  ! Do an mpi_allreduce over processor rows

  if (useNonBlockingCollectivesRows) then
    call obj%timer%start("mpi_nbc_communication")
    call mpi_iallreduce(z, tmp, int(n,kind=MPI_KIND), MPI_REAL4, MPI_PROD, int(mpi_comm_rows,kind=MPI_KIND), &
           allreduce_request1, mpierr)
    call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
    call obj%timer%stop("mpi_nbc_communication")
  else
    call obj%timer%start("mpi_communication")
    call mpi_allreduce(z, tmp, int(n,kind=MPI_KIND), MPI_REAL4, MPI_PROD, int(mpi_comm_rows,kind=MPI_KIND), &
           mpierr)
    call obj%timer%stop("mpi_communication")
  endif

  ! If only 1 processor column, we are done
  if (npc_n==1) then
    z(:) = tmp(:)
    return
  endif

  ! If all processor columns are involved, we can use mpi_allreduce
  if (npc_n==np_cols) then
    if (useNonBlockingCollectivesCols) then
      call obj%timer%start("mpi_nbc_communication")
      call mpi_iallreduce(tmp, z, int(n,kind=MPI_KIND), MPI_REAL4, MPI_PROD, int(mpi_comm_cols,kind=MPI_KIND), &
            allreduce_request2, mpierr)
      call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
      call obj%timer%stop("mpi_nbc_communication")
    else
      call obj%timer%start("mpi_communication")
      call mpi_allreduce(tmp, z, int(n,kind=MPI_KIND), MPI_REAL4, MPI_PROD, int(mpi_comm_cols,kind=MPI_KIND), &
             mpierr)
      call obj%timer%stop("mpi_communication")
    endif
    return
  endif

  ! We send all vectors to the first proc, do the product there
  ! and redistribute the result.

  if (my_pcol == npc_0) then
    z(1:n) = tmp(1:n)
    do np = npc_0+1, npc_0+npc_n-1
      call obj%timer%start("mpi_communication")
      call mpi_recv(tmp, int(n,kind=MPI_KIND), MPI_REAL4, int(np,kind=MPI_KIND), 1117_MPI_KIND, &
                    int(mpi_comm_cols,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
      call obj%timer%stop("mpi_communication")
      z(1:n) = z(1:n)*tmp(1:n)
    enddo
    do np = npc_0+1, npc_0+npc_n-1
      call obj%timer%start("mpi_communication")
      call mpi_send(z, int(n,kind=MPI_KIND), MPI_REAL4, int(np,kind=MPI_KIND), 1118_MPI_KIND, &
                    int(mpi_comm_cols,kind=MPI_KIND), mpierr)
      call obj%timer%stop("mpi_communication")
    enddo
  else
    call obj%timer%start("mpi_communication")
    call mpi_send(tmp, int(n,kind=MPI_KIND), MPI_REAL4, int(npc_0,kind=MPI_KIND), 1117_MPI_KIND, &
                  int(mpi_comm_cols,kind=MPI_KIND), mpierr)
    call mpi_recv(z, int(n,kind=MPI_KIND), MPI_REAL4, int(npc_0,kind=MPI_KIND), 1118_MPI_KIND, &
                  int(mpi_comm_cols,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
    call obj%timer%stop("mpi_communication")

  endif
end subroutine global_product_&
        &single

end module
