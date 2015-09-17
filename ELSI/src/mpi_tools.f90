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
  integer :: n_dim
  integer :: n_block
  integer, external :: numroc

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

subroutine elsi_initialize_blacs(n_dim_in,n_block_in)

   implicit none

   integer, intent(in) :: n_dim_in   !< dimension of matrix
   integer, intent(in) :: n_block_in !< dimension of sub blocks

   n_dim   = n_dim_in
   n_block = n_block_in

  ! Define blockcyclic setup
  do np_cols = NINT(SQRT(REAL(n_procs))),2,-1
     if(mod(n_procs,np_cols) == 0 ) exit
  enddo

  np_rows = n_procs / np_cols


  ! Set up BLACS and MPI communicators

  blacs_ctxt = mpi_comm_world
  call BLACS_Gridinit( blacs_ctxt, 'C', np_rows, np_cols )
  call BLACS_Gridinfo( blacs_ctxt, np_rows, np_cols, ip_row, ip_col )

  call mpi_comm_split(mpi_comm_world,ip_col,ip_row,mpi_comm_rows,mpierr)
  call mpi_comm_split(mpi_comm_world,ip_row,ip_col,mpi_comm_cols,mpierr)

  if (mpierr) then
      write(*,'(a)') "ELPA: Failed to get ELPA row / column communicators."
      stop
   end if


  n_rows = numroc(n_dim, n_block, ip_row, 0, np_rows)
  n_cols = numroc(n_dim, n_block, ip_col, 0, np_cols)

  if(myid == 0) print *, 'Global matrixsize: ',n_dim, ' x ', n_dim
  if(myid == 0) print *, 'Blocksize: ',n_block, ' x ', n_block
  if(myid == 0) print *, 'Processor grid: ',np_rows, ' x ', np_cols
  if(myid == 0) print *, 'Local Matrixsize: ',n_rows, ' x ', n_cols

  call MPI_BARRIER(mpi_comm_world,mpierr)

  call descinit( sc_desc, n_dim, n_dim, n_block, n_block, 0, 0, &
                 blacs_ctxt, MAX(1,n_rows), info )


end subroutine

logical function elsi_get_local_row (g_row, p_row, l_row)

   implicit none

   integer, intent(in) :: g_row   !< global row
   integer, intent(out) :: p_row   !< process row
   integer, intent(out) :: l_row   !< local row

   integer :: b_row !< local block row
   integer :: i_row !< local row index

   elsi_get_local_row = .False.

   p_row = FLOOR ( 1d0 * MOD(g_row-1,n_block*np_rows) / n_block)
   b_row = FLOOR ( 1d0 * (g_row - 1) / n_block / np_rows)
   i_row = MOD( (g_row - 1), n_block) + 1

   l_row = b_row * n_block + i_row
   
   if (p_row == ip_row) elsi_get_local_row = .True.

end function

subroutine elsi_get_global_row (g_row, p_row, l_row)

   implicit none

   integer, intent(out) :: g_row   !< global row
   integer, intent(in)  :: p_row   !< process row
   integer, intent(in)  :: l_row   !< local row

   integer :: b_row !< local block row
   integer :: i_row !< local row index


   b_row = FLOOR( 1d0 * (l_row - 1) / n_block) 
   i_row = l_row - b_row * n_block

   g_row = p_row * n_block + b_row * np_rows * n_block + i_row
  
end subroutine

logical function elsi_get_local_col (g_col, p_col, l_col)

   implicit none

   integer, intent(in) :: g_col   !< global col
   integer, intent(out) :: p_col   !< process col
   integer, intent(out) :: l_col   !< local col

   integer :: b_col !< local block col
   integer :: i_col !< local col index

   elsi_get_local_col = .False.

   p_col = FLOOR ( 1d0 * MOD(g_col-1,n_block*np_cols) / n_block)
   b_col = FLOOR ( 1d0 * (g_col - 1) / n_block / np_cols)
   i_col = MOD( (g_col - 1), n_block) + 1

   l_col = b_col * n_block + i_col

   if (p_col == ip_col) elsi_get_local_col = .True.

end function

subroutine elsi_get_global_col (g_col, p_col, l_col)

   implicit none

   integer, intent(out) :: g_col   !< global col
   integer, intent(in)  :: p_col   !< process col
   integer, intent(in)  :: l_col   !< local col

   integer :: b_col !< local block col
   integer :: i_col !< local col index


   b_col = FLOOR( 1d0 * (l_col - 1) / n_block) 
   i_col = l_col - b_col * n_block

   g_col = p_col * n_block + b_col * np_cols * n_block + i_col
  
end subroutine

end module MPI_TOOLS
