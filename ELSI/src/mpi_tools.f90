!> ELSI Interchange
!! 
!! Copyright of the original code rests with the authors inside the ELSI
!! consortium. The copyright of any additional modifications shall rest
!! with their original authors, but shall adhere to the licensing terms
!! distributed along with the original code in the file "COPYING".

module MPI_TOOLS

  use iso_c_binding
  use DIMENSIONS

  implicit none
  private

  integer, external :: numroc

  public :: elsi_initialize_mpi
  public :: elsi_initialize_blacs 
  public :: elsi_get_local_row
  public :: elsi_get_local_col
  public :: elsi_get_global_row
  public :: elsi_get_global_col
  public :: elsi_get_global_dimensions
  public :: elsi_get_local_dimensions
  public :: elsi_get_processor_grid
  public :: elsi_get_comm_grid
  public :: elsi_finalize_mpi
  public :: elsi_finalize_blacs

  contains

subroutine elsi_initialize_mpi(myid_out)

   implicit none
   include 'mpif.h'

   integer, intent(out) :: myid_out

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

   myid_out = myid

end subroutine

subroutine elsi_finalize_mpi()

   implicit none
   include 'mpif.h'

   ! MPI Finalization

   call mpi_finalize(mpierr)
   if (mpierr) then
      write(*,'(a)') "MPI: Failed to finalize MPI."
      stop
   end if

end subroutine


subroutine elsi_initialize_blacs(n_dim_in,n_block_rows,n_block_cols)

   implicit none

   include "mpif.h"

   integer, intent(in) :: n_dim_in   !< dimension of matrix
   integer, intent(in) :: n_block_rows !< dimension of sub blocks
   integer, intent(in) :: n_block_cols !< dimension of sub blocks

   n_g_rank = n_dim_in
   n_b_rows = n_block_rows
   n_b_cols = n_block_cols

  ! Define blockcyclic setup
  do n_p_cols = NINT(SQRT(REAL(n_procs))),2,-1
     if(mod(n_procs,n_p_cols) == 0 ) exit
  enddo

  n_p_rows = n_procs / n_p_cols


  ! Set up BLACS and MPI communicators

  blacs_ctxt = mpi_comm_world
  call BLACS_Gridinit( blacs_ctxt, 'C', n_p_rows, n_p_cols )
  call BLACS_Gridinfo( blacs_ctxt, n_p_rows, n_p_cols, my_p_row, my_p_col )

  call mpi_comm_split(mpi_comm_world,my_p_col,my_p_row,mpi_comm_row,mpierr)

   if (mpierr) then
      write(*,'(a)') "ELPA: Failed to get ELPA row communicators."
      stop
   end if
 
  call mpi_comm_split(mpi_comm_world,my_p_row,my_p_col,mpi_comm_col,mpierr)

  if (mpierr) then
      write(*,'(a)') "ELPA: Failed to get ELPA column communicators."
      stop
   end if

  n_l_rows = numroc(n_g_rank, n_b_rows, my_p_row, 0, n_p_rows)
  n_l_cols = numroc(n_g_rank, n_b_cols, my_p_col, 0, n_p_cols)

  if(myid == 0) print *, 'Global matrixsize: ',n_g_rank, ' x ', n_g_rank
  if(myid == 0) print *, 'Blocksize: ',n_b_rows, ' x ', n_b_cols
  if(myid == 0) print *, 'Processor grid: ',n_p_rows, ' x ', n_p_cols
  if(myid == 0) print *, 'Local Matrixsize: ',n_l_rows, ' x ', n_l_cols
  call MPI_BARRIER(mpi_comm_world,mpierr)

  print *, myid," Id MPI_COMM_WORLD: ", mpi_comm_world
  print *, myid," Id ip_row/col: ", my_p_row, "/", my_p_col
  print *, myid," Id MPI_COMM_rows/cols: ", mpi_comm_row, "/", mpi_comm_col
  call MPI_BARRIER(mpi_comm_world,mpierr)

  call descinit( sc_desc, n_g_rank, n_g_rank, n_b_rows, n_b_cols, 0, 0, &
                 blacs_ctxt, MAX(1,n_l_rows), blacs_info )


end subroutine

subroutine elsi_finalize_blacs()

   implicit none

   call blacs_gridexit(blacs_ctxt)

end subroutine


logical function elsi_get_local_row (g_row, p_row, l_row)

   implicit none

   integer, intent(in) :: g_row   !< global row
   integer, intent(out) :: p_row   !< process row
   integer, intent(out) :: l_row   !< local row

   integer :: b_row !< local block row
   integer :: i_row !< local row index

   elsi_get_local_row = .False.

   p_row = FLOOR ( 1d0 * MOD(g_row-1,n_b_rows * n_p_rows) / n_b_rows)
   b_row = FLOOR ( 1d0 * (g_row - 1) / n_b_rows / n_p_rows)
   i_row = MOD( (g_row - 1), n_b_rows) + 1

   l_row = b_row * n_b_rows + i_row
   
   if (p_row == my_p_row) elsi_get_local_row = .True.

end function

subroutine elsi_get_global_row (g_row, p_row, l_row)

   implicit none

   integer, intent(out) :: g_row   !< global row
   integer, intent(in)  :: p_row   !< process row
   integer, intent(in)  :: l_row   !< local row

   integer :: b_row !< local block row
   integer :: i_row !< local row index


   b_row = FLOOR( 1d0 * (l_row - 1) / n_b_rows) 
   i_row = l_row - b_row * n_b_rows

   g_row = p_row * n_b_rows + b_row * n_p_rows * n_b_rows + i_row
  
end subroutine

logical function elsi_get_local_col (g_col, p_col, l_col)

   implicit none

   integer, intent(in) :: g_col   !< global col
   integer, intent(out) :: p_col   !< process col
   integer, intent(out) :: l_col   !< local col

   integer :: b_col !< local block col
   integer :: i_col !< local col index

   elsi_get_local_col = .False.

   p_col = FLOOR ( 1d0 * MOD(g_col-1,n_b_cols * n_p_cols) / n_b_cols)
   b_col = FLOOR ( 1d0 * (g_col - 1) / n_b_cols / n_p_cols)
   i_col = MOD( (g_col - 1), n_b_cols) + 1

   l_col = b_col * n_b_cols + i_col

   if (p_col == my_p_col) elsi_get_local_col = .True.

end function

subroutine elsi_get_global_col (g_col, p_col, l_col)

   implicit none

   integer, intent(out) :: g_col   !< global col
   integer, intent(in)  :: p_col   !< process col
   integer, intent(in)  :: l_col   !< local col

   integer :: b_col !< local block col
   integer :: i_col !< local col index


   b_col = FLOOR( 1d0 * (l_col - 1) / n_b_cols) 
   i_col = l_col - b_col * n_b_cols

   g_col = p_col * n_b_cols + b_col * n_p_cols * n_b_cols + i_col
  
end subroutine

subroutine elsi_get_global_dimensions (g_dim, g_block_rows, g_block_cols)

   implicit none

   integer, intent(out) :: g_dim   !< global rank
   integer, intent(out) :: g_block_rows !< global blocksize
   integer, intent(out) :: g_block_cols !< global blocksize

   g_dim   = n_g_rank 
   g_block_rows = n_b_rows
   g_block_cols = n_b_cols

end subroutine

subroutine elsi_get_local_dimensions (rows, cols)

   implicit none

   integer, intent(out) :: rows   !< output local rows
   integer, intent(out) :: cols   !< output local cols

   rows  = n_l_rows 
   cols  = n_l_cols

end subroutine

subroutine elsi_get_processor_grid (rows, cols)

   implicit none

   integer, intent(out) :: rows   !< output processor rows
   integer, intent(out) :: cols   !< output processor cols

   rows  = n_p_rows 
   cols  = n_p_cols

end subroutine

subroutine elsi_get_comm_grid (rows, cols)

   implicit none

   integer, intent(out) :: rows   !< output row communicator
   integer, intent(out) :: cols   !< output col communicator

   rows  = mpi_comm_row 
   cols  = mpi_comm_col

end subroutine


end module MPI_TOOLS
