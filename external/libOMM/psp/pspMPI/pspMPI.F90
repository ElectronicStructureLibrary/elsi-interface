MODULE pspMPI

#ifdef MPI
  include 'mpif.h'
#endif

  private

  !**** INTERFACES ********************************!

  interface psp_gridinit
     module procedure psp_gridinit
  end interface psp_gridinit

  public :: psp_gridinit

contains

  !================================================!
  ! initialize grid topology for psp               !
  !================================================!
  subroutine psp_gridinit(mpi_size,nprow,order,bs_def_row,bs_def_col,icontxt)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: order ! ordering of processor grid: 'r/R' or other for row-major, 'c/C' for column-major
    !integer, intent(in) :: mpi_comm_world
    integer, intent(in) :: mpi_size ! total number of MPI processes for the processor grid
    integer, intent(in) :: nprow ! number of rows in the processor grid
    integer, intent(in) :: bs_def_row ! default block size in row
    integer, intent(in) :: bs_def_col ! default block size in column
    integer, intent(in), optional :: icontxt ! existing BLACS context handle in case ScaLAPACK is already initialized

    !**** INTERNAL ********************************!

    integer :: i, mpi_err
    logical :: remain_dims(2), periods(2), reorder
    integer :: mpi_rank, dims(2), coords(2)

    !**** GLOBAL **********************************!
#ifdef MPI
    character(1) :: psp_proc_order

    integer :: psp_mpi_comm_world
    integer :: psp_mpi_size
    integer :: psp_nprow
    integer :: psp_npcol
    integer :: psp_bs_def_row
    integer :: psp_bs_def_col
    integer :: psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    integer :: psp_bs_num
    integer :: psp_icontxt ! BLACS context handle used by psp
    integer :: psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col


    common /coeff/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /coeff/ psp_bs_def_row, psp_bs_def_col
    common /coeff/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /coeff/ psp_bs_num
    common /coeff/ psp_icontxt ! BLACS context handle used by psp
    common /coeff/ psp_mpi_comm_cart, psp_mpi_comm_row,  psp_mpi_comm_col
#endif

    !**********************************************!

    psp_mpi_comm_world=mpi_comm_world
    psp_mpi_size=mpi_size
    psp_nprow=nprow
    psp_npcol=mpi_size/nprow
    psp_proc_order=order
    psp_bs_def_row=bs_def_row
    psp_bs_def_col=bs_def_col
    if (bs_def_row.eq.bs_def_col) then
       psp_update_rank = bs_def_row
    else
       psp_update_rank=1
    end if 

    if (present(icontxt)) then
       psp_icontxt=icontxt
    else
       call blacs_get(-1,0,psp_icontxt)
       call blacs_gridinit(psp_icontxt,psp_proc_order,psp_nprow,psp_npcol)
    end if

    ! generate sub-topology for row and column boardcasting
    dims(1)=psp_nprow
    dims(2)=psp_npcol
    periods(1)=.true.
    periods(2)=.true.
    reorder=.false.
    call MPI_Cart_create(psp_mpi_comm_world, 2, dims, periods, reorder, psp_mpi_comm_cart,mpi_err)
    ! generate sub-topology for row and column boardcasting
    remain_dims(1) = .false.
    remain_dims(2) = .true.
    call MPI_Cart_Sub(psp_mpi_comm_cart,remain_dims,psp_mpi_comm_row,mpi_err)
    remain_dims(1) = .true.
    remain_dims(2) = .false.
    call MPI_Cart_Sub(psp_mpi_comm_cart,remain_dims,psp_mpi_comm_col,mpi_err)

  end subroutine psp_gridinit

END MODULE pspMPI
