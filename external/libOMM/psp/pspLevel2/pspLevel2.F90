MODULE pspLevel2
  use pspVariable
  use pspUtility
  use pspMPI
  use pspLevel1

#ifdef MPI
  include 'mpif.h'
#endif

  private

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**********************************************!

  integer, external :: numroc ! it is a function to compute local size

  !**** INTERFACES ********************************!


  interface psp_den2sp_m
     module procedure psp_den2sp_dm
     module procedure psp_den2sp_zm
  end interface psp_den2sp_m

  !interface psp_sp2den_m
  !  module procedure psp_sp2den_dm
  !  module procedure psp_sp2den_zm
  !end interface psp_sp2den_m

  interface die
     module procedure die
  end interface die

  !public :: psp_sp2den_m
  public :: psp_den2sp_m

contains

  subroutine psp_den2sp_dm(denMat,desc,spMat,fmt,thre)

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt ! storage format of sparse matrices, 'coo' or 'csc'
    real(dp), intent(in) ::   denMat(:,:) ! dense matrix
    integer, intent(in), target :: desc(9) ! BLACS array descriptor

    !**** OPTIONAL INPUT ***********************************!
    real(dp), optional :: thre ! non-negative threshold
    ! If thre=0, generate a sparse matrix with nonzero entries.
    ! If thre>0, generate a sparse matrix with entries with an absolute value >= thre.

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout) :: spMat ! sparse matrix

    !**** INTERNAL ********************************!

    integer :: nnz, i, j, cnt
    integer :: sz(2)
    real(dp), allocatable :: val(:)
    integer, allocatable :: idx1(:)
    integer, allocatable :: idx2(:)

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
    common /coeff/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
#endif

    sz=shape(denMat)

    if (.NOT. present(thre)) then
       thre=0.0_dp
    end if
    nnz=0
    if (thre<0.0_dp) then
       call die('psp_den2sp_dm: thre must be non-negative')
    else if (thre==0.0_dp) then
       do j=1,sz(2)
          do i=1,sz(1)
             if (denMat(i,j)/=0.0_dp) then
                nnz=nnz+1
             end if
          end do
       end do
    else
       do j=1,sz(2)
          do i=1,sz(1)
             if (abs(denMat(i,j))>=thre) then
                nnz=nnz+1
             end if
          end do
       end do
    end if

    allocate(val(nnz))
    allocate(idx1(nnz))
    allocate(idx2(nnz))
    cnt=1
    if (thre==0.0_dp) then
       do j=1,sz(2)
          do i=1,sz(1)
             if (denMat(i,j)/=0.0_dp) then
                idx1(cnt)=i
                idx2(cnt)=j
                val(cnt)=denMat(i,j)
                cnt=cnt+1
             end if
          end do
       end do
    else
       do j=1,sz(2)
          do i=1,sz(1)
             if (abs(denMat(i,j))>=thre) then
                idx1(cnt)=i
                idx2(cnt)=j
                val(cnt)=denMat(i,j)
                cnt=cnt+1
             end if
          end do
       end do
    end if

    call psp_register_spm(spMat,idx1,idx2,val,desc,'coo',sz,psp_nprow,psp_npcol)
    if (fmt.EQ.'csc') then
       call psp_coo2csc(spMat)
    endif

  end subroutine psp_den2sp_dm

  subroutine psp_den2sp_zm(denMat,desc,spMat,fmt,thre)

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt ! storage format of sparse matrices, 'coo' or 'csc'
    complex(dp), intent(in) ::   denMat(:,:) ! dense matrix
    integer, intent(in), target :: desc(9) ! BLACS array descriptor

    !**** OPTIONAL INPUT ***********************************!
    real(dp), optional :: thre ! non-negative threshold
    ! If thre=0, generate a sparse matrix with nonzero entries.
    ! If thre>0, generate a sparse matrix with entries with an absolute value >= thre.

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout) :: spMat ! sparse matrix

    !**** INTERNAL ********************************!

    integer :: sz(2)
    integer :: nnz, i, j, cnt
    complex(dp), allocatable :: val(:)
    integer, allocatable :: idx1(:)
    integer, allocatable :: idx2(:)

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
    common /coeff/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
#endif

    sz=shape(denMat)

    if (.NOT. present(thre)) then
       thre=0.0_dp
    end if
    nnz=0
    if (thre<0.0_dp) then
       call die('psp_den2sp_dm: thre must be non-negative')
    else if (thre==0.0_dp) then
       do j=1,sz(2)
          do i=1,sz(1)
             if (abs(denMat(i,j))/=0.0_dp) then
                nnz=nnz+1
             end if
          end do
       end do
    else
       do j=1,sz(2)
          do i=1,sz(1)
             if (abs(denMat(i,j))>=thre) then
                nnz=nnz+1
             end if
          end do
       end do
    end if

    allocate(val(nnz))
    allocate(idx1(nnz))
    allocate(idx2(nnz))
    cnt=1
    if (thre==0.0_dp) then
       do j=1,sz(2)
          do i=1,sz(1)
             if (abs(denMat(i,j))/=0.0_dp) then
                idx1(cnt)=i
                idx2(cnt)=j
                val(cnt)=denMat(i,j)
                cnt=cnt+1
             end if
          end do
       end do
    else
       do j=1,sz(2)
          do i=1,sz(1)
             if (abs(denMat(i,j))>=thre) then
                idx1(cnt)=i
                idx2(cnt)=j
                val(cnt)=denMat(i,j)
                cnt=cnt+1
             end if
          end do
       end do
    end if

    call psp_register_spm(spMat,idx1,idx2,val,desc,'coo',sz,psp_nprow,psp_npcol)
    if (fmt.EQ.'csc') then
       call psp_coo2csc(spMat)
    endif

  end subroutine psp_den2sp_zm


  subroutine die(message)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in), optional :: message

    !**** INTERNAL ********************************!

    logical, save :: log_start=.false.

    integer :: log_unit

    !**********************************************!

    if (log_start) then
       open(newunit=log_unit,file='MatrixSwitch.log',position='append')
    else
       open(newunit=log_unit,file='MatrixSwitch.log',status='replace')
       log_start=.true.
    end if
    write(log_unit,'(a)'), 'FATAL ERROR in matrix_switch!'
    write(log_unit,'(a)'), message
    close(log_unit)
    stop

  end subroutine die


END MODULE pspLevel2
