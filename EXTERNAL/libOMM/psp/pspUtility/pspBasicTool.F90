MODULE pspBasicTool
  use pspVariable
  ! This module contains routines for indexing, generating, copying, transforming (sparse) data

#ifdef MPI
  include 'mpif.h'
#endif

  private

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** INTERFACES ********************************!

  interface psp_spm_zeros
     module procedure psp_spm_zeros
  end interface psp_spm_zeros

  interface psp_den2sp_m
     module procedure psp_den2sp_dm
     module procedure psp_den2sp_zm
  end interface psp_den2sp_m

  interface psp_sp2den_m
     module procedure psp_sp2den_dm
     module procedure psp_sp2den_zm
  end interface psp_sp2den_m

  interface psp_sst_den2sp_m
     module procedure psp_sst_den2sp_dm
     module procedure psp_sst_den2sp_zm
  end interface psp_sst_den2sp_m

  interface psp_sst_sp2den_m
     module procedure psp_sst_sp2den_dm
     module procedure psp_sst_sp2den_zm
  end interface psp_sst_sp2den_m

  interface psp_sst_coo2csc
     module procedure psp_sst_coo2csc
  end interface psp_sst_coo2csc

  interface psp_sst_csc2coo
     module procedure psp_sst_csc2coo
  end interface psp_sst_csc2coo

  interface psp_sst_fmtCnvt
     module procedure psp_sst_fmtCnvt
  end interface psp_sst_fmtCnvt

  interface psp_coo2csc
     module procedure psp_coo2csc
  end interface psp_coo2csc

  interface psp_csc2coo
     module procedure psp_csc2coo
  end interface psp_csc2coo

  interface psp_copy_spm2st
     module procedure psp_copy_dspm2st
     module procedure psp_copy_zspm2st
  end interface psp_copy_spm2st

  interface psp_copy_m
     module procedure psp_copy_dm
     module procedure psp_copy_zm
  end interface psp_copy_m

  interface psp_copy_v
     module procedure psp_copy_dv
     module procedure psp_copy_zv
  end interface psp_copy_v

  interface psp_idx_glb2loc
     module procedure psp_idx_glb2loc
  end interface psp_idx_glb2loc

  interface psp_process_opM
     module procedure psp_process_lopM
     module procedure psp_process_iopM
  end interface psp_process_opM

  interface init_random_seed
     module procedure init_random_seed
  end interface init_random_seed

  interface die
     module procedure die
  end interface die

  public :: psp_sp2den_m
  public :: psp_den2sp_m
  public :: psp_sst_den2sp_m
  public :: psp_sst_sp2den_m
  public :: psp_sst_coo2csc
  public :: psp_sst_csc2coo
  public :: psp_sst_fmtCnvt
  public :: psp_coo2csc
  public :: psp_csc2coo
  public :: psp_copy_spm2st
  public :: psp_copy_m
  public :: psp_copy_v
  public :: psp_idx_glb2loc
  public :: psp_process_opM
  public :: init_random_seed
  public :: psp_spm_zeros

contains

  subroutine psp_spm_zeros(m_name,M,N,fmt,isReal)
    ! create a zero matrix of size M by N
    implicit none

    !**** INPUT ***********************************!
    character(3), intent(in), target :: fmt ! 'coo' or 'csc'
    integer, intent(in) :: M, N ! a zero matrix of size M by N
    logical, intent(in) :: isReal

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout) :: m_name ! sparse matrix to be allocated

    !**** LOCAL ***********************************!

    integer :: iprow, ipcol, info, cnt

    !**********************************************!

    integer, external :: numroc ! it is a function to compute local size

    !**** GLOBAL **********************************!
#ifdef MPI
    character(1) :: psp_proc_order

    integer :: psp_mpi_comm_world
    integer :: psp_mpi_size
    integer :: psp_nprow
    integer :: psp_npcol
    integer :: psp_bs_def_row
    integer :: psp_bs_def_col
    integer :: psp_update_rank ! In SUMMA paper for C=A*B, C is computed by rank psp_update_rank=1 local update
    integer :: psp_bs_num
    integer :: psp_icontxt ! BLACS context handle used by psp
    integer :: psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col


    common /psp_grid2D/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /psp_grid2D/ psp_bs_def_row, psp_bs_def_col
    common /psp_grid2D/ psp_update_rank ! In SUMMA paper for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /psp_grid2D/ psp_bs_num
    common /psp_grid2D/ psp_icontxt ! BLACS context handle used by psp
    common /psp_grid2D/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
#endif

    !**********************************************!

    call blacs_gridinfo(psp_icontxt,psp_nprow,psp_npcol,iprow,ipcol)

    m_name%nnz=0
    m_name%glb_dim1=M
    m_name%glb_dim2=N
    m_name%loc_dim1=numroc(M,psp_bs_def_row,iprow,0,psp_nprow)
    m_name%loc_dim2=numroc(N,psp_bs_def_col,ipcol,0,psp_npcol)
    call descinit(m_name%desc,M,N,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,m_name%loc_dim1,info)
    if (m_name%glb_dim1==m_name%glb_dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type=fmt
    m_name%is_serial=.false.
    m_name%is_real=isReal

    if (allocated(m_name%dval)) deallocate(m_name%dval)
    allocate(m_name%dval(m_name%nnz))

    if (allocated(m_name%zval)) deallocate(m_name%zval)
    allocate(m_name%zval(m_name%nnz))

    if (allocated(m_name%row_ind)) deallocate(m_name%row_ind)
    allocate(m_name%row_ind(m_name%nnz))

    if (fmt.EQ.'coo') then
       if (allocated(m_name%col_ind)) deallocate(m_name%col_ind)
       allocate(m_name%col_ind(m_name%nnz))
    else
       if (allocated(m_name%col_ptr)) deallocate(m_name%col_ptr)
       allocate(m_name%col_ptr(m_name%loc_dim2+1))
       do cnt=1,m_name%loc_dim2+1
          m_name%col_ptr(cnt) = 1
       end do
    end if

    m_name%is_initialized=.true.

  end subroutine psp_spm_zeros

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


    common /psp_grid2D/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /psp_grid2D/ psp_bs_def_row, psp_bs_def_col
    common /psp_grid2D/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /psp_grid2D/ psp_bs_num
    common /psp_grid2D/ psp_icontxt ! BLACS context handle used by psp
    common /psp_grid2D/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
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


    common /psp_grid2D/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /psp_grid2D/ psp_bs_def_row, psp_bs_def_col
    common /psp_grid2D/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /psp_grid2D/ psp_bs_num
    common /psp_grid2D/ psp_icontxt ! BLACS context handle used by psp
    common /psp_grid2D/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
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

  subroutine psp_sp2den_dm(spMat,denMat,desc)

    !**** INPUT ***********************************!
    type(psp_matrix_spm), intent(inout) :: spMat ! sparse matrix

    !**** INOUT ***********************************!
    real(dp), allocatable, intent(inout) ::   denMat(:,:) ! dense matrix
    integer, intent(inout), target :: desc(9) ! BLACS array descriptor

    !**** INTERNAL ********************************!

    integer :: cnt, cnt2

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


    common /psp_grid2D/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /psp_grid2D/ psp_bs_def_row, psp_bs_def_col
    common /psp_grid2D/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /psp_grid2D/ psp_bs_num
    common /psp_grid2D/ psp_icontxt ! BLACS context handle used by psp
    common /psp_grid2D/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
#endif

    if (allocated(denMat)) deallocate(denMat)
    allocate(denMat(spMat%loc_dim1,spMat%loc_dim2))
    denMat=0.0_dp
    desc=spMat%desc
    select case (spMat%str_type)
    case('coo')
       do cnt=1,spMat%nnz
          denMat(spMat%row_ind(cnt),spMat%col_ind(cnt))=spMat%dval(cnt)
       end do
    case('csc')
       do cnt=1,spMat%loc_dim2
          do cnt2=spMat%col_ptr(cnt),spMat%col_ptr(cnt+1)-1
             denMat(spMat%row_ind(cnt2),cnt)=spMat%dval(cnt2)
          end do
       end do
    end select

  end subroutine psp_sp2den_dm

  subroutine psp_sp2den_zm(spMat,denMat,desc)

    !**** INPUT ***********************************!
    type(psp_matrix_spm), intent(inout) :: spMat ! sparse matrix

    !**** INOUT ***********************************!
    complex(dp), allocatable, intent(inout) ::   denMat(:,:) ! dense matrix
    integer, intent(inout), target :: desc(9) ! BLACS array descriptor

    !**** INTERNAL ********************************!

    integer :: cnt, cnt2

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


    common /psp_grid2D/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /psp_grid2D/ psp_bs_def_row, psp_bs_def_col
    common /psp_grid2D/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /psp_grid2D/ psp_bs_num
    common /psp_grid2D/ psp_icontxt ! BLACS context handle used by psp
    common /psp_grid2D/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
#endif

    if (allocated(denMat)) deallocate(denMat)
    allocate(denMat(spMat%loc_dim1,spMat%loc_dim2))
    denMat=cmplx_0
    desc=spMat%desc
    select case (spMat%str_type)
    case('coo')
       do cnt=1,spMat%nnz
          denMat(spMat%row_ind(cnt),spMat%col_ind(cnt))=spMat%zval(cnt)
       end do
    case('csc')
       do cnt=1,spMat%loc_dim2
          do cnt2=spMat%col_ptr(cnt),spMat%col_ptr(cnt+1)-1
             denMat(spMat%row_ind(cnt2),cnt)=spMat%zval(cnt2)
          end do
       end do
    end select

  end subroutine psp_sp2den_zm


  subroutine psp_sst_sp2den_dm(m,n,idx1,idx2,val,fmt,denMat)

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt ! storage format of sparse matrices, 'coo' or 'csc'
    integer, intent(in)  :: m, n
    real(dp), allocatable , intent(in) :: val(:)
    integer, allocatable , intent(in) :: idx1(:)
    integer, allocatable , intent(in) :: idx2(:)

    !**** INOUT ***********************************!
    real(dp), allocatable, intent(inout) ::   denMat(:,:) ! dense matrix

    !**** INTERNAL ********************************!

    integer :: cnt, cnt2, nnz

    !**********************************************!
    nnz=size(val)
    if (allocated(denMat)) deallocate(denMat)
    allocate(denMat(m,n))
    denMat=0.0_dp
    if (nnz>0) then
       select case (fmt)
       case('coo')
          do cnt=1,nnz
             denMat(idx1(cnt),idx2(cnt))=val(cnt)
          end do
       case('csc')
          do cnt=1,n
             do cnt2=idx2(cnt),idx2(cnt+1)-1
                denMat(idx1(cnt2),cnt)=val(cnt2)
             end do
          end do
       end select
    end if

  end subroutine psp_sst_sp2den_dm

  subroutine psp_sst_sp2den_zm(m,n,idx1,idx2,val,fmt,denMat)

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt ! storage format of sparse matrices, 'coo' or 'csc'
    integer, intent(in)  :: m, n
    complex(dp), allocatable , intent(in) :: val(:)
    integer, allocatable , intent(in) :: idx1(:)
    integer, allocatable , intent(in) :: idx2(:)

    !**** INOUT ***********************************!
    complex(dp), allocatable, intent(inout) ::   denMat(:,:) ! dense matrix

    !**** INTERNAL ********************************!

    integer :: cnt, cnt2, nnz

    !**********************************************!
    nnz=size(val)
    if (allocated(denMat)) deallocate(denMat)
    allocate(denMat(m,n))
    denMat=cmplx_0
    if (nnz>0) then
       select case (fmt)
       case('coo')
          do cnt=1,nnz
             denMat(idx1(cnt),idx2(cnt))=val(cnt)
          end do
       case('csc')
          do cnt=1,n
             do cnt2=idx2(cnt),idx2(cnt+1)-1
                denMat(idx1(cnt2),cnt)=val(cnt2)
             end do
          end do
       end select
    end if

  end subroutine psp_sst_sp2den_zm

  subroutine psp_sst_den2sp_dm(denMat,idx1,idx2,val,fmt,thre)

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt ! storage format of sparse matrices, 'coo' or 'csc'
    real(dp), intent(in) ::   denMat(:,:) ! dense matrix

    !**** OPTIONAL INPUT ***********************************!
    real(dp), optional :: thre ! non-negative threshold
    ! If thre=0, generate a sparse matrix with nonzero entries.
    ! If thre>0, generate a sparse matrix with entries with an absolute value >= thre.

    !**** INOUT ***********************************!

    real(dp), allocatable , intent(inout) :: val(:)
    integer, allocatable , intent(inout) :: idx1(:)
    integer, allocatable , intent(inout) :: idx2(:)

    !**** INTERNAL ********************************!

    integer :: nnz, i, j, cnt
    integer :: sz(2)

    !**********************************************!

    sz=shape(denMat)

    if (.NOT. present(thre)) then
       thre=0.0_dp
    end if
    nnz=0
    if (thre<0.0_dp) then
       call die('psp_sst_den2sp_dm: thre must be non-negative')
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

    if (allocated(val)) deallocate(val)
    if (allocated(idx1)) deallocate(idx1)
    if (allocated(idx2)) deallocate(idx2)
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

    if (fmt.EQ.'csc') then
       call psp_sst_coo2csc(sz(1),sz(2),nnz,idx2)
    endif

  end subroutine psp_sst_den2sp_dm

  subroutine psp_sst_den2sp_zm(denMat,idx1,idx2,val,fmt,thre)

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt ! storage format of sparse matrices, 'coo' or 'csc'
    complex(dp), intent(in) ::   denMat(:,:) ! dense matrix

    !**** OPTIONAL INPUT ***********************************!
    real(dp), optional :: thre ! non-negative threshold
    ! If thre=0, generate a sparse matrix with nonzero entries.
    ! If thre>0, generate a sparse matrix with entries with an absolute value >= thre.

    !**** INOUT ***********************************!

    complex(dp), allocatable , intent(inout) :: val(:)
    integer, allocatable , intent(inout) :: idx1(:)
    integer, allocatable , intent(inout) :: idx2(:)

    !**** INTERNAL ********************************!

    integer :: nnz, i, j, cnt
    integer :: sz(2)

    !**********************************************!

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


    if (allocated(val)) deallocate(val)
    if (allocated(idx1)) deallocate(idx1)
    if (allocated(idx2)) deallocate(idx2)
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

    if (fmt.EQ.'csc') then
       call psp_sst_coo2csc(sz(1),sz(2),nnz,idx2)
    endif

  end subroutine psp_sst_den2sp_zm

  subroutine psp_sst_coo2csc(m,n,nnz,idx)

    !**** IN ***********************************!

    integer, intent(in) ::   m, n, nnz

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: idx(:)

    !**** INTERNAL ********************************!

    integer :: j, cnt, cnt2, num
    integer, allocatable :: col_ptr(:)

    !**********************************************!

    allocate(col_ptr(n+1))
    col_ptr(1)=1
    cnt=1
    do j=1,n
       num=0
       if (cnt<=nnz) then
          cnt2=0
          do while (idx(cnt+cnt2)==j)
             cnt2=cnt2+1
             num=num+1
             if (cnt+cnt2>nnz) then
                exit
             end if
          end do
       end if
       col_ptr(j+1)=col_ptr(j)+num
       cnt=cnt+num
    end do

    deallocate(idx)
    allocate(idx(n+1))
    do j=1,n+1
       idx(j)=col_ptr(j)
    end do
    deallocate(col_ptr)

  end subroutine psp_sst_coo2csc

  subroutine psp_sst_csc2coo(m,n,nnz,col_ptr)

    !**** IN ***********************************!

    integer, intent(in) ::   m, n, nnz

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: col_ptr(:)

    !**** INTERNAL ********************************!

    integer, allocatable :: col_ind(:)
    integer :: j, k, nst, jlo, jhi, klo, khi

    !**********************************************!

    allocate(col_ind(nnz))
    nst=0
    if (col_ptr(1)==0) then
       jlo=0
       jhi=n-1
       do j=jlo,jhi
          klo=col_ptr(j+1)
          khi=col_ptr(j+2)-1
          do k=klo,khi
             nst=nst+1
             col_ind(nst)=j
          end do
       end do
    else
       jlo=1
       jhi=n
       do j=jlo,jhi
          klo=col_ptr(j)
          khi=col_ptr(j+1)-1
          do k=klo,khi
             nst=nst+1
             col_ind(nst)=j
          end do
       end do
    end if

    deallocate(col_ptr)
    allocate(col_ptr(n+1))
    do j=1,n+1
       col_ptr(j)=col_ind(j)
    end do
    deallocate(col_ind)

  end subroutine psp_sst_csc2coo

  subroutine psp_sst_fmtCnvt(m,n,nnz,idx,fmt1,fmt2)

    !**** IN ***********************************!
    character(3), intent(in) ::   fmt1, fmt2
    integer, intent(in) ::   m, n, nnz

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: idx(:)

    !**********************************************!

    if (fmt1/=fmt2) then
       select case (fmt1)
       case('coo')
          call psp_sst_coo2csc(m,n,nnz,idx)
       case('csc')
          call psp_sst_csc2coo(m,n,nnz,idx)
       end select
    end if

  end subroutine psp_sst_fmtCnvt

  subroutine psp_coo2csc(spMat)

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout) :: spMat

    !**** INTERNAL ********************************!

    integer :: j, cnt, cnt2, num

    if (spMat%str_type.EQ.'coo') then
       ! Assuming that col_ind is increscent and row_ind is increscent in each column
       spMat%str_type='csc'
       if (allocated(spMat%col_ptr)) deallocate(spMat%col_ptr)
       allocate(spMat%col_ptr(spMat%loc_dim2+1))
       spMat%col_ptr(1)=1
       cnt=1
       do j=1,spMat%loc_dim2
          num=0
          if (cnt<=spMat%nnz) then
             cnt2=0
             do while (spMat%col_ind(cnt+cnt2)==j)
                cnt2=cnt2+1
                num=num+1
                if (cnt+cnt2>spMat%nnz) then
                   exit
                end if
             end do
          end if
          spMat%col_ptr(j+1)=spMat%col_ptr(j)+num
          cnt=cnt+num
       end do
       deallocate(spMat%col_ind)
    end if

  end subroutine psp_coo2csc

  subroutine psp_csc2coo(spMat)

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout) :: spMat

    !**** INTERNAL ********************************!

    integer :: j, k, nst, jlo, jhi, klo, khi

    if (spMat%str_type.EQ.'csc') then
       spMat%str_type='coo'
       if (allocated(spMat%col_ind)) deallocate(spMat%col_ind)
       allocate(spMat%col_ind(spMat%nnz))
       nst=0
       if (spMat%col_ptr(1)==0) then
          jlo=0
          jhi=spMat%loc_dim2-1
          do j=jlo,jhi
             klo=spMat%col_ptr(j+1)
             khi=spMat%col_ptr(j+2)-1
             do k=klo,khi
                nst=nst+1
                spMat%col_ind(nst)=j
             end do
          end do
       else
          jlo=1
          jhi=spMat%loc_dim2
          do j=jlo,jhi
             klo=spMat%col_ptr(j)
             khi=spMat%col_ptr(j+1)-1
             do k=klo,khi
                nst=nst+1
                spMat%col_ind(nst)=j
             end do
          end do
       end if

       deallocate(spMat%col_ptr)
    end if

  end subroutine psp_csc2coo

  subroutine psp_copy_dspm2st(M,N,A,IA,JA,B_idx1,B_idx2,B_val,B_dim1,B_dim2,IB,JB,beta)
    ! B(IB:IB+M-1,JB:JB+N-1) = A(IA:IA+M-1,JA:JA+N-1)
    ! other entries of B are set to be zero
    ! B in a csc format
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, IA, JA, IB, JB, B_dim1, B_dim2
    real(dp), intent(in) :: beta
    type(psp_matrix_spm), intent(in) :: A

    !**** INOUT ***********************************!

    integer, allocatable, intent(inout) :: B_idx1(:), B_idx2(:)
    real(dp), allocatable, intent(inout) :: B_val(:)

    !**** LOCAL ***********************************!

    integer :: i, j, nnz_loc, mk, L, totalMissed, numMissed, row, numOnes
    integer :: IB_1, JB_1, IA_1, JA_1, IA_IB
    integer, allocatable :: tmp_idx1(:)
    real(dp), allocatable :: tmp_val(:)

    IB_1=IB-1
    JB_1=JB-1
    IA_1=IA-1
    JA_1=JA-1
    IA_IB=IA-IB

    ! check sparse format
    if (A%str_type=='coo')  call die('psp_copyspm2st only works for csc format')
    nnz_loc=A%col_ptr(JA+N)-A%col_ptr(JA)
    if (IA==1 .and. M==A%loc_dim1) then
       if (allocated(B_idx1)) deallocate(B_idx1)
       if (allocated(B_val)) deallocate(B_val)
       if (allocated(B_idx2)) deallocate(B_idx2)      ! TODO: any efficient method to avoid this?
       allocate(B_idx1(nnz_loc))
       allocate(B_val(nnz_loc))
       allocate(B_idx2(B_dim2+1))
       mk=1
       numOnes=A%col_ptr(JA)-A%col_ptr(1)
       do i=JA,JA_1+N
          L=A%col_ptr(i+1)-A%col_ptr(i)
          do j=mk,mk+L-1
             B_idx1(j)=A%row_ind(A%col_ptr(i)+j-mk)-IA_IB
             B_val(j)=A%dval(A%col_ptr(i)+j-mk);
          end do
          B_idx2(i-JA_1)=A%col_ptr(i)-numOnes
          mk=mk+L
       end do
       B_idx2(B_dim2+1)=A%col_ptr(JA+N)-numOnes
    else ! In this case, the code is less efficient
       allocate(tmp_idx1(nnz_loc))
       allocate(tmp_val(nnz_loc))
       if (allocated(B_idx2)) deallocate(B_idx2)      ! TODO: any efficient method to avoid this?
       allocate(B_idx2(B_dim2+1))
       ! no efficient row indexing, not efficient in memory, 
       !but efficient in operation complexity
       numOnes=JB_1
       if(numOnes>0) B_idx2(1:numOnes)=1
       mk=1
       totalMissed=A%col_ptr(JA)-A%col_ptr(1)
       do i=JA,JA+N-1
          L=A%col_ptr(i+1)-A%col_ptr(i)
          numMissed=0
          do j=mk,mk+L-1
             row=A%row_ind(A%col_ptr(i)+j-mk)
             if (row<IA.or.row>IA_1+M) then
                numMissed=numMissed+1
             else
                tmp_idx1(j-numMissed)=row-IA_IB
                tmp_val(j-numMissed)=A%dval(A%col_ptr(i)+j-mk)
             end if
          end do
          B_idx2(numOnes+i-JA_1)=A%col_ptr(i)-totalMissed
          totalMissed=totalMissed+numMissed
          mk=mk+L-numMissed
       end do
       do j=1,B_dim2-(JB_1+N)+1
          B_idx2(j+N+JB_1)=A%col_ptr(JA+N)-totalMissed
       end do
       nnz_loc=B_idx2(B_dim2+1)-B_idx2(1)
       if (allocated(B_idx1)) deallocate(B_idx1)
       if (allocated(B_val)) deallocate(B_val)
       allocate(B_idx1(nnz_loc))
       allocate(B_val(nnz_loc))
       B_idx1(1:nnz_loc)=tmp_idx1(1:nnz_loc)
       B_val(1:nnz_loc)=tmp_val(1:nnz_loc)
    end if

    if (allocated(tmp_idx1)) deallocate(tmp_idx1)
    if (allocated(tmp_val)) deallocate(tmp_val)

  end subroutine psp_copy_dspm2st

  subroutine psp_copy_zspm2st(M,N,A,IA,JA,B_idx1,B_idx2,B_val,B_dim1,B_dim2,IB,JB,beta)
    ! B(IB:IB+M-1,JB:JB+N-1) = A(IA:IA+M-1,JA:JA+N-1)
    ! other entries of B are set to be zero
    ! B in a csc format
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, IA, JA, IB, JB, B_dim1, B_dim2
    complex(dp), intent(in) :: beta
    type(psp_matrix_spm), intent(in) :: A

    !**** INOUT ***********************************!

    integer, allocatable, intent(inout) :: B_idx1(:), B_idx2(:)
    complex(dp), allocatable, intent(inout) :: B_val(:)

    !**** LOCAL ***********************************!

    integer :: i, j, nnz_loc, mk, L, totalMissed, numMissed, row, numOnes
    integer :: IB_1, JB_1, IA_1, JA_1, IA_IB
    integer, allocatable :: tmp_idx1(:)
    complex(dp), allocatable :: tmp_val(:)

    IB_1=IB-1
    JB_1=JB-1
    IA_1=IA-1
    JA_1=JA-1
    IA_IB=IA-IB

    ! check sparse format
    if (A%str_type=='coo')  call die('psp_copyspm2st only works for csc format')
    nnz_loc=A%col_ptr(JA+N)-A%col_ptr(JA)
    if (IA==1 .and. M==A%loc_dim1) then
       if (allocated(B_idx1)) deallocate(B_idx1)
       if (allocated(B_val)) deallocate(B_val)
       if (allocated(B_idx2)) deallocate(B_idx2)      ! TODO: any efficient method to avoid this?
       allocate(B_idx1(nnz_loc))
       allocate(B_val(nnz_loc))
       allocate(B_idx2(B_dim2+1))
       mk=1
       numOnes=A%col_ptr(JA)-A%col_ptr(1)
       do i=JA,JA_1+N
          L=A%col_ptr(i+1)-A%col_ptr(i)
          do j=mk,mk+L-1
             B_idx1(j)=A%row_ind(A%col_ptr(i)+j-mk)-IA_IB
             B_val(j)=A%zval(A%col_ptr(i)+j-mk);
          end do
          B_idx2(i-JA_1)=A%col_ptr(i)-numOnes
          mk=mk+L
       end do
       B_idx2(B_dim2+1)=A%col_ptr(JA+N)-numOnes
    else ! In this case, the code is less efficient
       allocate(tmp_idx1(nnz_loc))
       allocate(tmp_val(nnz_loc))
       if (allocated(B_idx2)) deallocate(B_idx2)      ! TODO: any efficient method to avoid this?
       allocate(B_idx2(B_dim2+1))
       ! no efficient row indexing, not efficient in memory,
       !but efficient in operation complexity
       numOnes=JB_1
       if(numOnes>0) B_idx2(1:numOnes)=1
       mk=1
       totalMissed=A%col_ptr(JA)-A%col_ptr(1)
       do i=JA,JA_1+N
          L=A%col_ptr(i+1)-A%col_ptr(i)
          numMissed=0
          do j=mk,mk+L-1
             row=A%row_ind(A%col_ptr(i)+j-mk)
             if (row<IA.or.row>IA+M-1) then
                numMissed=numMissed+1
             else
                tmp_idx1(j-numMissed)=row-IA_IB
                tmp_val(j-numMissed)=A%zval(A%col_ptr(i)+j-mk)
             end if
          end do
          B_idx2(numOnes+i-JA_1)=A%col_ptr(i)-totalMissed
          totalMissed=totalMissed+numMissed
          mk=mk+L-numMissed
       end do
       do j=1,B_dim2-(JB_1+N)+1
          B_idx2(j+N+JB_1)=A%col_ptr(JA+N)-totalMissed
       end do
       nnz_loc=B_idx2(B_dim2+1)-B_idx2(1)
       if (allocated(B_idx1)) deallocate(B_idx1)
       if (allocated(B_val)) deallocate(B_val)
       allocate(B_idx1(nnz_loc))
       allocate(B_val(nnz_loc))
       B_idx1(1:nnz_loc)=tmp_idx1(1:nnz_loc)
       B_val(1:nnz_loc)=tmp_val(1:nnz_loc)
    end if

    if (allocated(tmp_idx1)) deallocate(tmp_idx1)
    if (allocated(tmp_val)) deallocate(tmp_val)

  end subroutine psp_copy_zspm2st

  subroutine psp_copy_dm(M,N,A,IA,JA,B,IB,JB,alpha,beta)
    ! B(IB:IB+M-1,JB:JB+N-1) = alpha*A(IA:IA+M-1,JA:JA+N-1)+beta*B(IB:IB+M-1,JB:JB+N-1)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, IA, JA, IB, JB
    real(dp), intent(in) :: A(:,:), alpha, beta

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: B(:,:)

    !**** LOCAL ***********************************!

    integer :: i, j, IB_1, JB_1, IA_1, JA_1

    IB_1=IB-1
    JB_1=JB-1
    IA_1=IA-1
    JA_1=JA-1

    if (alpha/=0.0_dp) then
       if (beta/=0.0_dp) then
          do j=1,N
             do i=1,M
                B(IB_1+i,JB_1+j)=alpha*A(IA_1+i,JA_1+j)+beta*B(IB_1+i,JB_1+j)
             enddo
          enddo
       else
          do j=1,N
             do i=1,M
                B(IB_1+i,JB_1+j)=alpha*A(IA_1+i,JA_1+j)
             enddo
          enddo
       end if
    else
       if (beta/=0.0_dp) then
          do j=1,N
             do i=1,M
                B(IB_1+i,JB_1+j)=beta*B(IB_1+i,JB_1+j)
             enddo
          enddo
       else
          do j=1,N
             do i=1,M
                B(IB_1+i,JB_1+j)=0.0_dp
             enddo
          enddo
       end if
    end if



  end subroutine psp_copy_dm

  subroutine psp_copy_zm(M,N,A,IA,JA,B,IB,JB,alpha,beta)
    ! B(IB:IB+M-1,JB:JB+N-1) = alpha*A(IA:IA+M-1,JA:JA+N-1)+beta*B(IB:IB+M-1,JB:JB+N-1)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, IA, JA, IB, JB
    complex(dp), intent(in) :: A(:,:), alpha, beta

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: B(:,:)

    !**** LOCAL ***********************************!

    integer :: i, j, IB_1, JB_1, IA_1, JA_1

    IB_1=IB-1
    JB_1=JB-1
    IA_1=IA-1
    JA_1=JA-1

    if (alpha/=cmplx_0) then
       if (beta/=cmplx_0) then
          do j=1,N
             do i=1,M
                B(IB_1+i,JB_1+j)=alpha*A(IA_1+i,JA_1+j)+beta*B(IB_1+i,JB_1+j)
             enddo
          enddo
       else
          do j=1,N
             do i=1,M
                B(IB_1+i,JB_1+j)=alpha*A(IA_1+i,JA_1+j)
             enddo
          enddo
       end if
    else
       if (beta/=cmplx_0) then
          do j=1,N
             do i=1,M
                B(IB_1+i,JB_1+j)=beta*B(IB_1+i,JB_1+j)
             enddo
          enddo
       else
          do j=1,N
             do i=1,M
                B(IB_1+i,JB_1+j)=cmplx_0
             enddo
          enddo
       end if
    end if

  end subroutine psp_copy_zm

  subroutine psp_copy_dv(M,A,IA,B,IB,beta)
    ! B(IB:IB+M-1) = A(IA:IA+M-1)+beta*B(IB:IB+M-1)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, IA, IB
    real(dp), intent(in) :: A(:), beta

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: B(:)

    !**** LOCAL ***********************************!

    integer :: i, IB_1, IA_1

    IB_1=IB-1
    IA_1=IA-1

    if (beta/=0.0_dp) then
       do i=1,M
          B(IB_1+i)=A(IA_1+i)+beta*B(IB_1+i)
       enddo
    else
       do i=1,M
          B(IB_1+i)=A(IA_1+i)
       enddo
    end if

  end subroutine psp_copy_dv

  subroutine psp_copy_zv(M,A,IA,B,IB,beta)
    ! B(IB:IB+M-1) = A(IA:IA+M-1)+beta*B(IB:IB+M-1)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, IA, IB
    complex(dp), intent(in) :: A(:), beta

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: B(:)

    !**** LOCAL ***********************************!

    integer :: i, IB_1, IA_1

    IB_1=IB-1
    IA_1=IA-1

    if (beta/=cmplx_0) then
       do i=1,M
          B(IB_1+i)=A(IA_1+i)+beta*B(IB_1+i)
       enddo
    else
       do i=1,M
          B(IB_1+i)=A(IA_1+i)
       enddo
    end if

  end subroutine psp_copy_zv

  subroutine psp_idx_glb2loc(glb,bs,npproc,loc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: glb, bs, npproc
    !**** INOUT ***********************************!

    integer, intent(inout) :: loc

    !**** LOCAL ***********************************!

    integer :: num_per_cycle, idx_cycle, idx_in_blk

    num_per_cycle = bs*npproc
    idx_cycle = ceiling(1.0_dp*glb/num_per_cycle) !glb is in the idx_cycle'th cycle
    idx_in_blk = mod(glb-1,num_per_cycle)+1
    idx_in_blk = mod(idx_in_blk-1,bs)+1 !glb is in the idx_in_blk'th position of the idx_cycle'th cycle
    loc = (idx_cycle-1)*bs + idx_in_blk

  end subroutine psp_idx_glb2loc

  subroutine psp_process_lopM(opM,trM)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opM

    !**** INOUT ***********************************!

    logical, intent(inout) :: trM

    !**********************************************!

    if ((opM .eq. 'T') .or. &
         (opM .eq. 't') .or. &
         (opM .eq. 'C') .or. &
         (opM .eq. 'c')) then
       trM=.true.
    else if ((opM .eq. 'N') .or. &
         (opM .eq. 'n')) then
       trM=.false.
    else
       call die('process_lopM: invalid opM')
    end if

  end subroutine psp_process_lopM

  subroutine psp_process_iopM(opM,tcM)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opM

    !**** INOUT ***********************************!

    integer, intent(inout) :: tcM

    !**********************************************!

    if ((opM .eq. 'T') .or. &
         (opM .eq. 't')) then
       tcM=2
    else if ((opM .eq. 'C') .or. &
         (opM .eq. 'c')) then
       tcM=1
    else if ((opM .eq. 'N') .or. &
         (opM .eq. 'n')) then
       tcM=0
    else
       call die('process_iopM: invalid opM')
    end if

  end subroutine psp_process_iopM

  subroutine init_random_seed()
    use iso_fortran_env, only: int64
    implicit none

#ifdef MPI
    include 'mpif.h'
#endif

    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid, mpi_err
    integer(int64) :: t

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
       call mpi_comm_rank(mpi_comm_world,pid,mpi_err)
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
  end subroutine init_random_seed


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


END MODULE pspBasicTool
