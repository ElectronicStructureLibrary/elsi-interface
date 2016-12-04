MODULE pspMspm_nn
  use pspVariable
  use pspUtility
  use pspMPI
  use pspLevel1
  use pspLevel2

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

  interface psp_gemspm_nn
     module procedure psp_dgemspm_nn
     module procedure psp_zgemspm_nn
  end interface psp_gemspm_nn

  interface die
     module procedure die
  end interface die

  public :: psp_gemspm_nn

contains

  !================================================!
  !    sparse pdgemm: mspm  C=alpha*A*B+beta*C     !
  !================================================!

  subroutine psp_dgemspm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, A of size M by K, B of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    real(dp), intent(in)  ::   alpha, beta ! scalar
    real(dp), intent(in) :: A(:,:)
    type(psp_matrix_spm), intent(in) ::   B ! matrix B

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: C(:,:)

    !**** LOCAL ***********************************!

    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz, i, idx_col
    integer :: trA, trB
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    real(dp), allocatable :: A_loc(:,:), C_loc(:,:), CC_loc(:,:)
    real(dp), allocatable :: B_loc_val(:)
    integer, allocatable :: B_loc_idx1(:), B_loc_idx2(:)
    integer :: A_loc_dim(2), B_loc_dim(2), C_loc_dim(2), coords(2)
    integer :: width, glb_st, glb_ed, loc_st, loc_ed, nnz_loc, L

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

    !**********************************************!

    integer, external :: numroc ! it is a function to compute local size

    !**********************************************!

    ! obtain grid information, working on the (iprow,ipcol) processor in a grid of size nprow by npcol
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)

    ! set up local variables
    ! allocate local matrices
    A_loc_dim(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
    A_loc_dim(2)=psp_update_rank
    allocate(A_loc(A_loc_dim(1),A_loc_dim(2)))
    A_loc=0.0_dp
    B_loc_dim(1)=psp_update_rank
    B_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    C_loc_dim(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
    C_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    allocate(C_loc(C_loc_dim(1),C_loc_dim(2)))
    C_loc=0.0_dp

    !**** Start the operation ********************!

    ! begin SUMMA
    do kloop = 1, K, psp_update_rank
       idx_k_row = (kloop-1)/psp_bs_def_row+1
       idx_k_col = (kloop-1)/psp_bs_def_col+1
       glb_st=kloop
       glb_ed=min(kloop+psp_update_rank,K+1)-1
       width=glb_ed-glb_st+1
       idx_pcol = mod(idx_k_col-1,npcol) ! identify the processor owning A(:,kth block), the cart coordinate
       if (ipcol==idx_pcol) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_col,npcol,loc_st)
          if (width<psp_update_rank) then
             A_loc=0.0_dp
          endif
          call psp_copy_m(A_loc_dim(1),width,A,1,loc_st,A_loc,1,1,1.0_dp,0.0_dp)
       end if

       ! boardcast in row
       call MPI_Bcast(A_loc, A_loc_dim(1)*A_loc_dim(2), MPI_DOUBLE, idx_pcol, psp_mpi_comm_row,mpi_err)

       idx_prow = mod(idx_k_row-1,nprow) ! identify the processor owing B(kth block,:), the cart coordinate
       if (iprow==idx_prow) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_row,nprow,loc_st)
          call psp_copy_spm2st(width,B_loc_dim(2),B,loc_st,1,B_loc_idx1,B_loc_idx2, &
               B_loc_val,width,B_loc_dim(2),1,1,0.0_dp)
       else
          if (allocated(B_loc_idx1)) deallocate(B_loc_idx1)
          if (allocated(B_loc_idx2)) deallocate(B_loc_idx2)
          if (allocated(B_loc_val)) deallocate(B_loc_val)
          allocate(B_loc_idx2(B_loc_dim(2)+1))
          allocate(B_loc_idx1(1))
          allocate(B_loc_val(1))
       end if

       ! boardcast in column
       if (psp_update_rank>1) then
          call MPI_Bcast(width, 1, MPI_INT, idx_prow, psp_mpi_comm_col,mpi_err)
       end if
       L=size(B_loc_idx2)
       ! data to be boardcast have the same size in every processor
       if (L/=B_loc_dim(2)+1) then
          if (allocated(B_loc_idx2)) deallocate(B_loc_idx2)
          allocate(B_loc_idx2(B_loc_dim(2)+1))
       end if
       call MPI_Bcast(B_loc_idx2, B_loc_dim(2)+1, MPI_INT, idx_prow, psp_mpi_comm_col,mpi_err)
       nnz_loc=B_loc_idx2(B_loc_dim(2)+1)-B_loc_idx2(1)
       L=size(B_loc_idx1)
       if (L/=nnz_loc) then
          if (allocated(B_loc_idx1)) deallocate(B_loc_idx1)
          if (allocated(B_loc_val)) deallocate(B_loc_val)
          allocate(B_loc_idx1(nnz_loc))
          allocate(B_loc_val(nnz_loc))
       end if
       call MPI_Bcast(B_loc_idx1, nnz_loc, MPI_INT, idx_prow, psp_mpi_comm_col,mpi_err)
       call MPI_Bcast(B_loc_val, nnz_loc, MPI_DOUBLE, idx_prow, psp_mpi_comm_col,mpi_err)

       ! compute local update of C
       ! C=MATMUL(A_loc,B_loc)+C
       do idx_col=1,C_loc_dim(2)
          do i=B_loc_idx2(idx_col),B_loc_idx2(idx_col+1)-1
             if (B_loc_idx1(i)<=width) then
                C_loc(1:C_loc_dim(1),idx_col)= B_loc_val(i)*A_loc(1:C_loc_dim(1),B_loc_idx1(i)) &
                     + C_loc(1:C_loc_dim(1),idx_col)
             end if
          end do
       end do
    enddo

    !C=beta*C+alpha*C_loc
    call psp_copy_m(C_loc_dim(1),C_loc_dim(2),C_loc,1,1,C,1,1,alpha,beta)

    if (allocated(A_loc)) deallocate(A_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)
    if (allocated(B_loc_val)) deallocate(B_loc_val)
    if (allocated(B_loc_idx1)) deallocate(B_loc_idx1)
    if (allocated(B_loc_idx2)) deallocate(B_loc_idx2)

  end subroutine psp_dgemspm_nn

  !================================================!
  !    sparse pzgemm: mspm  C=alpha*A*B+beta*C     !
  !================================================!

  subroutine psp_zgemspm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, A of size M by K, B of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    complex(dp), intent(in)  ::   alpha, beta ! scalar
    complex(dp), intent(in) :: A(:,:)
    type(psp_matrix_spm), intent(in) ::   B ! matrix B

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: C(:,:)

    !**** LOCAL ***********************************!

    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz, i, idx_col
    integer :: trA, trB
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    complex(dp), allocatable :: A_loc(:,:), C_loc(:,:), CC_loc(:,:)
    complex(dp), allocatable :: B_loc_val(:)
    integer, allocatable :: B_loc_idx1(:), B_loc_idx2(:)
    integer :: A_loc_dim(2), B_loc_dim(2), C_loc_dim(2), coords(2)
    integer :: width, glb_st, glb_ed, loc_st, loc_ed, nnz_loc, L

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

    !**********************************************!

    integer, external :: numroc ! it is a function to compute local size

    !**********************************************!

    ! obtain grid information, working on the (iprow,ipcol) processor in a grid of size nprow by npcol
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)

    ! set up local variables
    ! allocate local matrices
    A_loc_dim(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
    A_loc_dim(2)=psp_update_rank
    allocate(A_loc(A_loc_dim(1),A_loc_dim(2)))
    A_loc=cmplx_0
    B_loc_dim(1)=psp_update_rank
    B_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    C_loc_dim(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
    C_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    allocate(C_loc(C_loc_dim(1),C_loc_dim(2)))
    C_loc=cmplx_0

    !**** Start the operation ********************!

    ! begin SUMMA
    do kloop = 1, K, psp_update_rank
       idx_k_row = (kloop-1)/psp_bs_def_row+1
       idx_k_col = (kloop-1)/psp_bs_def_col+1
       glb_st=kloop
       glb_ed=min(kloop+psp_update_rank,K+1)-1
       width=glb_ed-glb_st+1
       idx_pcol = mod(idx_k_col-1,npcol) ! identify the processor owning A(:,kth block), the cart coordinate
       if (ipcol==idx_pcol) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_col,npcol,loc_st)
          if (width<psp_update_rank) then
             A_loc=cmplx_0
          endif
          call psp_copy_m(A_loc_dim(1),width,A,1,loc_st,A_loc,1,1,cmplx_1,cmplx_0)
       end if

       ! boardcast in row
       call MPI_Bcast(A_loc, A_loc_dim(1)*A_loc_dim(2), MPI_DOUBLE_COMPLEX, idx_pcol, psp_mpi_comm_row,mpi_err)

       idx_prow = mod(idx_k_row-1,nprow) ! identify the processor owing B(kth block,:), the cart coordinate
       if (iprow==idx_prow) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_row,nprow,loc_st)
          call psp_copy_spm2st(width,B_loc_dim(2),B,loc_st,1,B_loc_idx1,B_loc_idx2, &
               B_loc_val,width,B_loc_dim(2),1,1,cmplx_0)
       else
          if (allocated(B_loc_idx1)) deallocate(B_loc_idx1)
          if (allocated(B_loc_idx2)) deallocate(B_loc_idx2)
          if (allocated(B_loc_val)) deallocate(B_loc_val)
          allocate(B_loc_idx2(B_loc_dim(2)+1))
          allocate(B_loc_idx1(1))
          allocate(B_loc_val(1))
       end if

       ! boardcast in column
       if (psp_update_rank>1) then
          call MPI_Bcast(width, 1, MPI_INT, idx_prow, psp_mpi_comm_col,mpi_err)
       end if
       L=size(B_loc_idx2)
       ! data to be boardcast have the same size in every processor
       if (L/=B_loc_dim(2)+1) then
          if (allocated(B_loc_idx2)) deallocate(B_loc_idx2)
          allocate(B_loc_idx2(B_loc_dim(2)+1))
       end if
       call MPI_Bcast(B_loc_idx2, B_loc_dim(2)+1, MPI_INT, idx_prow, psp_mpi_comm_col,mpi_err)
       nnz_loc=B_loc_idx2(B_loc_dim(2)+1)-B_loc_idx2(1)
       L=size(B_loc_idx1)
       if (L/=nnz_loc) then
          if (allocated(B_loc_idx1)) deallocate(B_loc_idx1)
          if (allocated(B_loc_val)) deallocate(B_loc_val)
          allocate(B_loc_idx1(nnz_loc))
          allocate(B_loc_val(nnz_loc))
       end if
       call MPI_Bcast(B_loc_idx1, nnz_loc, MPI_INT, idx_prow, psp_mpi_comm_col,mpi_err)
       call MPI_Bcast(B_loc_val, nnz_loc, MPI_DOUBLE_COMPLEX, idx_prow, psp_mpi_comm_col,mpi_err)

       ! compute local update of C
       ! C=MATMUL(A_loc,B_loc)+C
       do idx_col=1,C_loc_dim(2)
          do i=B_loc_idx2(idx_col),B_loc_idx2(idx_col+1)-1
             if (B_loc_idx1(i)<=width) then
                C_loc(1:C_loc_dim(1),idx_col)= B_loc_val(i)*A_loc(1:C_loc_dim(1),B_loc_idx1(i)) &
                     + C_loc(1:C_loc_dim(1),idx_col)
             end if
          end do
       end do
    enddo

    !C=beta*C+alpha*C_loc
    call psp_copy_m(C_loc_dim(1),C_loc_dim(2),C_loc,1,1,C,1,1,alpha,beta)

    if (allocated(A_loc)) deallocate(A_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)
    if (allocated(B_loc_val)) deallocate(B_loc_val)
    if (allocated(B_loc_idx1)) deallocate(B_loc_idx1)
    if (allocated(B_loc_idx2)) deallocate(B_loc_idx2)

  end subroutine psp_zgemspm_nn

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


END MODULE pspMspm_nn
