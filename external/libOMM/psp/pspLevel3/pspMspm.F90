MODULE pspMspm
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

  interface psp_gemspm
     module procedure psp_dgemspm
     module procedure psp_zgemspm
  end interface psp_gemspm

  interface psp_gemspm_nn
     module procedure psp_dgemspm_nn
     module procedure psp_zgemspm_nn
  end interface psp_gemspm_nn

  interface psp_gemspm_nt
     module procedure psp_dgemspm_nt
     module procedure psp_zgemspm_nt
  end interface psp_gemspm_nt

  interface psp_gemspm_tn
     module procedure psp_dgemspm_tn
     module procedure psp_zgemspm_tn
  end interface psp_gemspm_tn

  interface psp_gemspm_tt
     module procedure psp_dgemspm_tt
     module procedure psp_zgemspm_tt
  end interface psp_gemspm_tt

  interface die
     module procedure die
  end interface die

  public :: psp_gemspm

contains

  subroutine psp_dgemspm(M,N,K,A,opA,B,opB,C,alpha,beta)
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

    integer :: trA, trB, ot

    !**********************************************!
    if (alpha/=0.0_dp) then
       call psp_process_opM(opA,trA)
       call psp_process_opM(opB,trB)
       ! operation table
       if (trA==0 .and. trB==0) then
          ot=1
       else if (trA==0 .and. trB>=1) then
          ot=2
       else if (trA>=1 .and. trB==0) then
          ot=3
       else if (trA>=1 .and. trB>=1) then
          ot=4
       else
          call die('mm_dmultiply: invalid implementation')
       end if

       select case (ot)
       case (1)
          call psp_gemspm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (2)
          call psp_gemspm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (3)
          call psp_gemspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (4)
          call psp_gemspm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
       end select
    else
       if (beta/=0.0_dp) C=beta*C
    end if

  end subroutine psp_dgemspm


  subroutine psp_zgemspm(M,N,K,A,opA,B,opB,C,alpha,beta)
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

    integer :: trA, trB, ot

    !**********************************************!
    if (alpha/=cmplx_0) then
       call psp_process_opM(opA,trA)
       call psp_process_opM(opB,trB)
       ! operation table
       if (trA==0 .and. trB==0) then
          ot=1
       else if (trA==0 .and. trB>=1) then
          ot=2
       else if (trA>=1 .and. trB==0) then
          ot=3
       else if (trA>=1 .and. trB>=1) then
          ot=4
       else
          call die('mm_dmultiply: invalid implementation')
       end if

       select case (ot)
       case (1)
          call psp_gemspm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (2)
          call psp_gemspm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (3)
          call psp_gemspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (4)
          call psp_gemspm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
       end select
    else
       if (beta/=cmplx_0) C=beta*C
    end if
  end subroutine psp_zgemspm

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


    common /coeff/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /coeff/ psp_bs_def_row, psp_bs_def_col
    common /coeff/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /coeff/ psp_bs_num
    common /coeff/ psp_icontxt ! BLACS context handle used by psp
    common /coeff/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
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
  !    sparse pdgemm: mspm  C=alpha*A*B'+beta*C    !
  !================================================!

  subroutine psp_dgemspm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)
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

    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz
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


    common /coeff/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /coeff/ psp_bs_def_row, psp_bs_def_col
    common /coeff/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /coeff/ psp_bs_num
    common /coeff/ psp_icontxt ! BLACS context handle used by psp
    common /coeff/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
#endif

    !**********************************************!

    integer, external :: numroc ! it is a function to compute local size

    !**********************************************!

    ! obtain grid information, working on the (iprow,ipcol) processor in a grid of size nprow by npcol
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
    ! set up local variables
    ! allocate local matrices
    A_loc_dim(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
    A_loc_dim(2)=numroc(K,psp_bs_def_col,ipcol,0,npcol)
    B_loc_dim(1)=psp_update_rank
    B_loc_dim(2)=numroc(K,psp_bs_def_col,ipcol,0,npcol)
    C_loc_dim(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
    C_loc_dim(2)=psp_update_rank
    allocate(C_loc(C_loc_dim(1),C_loc_dim(2)))
    C_loc=0.0_dp
    allocate(CC_loc(C_loc_dim(1),C_loc_dim(2)))
    CC_loc=0.0_dp

    !**** Start the operation ********************!

    ! begin SUMMA
    do kloop = 1, N, psp_update_rank
       idx_k_row = (kloop-1)/psp_bs_def_row+1
       idx_k_col = (kloop-1)/psp_bs_def_col+1
       glb_st=kloop
       glb_ed=min(kloop+psp_update_rank,N+1)-1
       width=glb_ed-glb_st+1
       idx_prow = mod(idx_k_row-1,nprow)
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
       ! C_loc = A*(B_loc^t)
       call psp_sst_gemspm(C_loc_dim(1),width,A_loc_dim(2),opA,opB, &
            1.0_dp,A,1,1,B_loc_idx1,B_loc_idx2,B_loc_val,C_loc,1,1,0.0_dp)

       idx_pcol = mod(idx_k_col-1,npcol ) ! identify the processor owing C(:,kth block,), the cart coordinate
       ! boardcast in row
       call MPI_REDUCE(C_loc, CC_loc, C_loc_dim(1)*C_loc_dim(2), MPI_DOUBLE, MPI_SUM, idx_pcol, psp_mpi_comm_row, mpi_err)

       if (ipcol==idx_pcol) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_col,npcol,loc_st)
          !C=beta*C+C_loc
          call psp_copy_m(C_loc_dim(1),width,CC_loc,1,1,C,1,loc_st,alpha,beta)
       end if
    enddo
    if (allocated(A_loc)) deallocate(A_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)
    if (allocated(B_loc_val)) deallocate(B_loc_val)
    if (allocated(B_loc_idx1)) deallocate(B_loc_idx1)
    if (allocated(B_loc_idx2)) deallocate(B_loc_idx2)

  end subroutine psp_dgemspm_nt

  !================================================!
  !    sparse pdgemm: mspm  C=alpha*A'*B+beta*C    !
  !================================================!

  subroutine psp_dgemspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
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

    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz
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


    common /coeff/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /coeff/ psp_bs_def_row, psp_bs_def_col
    common /coeff/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /coeff/ psp_bs_num
    common /coeff/ psp_icontxt ! BLACS context handle used by psp
    common /coeff/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
#endif

    !**********************************************!

    integer, external :: numroc ! it is a function to compute local size

    !**********************************************!

    ! obtain grid information, working on the (iprow,ipcol) processor in a grid of size nprow by npcol
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
    ! set up local variables
    ! allocate local matrices
    A_loc_dim(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
    A_loc_dim(2)=psp_update_rank
    allocate(A_loc(A_loc_dim(1),A_loc_dim(2)))
    A_loc=0.0_dp
    B_loc_dim(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
    B_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    C_loc_dim(1)=psp_update_rank
    C_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    allocate(C_loc(C_loc_dim(1),C_loc_dim(2)))
    C_loc=0.0_dp
    allocate(CC_loc(C_loc_dim(1),C_loc_dim(2)))
    CC_loc=0.0_dp

    !**** Start the operation ********************!

    ! begin SUMMA
    do kloop = 1, M, psp_update_rank
       idx_k_row = (kloop-1)/psp_bs_def_row+1
       idx_k_col = (kloop-1)/psp_bs_def_col+1
       glb_st=kloop
       glb_ed=min(kloop+psp_update_rank,M+1)-1
       width=glb_ed-glb_st+1
       idx_pcol = mod(idx_k_col-1,npcol)
       if (ipcol==idx_pcol) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_col,npcol,loc_st)
          if (width<psp_update_rank) then
             A_loc=0.0_dp
          endif
          call psp_copy_m(A_loc_dim(1),width,A,1,loc_st,A_loc,1,1,1.0_dp,0.0_dp)
       end if

       ! boardcast in row
       call MPI_Bcast(A_loc, A_loc_dim(1)*A_loc_dim(2), MPI_DOUBLE, idx_pcol, psp_mpi_comm_row,mpi_err)

       ! compute local update of C
       ! C_loc = (A_loc)^t*B
       call psp_sst_gemspm(width,C_loc_dim(2),B_loc_dim(1),opA,opB, &
            1.0_dp,A_loc,1,1,B%row_ind,B%col_ptr,B%dval,C_loc,1,1,0.0_dp)

       idx_prow = mod(idx_k_row-1,nprow ) ! identify the processor owing C(kth block,:), the cart coordinate
       ! boardcast in column
       call MPI_REDUCE(C_loc, CC_loc, C_loc_dim(1)*C_loc_dim(2), MPI_DOUBLE, MPI_SUM, idx_prow, psp_mpi_comm_col, mpi_err)

       if (iprow==idx_prow) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_row,nprow,loc_st)
          !C=beta*C+C_loc
          call psp_copy_m(width,C_loc_dim(2),CC_loc,1,1,C,loc_st,1,alpha,beta)
       end if
    enddo

    if (allocated(A_loc)) deallocate(A_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)
    if (allocated(B_loc_val)) deallocate(B_loc_val)
    if (allocated(B_loc_idx1)) deallocate(B_loc_idx1)
    if (allocated(B_loc_idx2)) deallocate(B_loc_idx2)

  end subroutine psp_dgemspm_tn

  !================================================!
  !   sparse pdgemm: mspm  C=alpha*A'*B'+beta*C    !
  !================================================!

  subroutine psp_dgemspm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
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

    integer :: iprow, ipcol, nprow, npcol, mpi_err
    real(dp), allocatable :: tmp(:,:)
    integer :: dims_before(2), dims_after(2)
    integer :: desc_before(9), desc_after(9)

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

    !**********************************************!

    integer, external :: numroc ! it is a function to compute local size

    !**********************************************!

    ! obtain grid information, working on the (iprow,ipcol) processor in a grid of size nprow by npcol
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)

    ! TODO: optimize the following code by introducing a transpose operator for sparse matrices
    if (.true.) then
       ! tmp = transpose(A)
       dims_before(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
       dims_before(2)=numroc(M,psp_bs_def_col,ipcol,0,npcol)
       call descinit(desc_before,K,M,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_before(1),mpi_err)
       dims_after(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
       dims_after(2)=numroc(K,psp_bs_def_col,ipcol,0,npcol)
       call descinit(desc_after,M,K,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_after(1),mpi_err)
       allocate(tmp(dims_after(1),dims_after(2)))
       call pdtran(M,K,1.0_dp,A,1,1,desc_before, 0.0_dp,tmp,1,1,desc_after)

       ! compute C = alpha*tmp*transpose(B) + beta*C
       call psp_dgemspm_nt(M,N,K,tmp,'n',B,opB,C,alpha,beta)
    end if

    if (allocated(tmp)) deallocate(tmp)

  end subroutine psp_dgemspm_tt

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


    common /coeff/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /coeff/ psp_bs_def_row, psp_bs_def_col
    common /coeff/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /coeff/ psp_bs_num
    common /coeff/ psp_icontxt ! BLACS context handle used by psp
    common /coeff/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
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
       do idx_col=1,C_loc_dim(2)
          if (B_loc_idx2(idx_col)<B_loc_idx2(idx_col+1)) then
             C_loc(1:C_loc_dim(1),idx_col)= B_loc_val(B_loc_idx2(idx_col))*A_loc(1:C_loc_dim(1),1) &
                  + C_loc(1:C_loc_dim(1),idx_col)
          end if
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

  !================================================!
  !    sparse pzgemm: mspm  C=alpha*A*B'+beta*C    !
  !================================================!

  subroutine psp_zgemspm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)
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

    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz
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


    common /coeff/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /coeff/ psp_bs_def_row, psp_bs_def_col
    common /coeff/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /coeff/ psp_bs_num
    common /coeff/ psp_icontxt ! BLACS context handle used by psp
    common /coeff/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
#endif

    !**********************************************!

    integer, external :: numroc ! it is a function to compute local size

    !**********************************************!

    ! obtain grid information, working on the (iprow,ipcol) processor in a grid of size nprow by npcol
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
    ! set up local variables
    ! allocate local matrices
    A_loc_dim(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
    A_loc_dim(2)=numroc(K,psp_bs_def_col,ipcol,0,npcol)
    B_loc_dim(1)=psp_update_rank
    B_loc_dim(2)=numroc(K,psp_bs_def_col,ipcol,0,npcol)
    C_loc_dim(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
    C_loc_dim(2)=psp_update_rank
    allocate(C_loc(C_loc_dim(1),C_loc_dim(2)))
    C_loc=cmplx_0
    allocate(CC_loc(C_loc_dim(1),C_loc_dim(2)))
    CC_loc=cmplx_0

    !**** Start the operation ********************!

    ! begin SUMMA
    do kloop = 1, N, psp_update_rank
       idx_k_row = (kloop-1)/psp_bs_def_row+1
       idx_k_col = (kloop-1)/psp_bs_def_col+1
       glb_st=kloop
       glb_ed=min(kloop+psp_update_rank,N+1)-1
       width=glb_ed-glb_st+1
       idx_prow = mod(idx_k_row-1,nprow)
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
       ! C_loc = A*(B_loc^t)
       call psp_sst_gemspm(C_loc_dim(1),width,A_loc_dim(2),opA,opB, &
            cmplx_1,A,1,1,B_loc_idx1,B_loc_idx2,B_loc_val,C_loc,1,1,cmplx_0)

       idx_pcol = mod(idx_k_col-1,npcol ) ! identify the processor owing C(:,kth block,), the cart coordinate
       ! boardcast in row
       !CC_loc=cmplx_0 ! not necessary
       call MPI_REDUCE(C_loc, CC_loc, C_loc_dim(1)*C_loc_dim(2), MPI_DOUBLE_COMPLEX, MPI_SUM, idx_pcol, psp_mpi_comm_row, mpi_err)
       ! Use CC_loc because cannot reduce to the same C_loc.

       if (ipcol==idx_pcol) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_col,npcol,loc_st)
          !C=beta*C+C_loc
          call psp_copy_m(C_loc_dim(1),width,CC_loc,1,1,C,1,loc_st,alpha,beta)
       end if
    enddo
    if (allocated(A_loc)) deallocate(A_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)
    if (allocated(B_loc_val)) deallocate(B_loc_val)
    if (allocated(B_loc_idx1)) deallocate(B_loc_idx1)
    if (allocated(B_loc_idx2)) deallocate(B_loc_idx2)

  end subroutine psp_zgemspm_nt

  !================================================!
  !    sparse pzgemm: mspm  C=alpha*A'*B+beta*C    !
  !================================================!

  subroutine psp_zgemspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
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

    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz
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


    common /coeff/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /coeff/ psp_bs_def_row, psp_bs_def_col
    common /coeff/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /coeff/ psp_bs_num
    common /coeff/ psp_icontxt ! BLACS context handle used by psp
    common /coeff/ psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col
#endif

    !**********************************************!

    integer, external :: numroc ! it is a function to compute local size

    !**********************************************!

    ! obtain grid information, working on the (iprow,ipcol) processor in a grid of size nprow by npcol
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
    ! set up local variables
    ! allocate local matrices
    A_loc_dim(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
    A_loc_dim(2)=psp_update_rank
    allocate(A_loc(A_loc_dim(1),A_loc_dim(2)))
    A_loc=cmplx_0
    B_loc_dim(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
    B_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    C_loc_dim(1)=psp_update_rank
    C_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    allocate(C_loc(C_loc_dim(1),C_loc_dim(2)))
    C_loc=cmplx_0
    allocate(CC_loc(C_loc_dim(1),C_loc_dim(2)))
    CC_loc=cmplx_0

    !**** Start the operation ********************!

    ! begin SUMMA
    do kloop = 1, M, psp_update_rank
       idx_k_row = (kloop-1)/psp_bs_def_row+1
       idx_k_col = (kloop-1)/psp_bs_def_col+1
       glb_st=kloop
       glb_ed=min(kloop+psp_update_rank,M+1)-1
       width=glb_ed-glb_st+1
       idx_pcol = mod(idx_k_col-1,npcol)
       if (ipcol==idx_pcol) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_col,npcol,loc_st)
          if (width<psp_update_rank) then
             A_loc=cmplx_0
          endif
          call psp_copy_m(A_loc_dim(1),width,A,1,loc_st,A_loc,1,1,cmplx_1,cmplx_0)
       end if

       ! boardcast in row
       call MPI_Bcast(A_loc, A_loc_dim(1)*A_loc_dim(2), MPI_DOUBLE_COMPLEX, idx_pcol, psp_mpi_comm_row,mpi_err)

       ! compute local update of C
       ! C_loc = (A_loc)^t*B
       call psp_sst_gemspm(width,C_loc_dim(2),B_loc_dim(1),opA,opB, &
            cmplx_1,A_loc,1,1,B%row_ind,B%col_ptr,B%zval,C_loc,1,1,cmplx_0)

       idx_prow = mod(idx_k_row-1,nprow ) ! identify the processor owing C(kth block,:), the cart coordinate
       ! boardcast in column
       call MPI_REDUCE(C_loc, CC_loc, C_loc_dim(1)*C_loc_dim(2), MPI_DOUBLE_COMPLEX, MPI_SUM, idx_prow, psp_mpi_comm_col, mpi_err)

       if (iprow==idx_prow) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_row,nprow,loc_st)
          !C=beta*C+C_loc
          call psp_copy_m(width,C_loc_dim(2),CC_loc,1,1,C,loc_st,1,alpha,beta)
       end if
    enddo

    if (allocated(A_loc)) deallocate(A_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)
    if (allocated(B_loc_val)) deallocate(B_loc_val)
    if (allocated(B_loc_idx1)) deallocate(B_loc_idx1)
    if (allocated(B_loc_idx2)) deallocate(B_loc_idx2)

  end subroutine psp_zgemspm_tn

  !================================================!
  !   sparse pzgemm: mspm  C=alpha*A'*B'+beta*C    !
  !================================================!

  subroutine psp_zgemspm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
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

    integer :: iprow, ipcol, nprow, npcol, mpi_err
    complex(dp), allocatable :: tmp(:,:), tmpConj(:,:)
    integer :: dims_before(2), dims_after(2), dims_conj(2)
    integer :: desc_before(9), desc_after(9)
    integer :: trA, trB

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

    !**********************************************!

    integer, external :: numroc ! it is a function to compute local size

    !**********************************************!

    ! obtain grid information, working on the (iprow,ipcol) processor in a grid of size nprow by npcol
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)

    ! TODO: optimize the following code by introducing a transpose operator for sparse matrices
    if (.true.) then
       ! tmp = transpose(A)
       dims_before(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
       dims_before(2)=numroc(M,psp_bs_def_col,ipcol,0,npcol)
       call descinit(desc_before,K,M,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_before(1),mpi_err)
       dims_after(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
       dims_after(2)=numroc(K,psp_bs_def_col,ipcol,0,npcol)
       call descinit(desc_after,M,K,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_after(1),mpi_err)
       allocate(tmp(dims_after(1),dims_after(2)))
       call psp_process_opM(opA,trA)
       if (trA==1) then
          call pztranc(M,K,cmplx_1,A,1,1,desc_before,cmplx_0,tmp,1,1,desc_after)
       else
          call pztranu(M,K,cmplx_1,A,1,1,desc_before,cmplx_0,tmp,1,1,desc_after)
       end if

       ! compute C = alpha*tmp*transpose(B) + beta*C
       call psp_zgemspm_nt(M,N,K,tmp,'n',B,opB,C,alpha,beta)
    end if

    if (allocated(tmp)) deallocate(tmp)

  end subroutine psp_zgemspm_tt


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


END MODULE pspMspm
