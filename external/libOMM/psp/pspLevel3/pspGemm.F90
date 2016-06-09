MODULE pspGemm
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

  interface psp_gemm
     module procedure psp_dgemm
     module procedure psp_zgemm
  end interface psp_gemm

  interface psp_gemm_nn
     module procedure psp_dgemm_nn
     module procedure psp_zgemm_nn
  end interface psp_gemm_nn

  interface psp_gemm_nt
     module procedure psp_dgemm_nt
     module procedure psp_zgemm_nt
  end interface psp_gemm_nt

  interface psp_gemm_tn
     module procedure psp_dgemm_tn
     module procedure psp_zgemm_tn
  end interface psp_gemm_tn

  interface psp_gemm_tt
     module procedure psp_dgemm_tt
     module procedure psp_zgemm_tt
  end interface psp_gemm_tt

  interface die
     module procedure die
  end interface die

  public :: psp_gemm

contains

  !================================================!
  !             pdgemm for SUMMA                   !
  !================================================!
  subroutine psp_dgemm(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, opA(A) of size M by K, opB(B) of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    real(dp), intent(in)  ::   alpha, beta ! scalar
    real(dp), intent(in) ::   A(:,:) ! matrix A
    real(dp), intent(in) ::   B(:,:) ! matrix B

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: C(:,:) ! matrix C

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
          call psp_gemm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (2)
          call psp_gemm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (3)
          call psp_gemm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (4)
          call psp_gemm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
       end select
    else
       if (beta/=0.0_dp) C=beta*C
    end if

  end subroutine psp_dgemm

  !================================================!
  !             pzgemm for SUMMA                   !
  !================================================!

  subroutine psp_zgemm(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, opA(A) of size M by K, opB(B) of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    complex(dp), intent(in)  ::   alpha, beta ! scalar
    complex(dp), intent(in) ::   A(:,:) ! matrix A
    complex(dp), intent(in) ::   B(:,:) ! matrix B

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: C(:,:) ! matrix C

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
          call psp_gemm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (2)
          call psp_gemm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (3)
          call psp_gemm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (4)
          call psp_gemm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
       end select
    else
       if (beta/=cmplx_0) C=beta*C
    end if

  end subroutine psp_zgemm

  !================================================!
  !     pdgemm for SUMMA  C=alpha*A*B+beta*C       !
  !================================================!

  subroutine psp_dgemm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, opA(A) of size M by K, opB(B) of size K by N
    real(dp), intent(in)  ::   alpha, beta ! scalar
    real(dp), intent(in) ::   A(:,:) ! matrix A
    real(dp), intent(in) ::   B(:,:) ! matrix B
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: C(:,:) ! matrix C

    !**** LOCAL ***********************************!

    !integer ::            M, N, K
    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    real(dp), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:), CC_loc(:,:)
    integer :: A_loc_dim(2), B_loc_dim(2), C_loc_dim(2), coords(2)
    integer :: width, glb_st, glb_ed, loc_st, loc_ed

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
    allocate(B_loc(B_loc_dim(1),B_loc_dim(2)))
    B_loc=0.0_dp
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
          if (width<psp_update_rank) then
             B_loc=0.0_dp
          end if
          call psp_copy_m(width,B_loc_dim(2),B,loc_st,1,B_loc,1,1,1.0_dp,0.0_dp)
       end if
       ! boardcast in column
       call MPI_Bcast(B_loc, B_loc_dim(1)*B_loc_dim(2), MPI_DOUBLE, idx_prow, psp_mpi_comm_col,mpi_err)

       ! compute local update of C
       !C=MATMUL(A_loc,B_loc)+C
       call dgemm('n','n',C_loc_dim(1),C_loc_dim(2),width,1.0_dp,A_loc,A_loc_dim(1),B_loc,psp_update_rank, &
            1.0_dp,C_loc,C_loc_dim(1))
    enddo
    !C=beta*C+alpha*C_loc
    call psp_copy_m(C_loc_dim(1),C_loc_dim(2),C_loc,1,1,C,1,1,alpha,beta)

    if (allocated(A_loc)) deallocate(A_loc)
    if (allocated(B_loc)) deallocate(B_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)

  end subroutine psp_dgemm_nn

  !================================================!
  !     pdgemm for SUMMA  C=alpha*A*B'+beta*C      !
  !================================================!

  subroutine psp_dgemm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, opA(A) of size M by K, opB(B) of size K by N
    real(dp), intent(in)  ::   alpha, beta ! scalar
    real(dp), intent(in) ::   A(:,:) ! matrix A
    real(dp), intent(in) ::   B(:,:) ! matrix B
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: C(:,:) ! matrix C

    !**** LOCAL ***********************************!

    !integer ::            M, N, K
    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    real(dp), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:), CC_loc(:,:)
    integer :: A_loc_dim(2), B_loc_dim(2), C_loc_dim(2), coords(2)
    integer :: width, glb_st, glb_ed, loc_st, loc_ed

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
    allocate(B_loc(B_loc_dim(1),B_loc_dim(2)))
    B_loc=0.0_dp
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
          if (width<psp_update_rank) then
             B_loc=0.0_dp
          endif
          call psp_copy_m(width,B_loc_dim(2),B,loc_st,1,B_loc,1,1,1.0_dp,0.0_dp)
       end if

       ! boardcast in column
       call MPI_Bcast(B_loc, B_loc_dim(1)*B_loc_dim(2), MPI_DOUBLE, idx_prow, psp_mpi_comm_col,mpi_err)

       ! compute local update of C
       ! C_loc = A*(B_loc^t)
       call dgemm(opA,opB,C_loc_dim(1),C_loc_dim(2),A_loc_dim(2),1.0_dp,A,A_loc_dim(1),B_loc,B_loc_dim(1), &
            0.0_dp,C_loc,C_loc_dim(1))

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
    if (allocated(B_loc)) deallocate(B_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)

  end subroutine psp_dgemm_nt

  !================================================!
  !     pdgemm for SUMMA  C=alpha*A'*B+beta*C      !
  !================================================!

  subroutine psp_dgemm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, opA(A) of size M by K, opB(B) of size K by N
    real(dp), intent(in)  ::   alpha, beta ! scalar
    real(dp), intent(in) ::   A(:,:) ! matrix A
    real(dp), intent(in) ::   B(:,:) ! matrix B
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: C(:,:) ! matrix C

    !**** LOCAL ***********************************!

    !integer ::            M, N, K
    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    real(dp), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:), CC_loc(:,:)
    integer :: A_loc_dim(2), B_loc_dim(2), C_loc_dim(2), coords(2)
    integer :: width, glb_st, glb_ed, loc_st, loc_ed

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
       call dgemm(opA,opB,C_loc_dim(1),C_loc_dim(2),B_loc_dim(1),1.0_dp,A_loc,A_loc_dim(1),B,B_loc_dim(1), &
            0.0_dp,C_loc,C_loc_dim(1))

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
    if (allocated(B_loc)) deallocate(B_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)

  end subroutine psp_dgemm_tn

  !================================================!
  !    pdgemm for SUMMA  C=alpha*A'*B'+beta*C      !
  !================================================!

  subroutine psp_dgemm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, opA(A) of size M by K, opB(B) of size K by N
    real(dp), intent(in)  ::   alpha, beta ! scalar
    real(dp), intent(in) ::   A(:,:) ! matrix A
    real(dp), intent(in) ::   B(:,:) ! matrix B
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: C(:,:) ! matrix C

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

    integer, external :: numroc

    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)

    if (M*N<=MIN(M*K,N*K)) then
       ! compute tmp = alpha*B*A
       dims_before(1)=numroc(N,psp_bs_def_row,iprow,0,nprow)
       dims_before(2)=numroc(M,psp_bs_def_col,ipcol,0,npcol)
       call descinit(desc_before,N,M,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_before(1),mpi_err)
       allocate(tmp(dims_before(1),dims_before(2)))
       call psp_dgemm_nn(N,M,K,B,'n',A,'n',tmp,alpha,0.0_dp)

       ! C=transpose(tmp) + beta*C
       dims_after(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
       dims_after(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
       call descinit(desc_after,M,N,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_after(1),mpi_err)
       call pdtran(M,N,1.0_dp,tmp,1,1,desc_before, beta,C,1,1,desc_after)

    else if (M*K<=MIN(M*N,N*K)) then
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
       call psp_dgemm_nt(M,N,K,tmp,'n',B,opB,C,alpha,beta)
    else
       ! tmp = transpose(B)
       dims_before(1)=numroc(N,psp_bs_def_row,iprow,0,nprow)
       dims_before(2)=numroc(K,psp_bs_def_col,ipcol,0,npcol)
       call descinit(desc_before,N,K,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_before(1),mpi_err)
       dims_after(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
       dims_after(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
       call descinit(desc_after,K,N,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_after(1),mpi_err)
       allocate(tmp(dims_after(1),dims_after(2)))
       call pdtran(K,N,1.0_dp,B,1,1,desc_before, 0.0_dp,tmp,1,1,desc_after)

       ! compute C = alpha*transpose(A)*tmp + beta*C
       call psp_dgemm_tn(M,N,K,A,opA,tmp,'n',C,alpha,beta)
    end if
    if (allocated(tmp)) deallocate(tmp)

  end subroutine psp_dgemm_tt

  !================================================!
  !     pzgemm for SUMMA  C=alpha*A*B+beta*C       !
  !================================================!

  subroutine psp_zgemm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, opA(A) of size M by K, opB(B) of size K by N
    complex(dp), intent(in)  ::   alpha, beta ! scalar
    complex(dp), intent(in) ::   A(:,:) ! matrix A
    complex(dp), intent(in) ::   B(:,:) ! matrix B
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: C(:,:) ! matrix C

    !**** LOCAL ***********************************!

    !integer ::            M, N, K
    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    complex(dp), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:), CC_loc(:,:)
    integer :: A_loc_dim(2), B_loc_dim(2), C_loc_dim(2), coords(2)
    integer :: width, glb_st, glb_ed, loc_st, loc_ed

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
    allocate(B_loc(B_loc_dim(1),B_loc_dim(2)))
    B_loc=cmplx_0
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
       call MPI_Bcast(A_loc, A_loc_dim(1)*A_loc_dim(2), MPI_DOUBLE_COMPLEX,  idx_pcol, psp_mpi_comm_row,mpi_err)

       idx_prow = mod(idx_k_row-1,nprow) ! identify the processor owing B(kth block,:), the cart coordinate
       if (iprow==idx_prow) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_row,nprow,loc_st)
          if (width<psp_update_rank) then
             B_loc=cmplx_0
          end if
          call psp_copy_m(width,B_loc_dim(2),B,loc_st,1,B_loc,1,1,cmplx_1,cmplx_0)
       end if
       ! boardcast in column
       call MPI_Bcast(B_loc, B_loc_dim(1)*B_loc_dim(2), MPI_DOUBLE_COMPLEX,  idx_prow, psp_mpi_comm_col,mpi_err)

       ! compute local update of C
       !C=MATMUL(A_loc,B_loc)+C
       call zgemm('n','n',C_loc_dim(1),C_loc_dim(2),width,cmplx_1,A_loc,A_loc_dim(1),B_loc,psp_update_rank, &
            cmplx_1,C_loc,C_loc_dim(1))
    enddo

    !C=beta*C+C_loc
    call psp_copy_m(C_loc_dim(1),C_loc_dim(2),C_loc,1,1,C,1,1,alpha,beta)

    if (allocated(A_loc)) deallocate(A_loc)
    if (allocated(B_loc)) deallocate(B_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)

  end subroutine psp_zgemm_nn

  !================================================!
  !     pzgemm for SUMMA  C=alpha*A*B'+beta*C      !
  !================================================!

  subroutine psp_zgemm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, opA(A) of size M by K, opB(B) of size K by N
    complex(dp), intent(in)  ::   alpha, beta ! scalar
    complex(dp), intent(in) ::   A(:,:) ! matrix A
    complex(dp), intent(in) ::   B(:,:) ! matrix B
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: C(:,:) ! matrix C

    !**** LOCAL ***********************************!

    !integer ::            M, N, K
    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    complex(dp), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:), CC_loc(:,:)
    integer :: A_loc_dim(2), B_loc_dim(2), C_loc_dim(2), coords(2)
    integer :: width, glb_st, glb_ed, loc_st, loc_ed

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
    allocate(B_loc(B_loc_dim(1),B_loc_dim(2)))
    B_loc=cmplx_0
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
          if (width<psp_update_rank) then
             B_loc=cmplx_0
          endif
          call psp_copy_m(width,B_loc_dim(2),B,loc_st,1,B_loc,1,1,cmplx_1,cmplx_0)
       end if

       ! boardcast in column
       call MPI_Bcast(B_loc, B_loc_dim(1)*B_loc_dim(2), MPI_DOUBLE_COMPLEX,  idx_prow, psp_mpi_comm_col,mpi_err)

       ! compute local update of C
       ! C_loc = A*(B_loc^t)
       call zgemm(opA,opB,C_loc_dim(1),C_loc_dim(2),A_loc_dim(2),cmplx_1,A,A_loc_dim(1),B_loc,B_loc_dim(1), &
            cmplx_0,C_loc,C_loc_dim(1))

       idx_pcol = mod(idx_k_col-1,npcol ) ! identify the processor owing C(:,kth block,), the cart coordinate
       ! boardcast in row
       !CC_loc=0.0_dp ! not necessary
       call MPI_REDUCE(C_loc, CC_loc, C_loc_dim(1)*C_loc_dim(2), MPI_DOUBLE_COMPLEX,  MPI_SUM, idx_pcol, psp_mpi_comm_row, mpi_err)

       if (ipcol==idx_pcol) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_col,npcol,loc_st)
          !C=beta*C+C_loc
          call psp_copy_m(C_loc_dim(1),width,CC_loc,1,1,C,1,loc_st,alpha,beta)
       end if
    enddo

    if (allocated(A_loc)) deallocate(A_loc)
    if (allocated(B_loc)) deallocate(B_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)

  end subroutine psp_zgemm_nt

  !================================================!
  !     pzgemm for SUMMA  C=alpha*A'*B+beta*C      !
  !================================================!

  subroutine psp_zgemm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, opA(A) of size M by K, opB(B) of size K by N
    complex(dp), intent(in)  ::   alpha, beta ! scalar
    complex(dp), intent(in) ::   A(:,:) ! matrix A
    complex(dp), intent(in) ::   B(:,:) ! matrix B
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: C(:,:) ! matrix C

    !**** LOCAL ***********************************!

    !integer ::            M, N, K
    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2, chunk_sz
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    complex(dp), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:), CC_loc(:,:)
    integer :: A_loc_dim(2), B_loc_dim(2), C_loc_dim(2), coords(2)
    integer :: width, glb_st, glb_ed, loc_st, loc_ed

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
          !call psp_idx_glb2loc(glb_ed,psp_bs_def_col,npcol,loc_ed)
          if (width<psp_update_rank) then
             A_loc=cmplx_0
          endif
          call psp_copy_m(A_loc_dim(1),width,A,1,loc_st,A_loc,1,1,cmplx_1,cmplx_0)
       end if

       ! boardcast in row
       call MPI_Bcast(A_loc, A_loc_dim(1)*A_loc_dim(2), MPI_DOUBLE_COMPLEX,  idx_pcol, psp_mpi_comm_row,mpi_err)

       ! compute local update of C
       ! C_loc = (A_loc)^t*B
       call zgemm(opA,opB,C_loc_dim(1),C_loc_dim(2),B_loc_dim(1),cmplx_1,A_loc,A_loc_dim(1),B,B_loc_dim(1), &
            cmplx_0,C_loc,C_loc_dim(1))

       idx_prow = mod(idx_k_row-1,nprow ) ! identify the processor owing C(kth block,:), the cart coordinate
       ! boardcast in column
       call MPI_REDUCE(C_loc, CC_loc, C_loc_dim(1)*C_loc_dim(2), MPI_DOUBLE_COMPLEX,  MPI_SUM, idx_prow, psp_mpi_comm_col, mpi_err)

       if (iprow==idx_prow) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_row,nprow,loc_st)
          !C=beta*C+C_loc
          call psp_copy_m(width,C_loc_dim(2),CC_loc,1,1,C,loc_st,1,alpha,beta)
       end if
    enddo

    if (allocated(A_loc)) deallocate(A_loc)
    if (allocated(B_loc)) deallocate(B_loc)
    if (allocated(C_loc)) deallocate(C_loc)
    if (allocated(CC_loc)) deallocate(CC_loc)

  end subroutine psp_zgemm_tn

  !================================================!
  !    pzgemm for SUMMA  C=alpha*A'*B'+beta*C      !
  !================================================!

  subroutine psp_zgemm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, opA(A) of size M by K, opB(B) of size K by N
    complex(dp), intent(in)  ::   alpha, beta ! scalar
    complex(dp), intent(in) ::   A(:,:) ! matrix A
    complex(dp), intent(in) ::   B(:,:) ! matrix B
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    complex(dp), intent(inout) ::   C(:,:) ! matrix C

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

    integer, external :: numroc

    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)

    if (M*N<=MIN(M*K,N*K)) then
       ! compute tmp = alpha*B*A
       call psp_process_opM(opA,trA)
       call psp_process_opM(opB,trB)
       dims_before(1)=numroc(N,psp_bs_def_row,iprow,0,nprow)
       dims_before(2)=numroc(M,psp_bs_def_col,ipcol,0,npcol)
       call descinit(desc_before,N,M,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_before(1),mpi_err)
       allocate(tmp(dims_before(1),dims_before(2)))
       if (trA==1 .and. trB==2) then
          dims_conj(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
          dims_conj(2)=numroc(M,psp_bs_def_col,ipcol,0,npcol)
          allocate(tmpConj(dims_conj(1),dims_conj(2)))
          tmpConj=CONJG(A)
          call psp_zgemm_nn(N,M,K,B,'n',tmpConj,'n',tmp,alpha,cmplx_0)
       else if (tra==2 .and. trB==1) then
          dims_conj(1)=numroc(N,psp_bs_def_row,iprow,0,nprow)
          dims_conj(2)=numroc(K,psp_bs_def_col,ipcol,0,npcol)
          allocate(tmpConj(dims_conj(1),dims_conj(2)))
          tmpConj=CONJG(B)
          call psp_zgemm_nn(N,M,K,tmpConj,'n',A,'n',tmp,alpha,cmplx_0)
       else
          call psp_zgemm_nn(N,M,K,B,'n',A,'n',tmp,alpha,cmplx_0)
       end if

       ! C=transpose(tmp) + beta*C
       dims_after(1)=numroc(M,psp_bs_def_row,iprow,0,nprow)
       dims_after(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
       call descinit(desc_after,M,N,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_after(1),mpi_err)
       if (trA==1 .and. trB==1) then
          call pztranc(M,N,cmplx_1,tmp,1,1,desc_before, beta,C,1,1,desc_after)
       else
          call pztranu(M,N,cmplx_1,tmp,1,1,desc_before, beta,C,1,1,desc_after)
       end if
    else if (M*K<=MIN(M*N,N*K)) then
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
       call psp_zgemm_nt(M,N,K,tmp,'n',B,opB,C,alpha,beta)
    else
       ! tmp = transpose(B)
       dims_before(1)=numroc(N,psp_bs_def_row,iprow,0,nprow)
       dims_before(2)=numroc(K,psp_bs_def_col,ipcol,0,npcol)
       call descinit(desc_before,N,K,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_before(1),mpi_err)
       dims_after(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
       dims_after(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
       call descinit(desc_after,K,N,psp_bs_def_row,psp_bs_def_col,0,0,psp_icontxt,dims_after(1),mpi_err)
       allocate(tmp(dims_after(1),dims_after(2)))
       call psp_process_opM(opB,trB)
       if (trB==1) then
          call pztranc(K,N,cmplx_1,B,1,1,desc_before,cmplx_0,tmp,1,1,desc_after)
       else
          call pztranu(K,N,cmplx_1,B,1,1,desc_before,cmplx_0,tmp,1,1,desc_after)
       end if

       ! compute C = alpha*transpose(A)*tmp + beta*C
       call psp_zgemm_tn(M,N,K,A,opA,tmp,'n',C,alpha,beta)
    end if

    if (allocated(tmp)) deallocate(tmp)
    if (allocated(tmpConj)) deallocate(tmpConj)

  end subroutine psp_zgemm_tt


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


END MODULE pspGemm
