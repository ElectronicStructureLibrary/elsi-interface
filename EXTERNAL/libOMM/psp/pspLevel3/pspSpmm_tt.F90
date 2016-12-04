MODULE pspSpmm_tt
  use pspVariable
  use pspUtility
  use pspMPI
  use pspLevel1
  use pspLevel2
  use pspSpmm_nt, only: psp_gespmm_nt
  use pspSpmm_tn, only: psp_gespmm_tn

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

  interface psp_gespmm_tt
     module procedure psp_dgespmm_tt
     module procedure psp_zgespmm_tt
  end interface psp_gespmm_tt

  interface die
     module procedure die
  end interface die

  public :: psp_gespmm_tt

contains

  !================================================!
  !   sparse pdgemm: spmm  C=alpha*A'*B'+beta*C    !
  !================================================!

  subroutine psp_dgespmm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, A of size M by K, B of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    real(dp), intent(in)  ::   alpha, beta ! scalar
    type(psp_matrix_spm), intent(in) :: A ! matrix A
    real(dp), intent(in) :: B(:,:)

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

    ! TODO: optimize the following code by introducing a transpose operator for sparse matrices
    if (.true.) then
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
       call psp_gespmm_tn(M,N,K,A,opA,tmp,'n',C,alpha,beta)
    end if

    if (allocated(tmp)) deallocate(tmp)

  end subroutine psp_dgespmm_tt

  !================================================!
  !   sparse pzgemm: spmm  C=alpha*A'*B'+beta*C    !
  !================================================!

  subroutine psp_zgespmm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, A of size M by K, B of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    complex(dp), intent(in)  ::   alpha, beta ! scalar
    type(psp_matrix_spm), intent(in) :: A ! matrix A
    complex(dp), intent(in) :: B(:,:)

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

    ! TODO: optimize the following code by introducing a transpose operator for sparse matrices
    if (.true.) then
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
       call psp_gespmm_tn(M,N,K,A,opA,tmp,'n',C,alpha,beta)
    end if

    if (allocated(tmp)) deallocate(tmp)

  end subroutine psp_zgespmm_tt


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


END MODULE pspSpmm_tt
