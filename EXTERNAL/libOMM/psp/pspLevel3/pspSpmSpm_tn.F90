MODULE pspSpmSpm_tn
  use pspVariable
  use pspUtility
  use pspMPI
  use pspLevel1
  use pspLevel2
  use pspMatSum

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

  interface psp_gespmspm_tn
     module procedure psp_dgespmspm_tn
     module procedure psp_zgespmspm_tn
  end interface psp_gespmspm_tn

  interface die
     module procedure die
  end interface die

  public :: psp_gespmspm_tn

contains

  !================================================!
  !    sparse pdgemm: spmspm  C=alpha*A'*B+beta*C  !
  !================================================!

  subroutine psp_dgespmspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, A of size M by K, B of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    real(dp), intent(in)  ::   alpha, beta ! scalar
    type(psp_matrix_spm), intent(in) :: A, B ! matrix A and B

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout) :: C

    !**** LOCAL ***********************************!

    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    real(dp), allocatable :: A_loc_val(:)
    integer, allocatable :: A_loc_idx1(:), A_loc_idx2(:)
    integer :: A_loc_dim(2), C_loc_dim(2), coords(2), disp(2)
    integer :: width, glb_st, glb_ed, loc_st, loc_ed, nnz_loc, L, C_loc_total_nnz
    type(psp_MPI_dspm) :: C_loc
    type(dList), pointer :: list
    logical, allocatable :: listCreated(:), listCreatedLoc(:)
    type(dListPtrArray), allocatable :: listArraySt(:), listArrayEd(:)
    type(dListPtrArray), allocatable :: listArrayStLoc(:), listArrayEdLoc(:)

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
    A_loc_dim(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
    A_loc_dim(2)=psp_update_rank
    C_loc_dim(1)=psp_update_rank
    C_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    C_loc_total_nnz=0
    allocate(listArraySt(C_loc_dim(2)))
    allocate(listArrayEd(C_loc_dim(2)))
    allocate(listArrayStLoc(C_loc_dim(2)))
    allocate(listArrayEdLoc(C_loc_dim(2)))
    allocate(listCreated(C_loc_dim(2)))
    do cnt1=1,C_loc_dim(2)
       listCreated(cnt1)=.false.
    end do
    allocate(listCreatedLoc(C_loc_dim(2)))

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
          ! copy the sparse local submatrix
          ! A_loc=A_loc(1:A_loc_dim(1),1:width)=A(1:A_loc_dim(1),loc_st:loc_st+width-1)
          call psp_copy_spm2st(A_loc_dim(1),width,A,1,loc_st,A_loc_idx1,A_loc_idx2, &
               A_loc_val,A_loc_dim(1),width,1,1,0.0_dp)
       else
          if (allocated(A_loc_idx1)) deallocate(A_loc_idx1)
          if (allocated(A_loc_idx2)) deallocate(A_loc_idx2)
          if (allocated(A_loc_val)) deallocate(A_loc_val)
          allocate(A_loc_idx2(width+1))
          allocate(A_loc_idx1(1))
          allocate(A_loc_val(1))
       end if

       ! boardcast in row
       if (psp_update_rank>1) then
          call MPI_Bcast(width, 1, MPI_INT, idx_pcol, psp_mpi_comm_row,mpi_err)
       end if
       L=size(A_loc_idx2)
       ! data to be boardcast have the same size in every processor
       if (L/=width+1) then
          if (allocated(A_loc_idx2)) deallocate(A_loc_idx2)
          allocate(A_loc_idx2(width+1))
       end if
       call MPI_Bcast(A_loc_idx2, width+1, MPI_INT, idx_pcol, psp_mpi_comm_row,mpi_err)
       nnz_loc=A_loc_idx2(width+1)-A_loc_idx2(1)
       L=size(A_loc_idx1)
       if (L/=nnz_loc) then
          if (allocated(A_loc_idx1)) deallocate(A_loc_idx1)
          if (allocated(A_loc_val)) deallocate(A_loc_val)
          allocate(A_loc_idx1(nnz_loc))
          allocate(A_loc_val(nnz_loc))
       end if
       call MPI_Bcast(A_loc_idx1, nnz_loc, MPI_INT, idx_pcol, psp_mpi_comm_row,mpi_err)
       call MPI_Bcast(A_loc_val, nnz_loc, MPI_DOUBLE, idx_pcol, psp_mpi_comm_row,mpi_err)

       ! compute local update of C
       ! C_loc = (A_loc)^t*B, where A_loc and B are sparse matrices
       call psp_sst_gespmspm(width,C_loc_dim(2),A_loc_dim(1), opA,opB,1.0_dp,A_loc_idx1,A_loc_idx2,A_loc_val, &
            B%row_ind,B%col_ptr,B%dval,C_loc%idx1,C_loc%idx2,C_loc%val,0.0_dp)

       idx_prow = mod(idx_k_row-1,nprow ) ! identify the processor owing C(kth block,:), the cart coordinate
       C_loc%nnz=size(C_loc%val)
       C_loc%loc_dim1=width
       C_loc%loc_dim2=C_loc_dim(2)
       ! reduce in column
       call psp_MPI_REDUCE_spm_packed(C_loc, idx_prow, .false.)

       if (iprow==idx_prow) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_row,nprow,loc_st)
          ! stack C_loc together into a larger C_loc=[C_loc_1;C_loc_2;...;C_loc_(M/psp_update_rank)] using linked list
          if (C_loc%nnz>0) then
             disp(1)=loc_st-1
             disp(2)=0
             C_loc_total_nnz=C_loc_total_nnz+C_loc%nnz
             call psp_spm2lists_shift(width,C_loc_dim(2),C_loc%idx1,C_loc%idx2,C_loc%val,&
                  C_loc%nnz,'csc',listArrayStLoc,listArrayEdLoc,listCreatedLoc,disp)
             do cnt1=1,C_Loc_dim(2)
                if (listCreated(cnt1)) then
                   if (listCreatedLoc(cnt1)) then
                      listArrayEd(cnt1)%ptr%next=>listArrayStLoc(cnt1)%ptr
                      listArrayEd(cnt1)%ptr=>listArrayEdLoc(cnt1)%ptr
                   end if
                else
                   if (listCreatedLoc(cnt1)) then
                      ! create the cnt'th list
                      listArraySt(cnt1)%ptr=>listArrayStLoc(cnt1)%ptr
                      listArrayEd(cnt1)%ptr=>listArrayEdLoc(cnt1)%ptr
                      listCreated(cnt1)=.true.
                   end if
                end if
             end do
          end if
       end if
    enddo

    if (C_loc_total_nnz>0) then
       ! create list
       cnt1=1
       do while (cnt1<=C_loc_dim(2))
          if (listCreated(cnt1)) then
             list=>listArraySt(cnt1)%ptr
             exit
          else
             cnt1=cnt1+1
          end if
       end do
       ! link lists for all columns
       if (C_loc_dim(2)>1) then
          do cnt1=1,C_loc_dim(2)-1
             if (listCreated(cnt1)) then
                if (listCreated(cnt1+1)) then
                   listArrayEd(cnt1)%ptr%next=>listArraySt(cnt1+1)%ptr
                else
                   listArrayEd(cnt1+1)%ptr=>listArrayEd(cnt1)%ptr
                   listArraySt(cnt1+1)%ptr=>listArraySt(cnt1)%ptr
                end if
             else
                if (listCreated(cnt1+1)) then
                   listArrayEd(cnt1)%ptr%next=>listArraySt(cnt1+1)%ptr
                else
                   listArrayEd(cnt1+1)%ptr=>listArrayEd(cnt1)%ptr
                   listArraySt(cnt1+1)%ptr=>listArraySt(cnt1)%ptr
                end if
             end if
          end do
       end if
    end if

    ! list to spm
    call psp_list2spm(C%loc_dim1,C%loc_dim2,C_loc%idx1,C_loc%idx2,C_loc%val,&
         'csc',list,C_loc_total_nnz,.false.)
    ! compute C=alpha*C_loc+beta*C
    call psp_sst_sum_spmspm(C%loc_dim1,C%loc_dim2,alpha,C_loc%idx1,C_loc%idx2,C_loc%val,'csc',beta,&
         C%row_ind,C%col_ptr,C%dval,'csc',C%row_ind,C%col_ptr,C%dval,'csc',C_loc_total_nnz,C%nnz)


    if (allocated(C_loc%idx1)) deallocate(C_loc%idx1)
    if (allocated(C_loc%idx2)) deallocate(C_loc%idx2)
    if (allocated(C_loc%val)) deallocate(C_loc%val)
    if (allocated(A_loc_val)) deallocate(A_loc_val)
    if (allocated(A_loc_idx1)) deallocate(A_loc_idx1)
    if (allocated(A_loc_idx2)) deallocate(A_loc_idx2)
    if (allocated(listArraySt)) deallocate(listArraySt)
    if (allocated(listArrayEd)) deallocate(listArrayEd)
    if (allocated(listArrayStLoc)) deallocate(listArrayStLoc)
    if (allocated(listArrayEdLoc)) deallocate(listArrayEdLoc)
    if (allocated(listCreated)) deallocate(listCreated)
    if (allocated(listCreatedLoc)) deallocate(listCreatedLoc)
    call list_destroy(list)

  end subroutine psp_dgespmspm_tn


  !================================================!
  !    sparse pzgemm: spmspm  C=alpha*A'*B+beta*C  !
  !================================================!


  subroutine psp_zgespmspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, A of size M by K, B of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    complex(dp), intent(in)  ::   alpha, beta ! scalar
    type(psp_matrix_spm), intent(in) :: A, B ! matrix A and B

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout) :: C

    !**** LOCAL ***********************************!

    integer :: kloop, idx_k_row, idx_k_col, idx_pcol, idx_prow, cnt1, cnt2
    integer :: iprow, ipcol, nprow, npcol, mpi_err
    complex(dp), allocatable :: A_loc_val(:)
    integer, allocatable :: A_loc_idx1(:), A_loc_idx2(:)
    integer :: A_loc_dim(2), C_loc_dim(2), coords(2), disp(2)
    integer :: width, glb_st, glb_ed, loc_st, loc_ed, nnz_loc, L, C_loc_total_nnz
    type(psp_MPI_zspm) :: C_loc
    type(zList), pointer :: list
    logical, allocatable :: listCreated(:), listCreatedLoc(:)
    type(zListPtrArray), allocatable :: listArraySt(:), listArrayEd(:)
    type(zListPtrArray), allocatable :: listArrayStLoc(:), listArrayEdLoc(:)

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
    A_loc_dim(1)=numroc(K,psp_bs_def_row,iprow,0,nprow)
    A_loc_dim(2)=psp_update_rank
    C_loc_dim(1)=psp_update_rank
    C_loc_dim(2)=numroc(N,psp_bs_def_col,ipcol,0,npcol)
    C_loc_total_nnz=0
    allocate(listArraySt(C_loc_dim(2)))
    allocate(listArrayEd(C_loc_dim(2)))
    allocate(listArrayStLoc(C_loc_dim(2)))
    allocate(listArrayEdLoc(C_loc_dim(2)))
    allocate(listCreated(C_loc_dim(2)))
    do cnt1=1,C_loc_dim(2)
       listCreated(cnt1)=.false.
    end do
    allocate(listCreatedLoc(C_loc_dim(2)))

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
          ! copy the sparse local submatrix
          ! A_loc=A_loc(1:A_loc_dim(1),1:width)=A(1:A_loc_dim(1),loc_st:loc_st+width-1)
          call psp_copy_spm2st(A_loc_dim(1),width,A,1,loc_st,A_loc_idx1,A_loc_idx2, &
               A_loc_val,A_loc_dim(1),width,1,1,cmplx_0)
       else
          if (allocated(A_loc_idx1)) deallocate(A_loc_idx1)
          if (allocated(A_loc_idx2)) deallocate(A_loc_idx2)
          if (allocated(A_loc_val)) deallocate(A_loc_val)
          allocate(A_loc_idx2(width+1))
          allocate(A_loc_idx1(1))
          allocate(A_loc_val(1))
       end if

       ! boardcast in row
       if (psp_update_rank>1) then
          call MPI_Bcast(width, 1, MPI_INT, idx_pcol, psp_mpi_comm_row,mpi_err)
       end if
       L=size(A_loc_idx2)
       ! data to be boardcast have the same size in every processor
       if (L/=width+1) then
          if (allocated(A_loc_idx2)) deallocate(A_loc_idx2)
          allocate(A_loc_idx2(width+1))
       end if
       call MPI_Bcast(A_loc_idx2, width+1, MPI_INT, idx_pcol, psp_mpi_comm_row,mpi_err)
       nnz_loc=A_loc_idx2(width+1)-A_loc_idx2(1)
       L=size(A_loc_idx1)
       if (L/=nnz_loc) then
          if (allocated(A_loc_idx1)) deallocate(A_loc_idx1)
          if (allocated(A_loc_val)) deallocate(A_loc_val)
          allocate(A_loc_idx1(nnz_loc))
          allocate(A_loc_val(nnz_loc))
       end if
       call MPI_Bcast(A_loc_idx1, nnz_loc, MPI_INT, idx_pcol, psp_mpi_comm_row,mpi_err)
       call MPI_Bcast(A_loc_val, nnz_loc, MPI_DOUBLE_COMPLEX, idx_pcol, psp_mpi_comm_row,mpi_err)

       ! compute local update of C
       ! C_loc = (A_loc)^t*B, where A_loc and B are sparse matrices
       call psp_sst_gespmspm(width,C_loc_dim(2),A_loc_dim(1), opA,opB,&
            cmplx_1,A_loc_idx1,A_loc_idx2,A_loc_val, &
            B%row_ind,B%col_ptr,B%zval,C_loc%idx1,C_loc%idx2,C_loc%val,cmplx_0)

       idx_prow = mod(idx_k_row-1,nprow ) ! identify the processor owing C(kth block,:), the cart coordinate
       C_loc%nnz=size(C_loc%val)
       C_loc%loc_dim1=width
       C_loc%loc_dim2=C_loc_dim(2)
       ! reduce in column
       call psp_MPI_REDUCE_spm_packed(C_loc, idx_prow, .false.)

       if (iprow==idx_prow) then
          call psp_idx_glb2loc(glb_st,psp_bs_def_row,nprow,loc_st)
          ! stack C_loc together into a larger C_loc=[C_loc_1;C_loc_2;...;C_loc_(M/psp_update_rank)] using linked list
          if (C_loc%nnz>0) then
             disp(1)=loc_st-1
             disp(2)=0
             C_loc_total_nnz=C_loc_total_nnz+C_loc%nnz
             call psp_spm2lists_shift(width,C_loc_dim(2),C_loc%idx1,C_loc%idx2,C_loc%val,&
                  C_loc%nnz,'csc',listArrayStLoc,listArrayEdLoc,listCreatedLoc,disp)
             do cnt1=1,C_Loc_dim(2)
                if (listCreated(cnt1)) then
                   if (listCreatedLoc(cnt1)) then
                      listArrayEd(cnt1)%ptr%next=>listArrayStLoc(cnt1)%ptr
                      listArrayEd(cnt1)%ptr=>listArrayEdLoc(cnt1)%ptr
                   end if
                else
                   if (listCreatedLoc(cnt1)) then
                      ! create the cnt'th list
                      listArraySt(cnt1)%ptr=>listArrayStLoc(cnt1)%ptr
                      listArrayEd(cnt1)%ptr=>listArrayEdLoc(cnt1)%ptr
                      listCreated(cnt1)=.true.
                   end if
                end if
             end do
          end if
       end if
    enddo

    if (C_loc_total_nnz>0) then
       ! create list
       cnt1=1
       do while (cnt1<=C_loc_dim(2))
          if (listCreated(cnt1)) then
             list=>listArraySt(cnt1)%ptr
             exit
          else
             cnt1=cnt1+1
          end if
       end do
       ! link lists for all columns
       if (C_loc_dim(2)>1) then
          do cnt1=1,C_loc_dim(2)-1
             if (listCreated(cnt1)) then
                if (listCreated(cnt1+1)) then
                   listArrayEd(cnt1)%ptr%next=>listArraySt(cnt1+1)%ptr
                else
                   listArrayEd(cnt1+1)%ptr=>listArrayEd(cnt1)%ptr
                   listArraySt(cnt1+1)%ptr=>listArraySt(cnt1)%ptr
                end if
             else
                if (listCreated(cnt1+1)) then
                   listArrayEd(cnt1)%ptr%next=>listArraySt(cnt1+1)%ptr
                else
                   listArrayEd(cnt1+1)%ptr=>listArrayEd(cnt1)%ptr
                   listArraySt(cnt1+1)%ptr=>listArraySt(cnt1)%ptr
                end if
             end if
          end do
       end if
    end if

    ! list to spm
    call psp_list2spm(C%loc_dim1,C%loc_dim2,C_loc%idx1,C_loc%idx2,C_loc%val,&
         'csc',list,C_loc_total_nnz,.false.)
    ! compute C=alpha*C_loc+beta*C
    call psp_sst_sum_spmspm(C%loc_dim1,C%loc_dim2,alpha,C_loc%idx1,C_loc%idx2,C_loc%val,'csc',beta,&
         C%row_ind,C%col_ptr,C%zval,'csc',C%row_ind,C%col_ptr,C%zval,'csc',C_loc_total_nnz,C%nnz)


    if (allocated(C_loc%idx1)) deallocate(C_loc%idx1)
    if (allocated(C_loc%idx2)) deallocate(C_loc%idx2)
    if (allocated(C_loc%val)) deallocate(C_loc%val)
    if (allocated(A_loc_val)) deallocate(A_loc_val)
    if (allocated(A_loc_idx1)) deallocate(A_loc_idx1)
    if (allocated(A_loc_idx2)) deallocate(A_loc_idx2)
    if (allocated(listArraySt)) deallocate(listArraySt)
    if (allocated(listArrayEd)) deallocate(listArrayEd)
    if (allocated(listArrayStLoc)) deallocate(listArrayStLoc)
    if (allocated(listArrayEdLoc)) deallocate(listArrayEdLoc)
    if (allocated(listCreated)) deallocate(listCreated)
    if (allocated(listCreatedLoc)) deallocate(listCreatedLoc)
    call list_destroy(list)

  end subroutine psp_zgespmspm_tn

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



END MODULE pspSpmSpm_tn
