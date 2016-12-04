MODULE psp_spBLAS_Level3
  use pspVariable
  use pspListTool
  use pspBasicTool
  use psp_spBLAS_Level1
  use psp_spBLAS_Level2

  ! This module contains sequential sparse BLAS

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

  interface psp_sst_sum_spmm
     module procedure psp_sst_dsum_spmm
     module procedure psp_sst_zsum_spmm
  end interface psp_sst_sum_spmm

  interface psp_sst_sum_spmspm
     module procedure psp_sst_dsum_spmspm
     module procedure psp_sst_zsum_spmspm
  end interface psp_sst_sum_spmspm

  interface psp_sst_gespmm
     module procedure psp_sst_dgespmm
     module procedure psp_sst_zgespmm
  end interface psp_sst_gespmm

  interface psp_sst_gemspm
     module procedure psp_sst_dgemspm
     module procedure psp_sst_zgemspm
  end interface psp_sst_gemspm

  interface psp_sst_gespmspm
     module procedure psp_sst_dgespmspm
     module procedure psp_sst_zgespmspm
  end interface psp_sst_gespmspm

  interface die
     module procedure die
  end interface die

  public :: psp_sst_sum_spmm
  public :: psp_sst_sum_spmspm
  public :: psp_sst_gespmm
  public :: psp_sst_gemspm
  public :: psp_sst_gespmspm

contains

  subroutine psp_sst_dsum_spmm(m,n,alpha,idx1,idx2,val,fmt,beta,B,C)
    ! C = alpha*A+beta*B

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt ! storage format of sparse matrices, 'coo' or 'csc'
    integer, intent(in)  :: m, n
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(in) :: val(:)
    integer, intent(in) :: idx1(:)
    integer, intent(in) :: idx2(:)
    real(dp), intent(in) ::   B(:,:) ! dense matrix

    !**** INOUT ***********************************!
    real(dp), intent(inout) ::   C(:,:) ! dense matrix

    !**** INTERNAL ********************************!

    integer :: cnt, cnt2, nnz

    !**********************************************!
    do cnt2=1,n
       do cnt=1,m
          C(cnt,cnt2) = beta*B(cnt,cnt2)
       end do
    end do
    nnz=size(val)
    select case (fmt)
    case('coo')
       do cnt=1,nnz
          C(idx1(cnt),idx2(cnt))=C(idx1(cnt),idx2(cnt))+alpha*val(cnt)
       end do
    case('csc')
       do cnt=1,n
          do cnt2=idx2(cnt),idx2(cnt+1)-1
             C(idx1(cnt2),cnt)=C(idx1(cnt2),cnt)+alpha*val(cnt2)
          end do
       end do
    end select

  end subroutine psp_sst_dsum_spmm

  subroutine psp_sst_zsum_spmm(m,n,alpha,idx1,idx2,val,fmt,beta,B,C)
    ! C = alpha*A+beta*B

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt ! storage format of sparse matrices, 'coo' or 'csc'
    integer, intent(in)  :: m, n
    complex(dp), intent(in) :: alpha, beta
    complex(dp), intent(in) :: val(:)
    integer, intent(in) :: idx1(:)
    integer, intent(in) :: idx2(:)
    complex(dp), intent(in) ::   B(:,:) ! dense matrix

    !**** INOUT ***********************************!
    complex(dp), intent(inout) ::   C(:,:) ! dense matrix

    !**** INTERNAL ********************************!

    integer :: cnt, cnt2, nnz

    !**********************************************!
    do cnt2=1,n
       do cnt=1,m
          C(cnt,cnt2) = beta*B(cnt,cnt2)
       end do
    end do
    nnz=size(val)
    select case (fmt)
    case('coo')
       do cnt=1,nnz
          C(idx1(cnt),idx2(cnt))=C(idx1(cnt),idx2(cnt))+alpha*val(cnt)
       end do
    case('csc')
       do cnt=1,n
          do cnt2=idx2(cnt),idx2(cnt+1)-1
             C(idx1(cnt2),cnt)=C(idx1(cnt2),cnt)+alpha*val(cnt2)
          end do
       end do
    end select

  end subroutine psp_sst_zsum_spmm

  subroutine psp_sst_dsum_spmspm(m,n,alpha,idx1A,idx2A,valA,fmtA,beta,&
       idx1B,idx2B,valB,fmtB,idx1C,idx2C,valC,fmtC,nnzA,nnzB)
    ! C = alpha*A+beta*B

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmtA, fmtB, fmtC ! storage format of sparse matrices, 'coo' or 'csc'
    integer, intent(in)  :: m, n
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(in) :: valA(:), valB(:)
    integer, intent(in) :: idx1A(:), idx1B(:)
    integer, intent(in) :: idx2A(:), idx2B(:)
    integer, intent(in) :: nnzA, nnzB

    !**** INOUT ***********************************!
    real(dp), allocatable, intent(inout) :: valC(:)
    integer, allocatable, intent(inout) :: idx1C(:)
    integer, allocatable, intent(inout) :: idx2C(:)

    !**** INTERNAL ********************************!
    type(dList), pointer :: list, lastElem
    integer :: cnt, cnt2, nnz
    real(dp), allocatable :: valCLoc(:)
    integer, allocatable :: idx1CLoc(:)
    integer, allocatable :: idx2CLoc(:)
    !**********************************************!

    ! let CLoc=alpha*A
    if (fmtA=='csc') then
       nnz=n+1
       allocate(idx2CLoc(nnz))
       do cnt=1,nnz
          idx2CLoc(cnt)=idx2A(cnt)
       end do
    else
       nnz=nnzA
       allocate(idx2CLoc(nnz))
       do cnt=1,nnz
          idx2CLoc(cnt)=idx2A(cnt)
       end do
    end if
    nnz=nnzA!size(idx1A)
    allocate(idx1CLoc(nnz))
    allocate(valCLoc(nnz))
    do cnt=1,nnz
       idx1CLoc(cnt)=idx1A(cnt)
       valCLoc(cnt)=alpha*valA(cnt)
    end do

    ! convert the format of CLoc to fmtC
    call psp_sst_fmtCnvt(m,n,nnz,idx2CLoc,fmtA,fmtC)
    ! create a list storage for CLoc
    call psp_spm2list(m,n,idx1CLoc,idx2CLoc,valCLoc,nnz,fmtC,list)
    ! add B to list and update nnz, essentially computing CLoc=CLoc+beta*B in a form of lists
    if (nnz==0) then
       call psp_list_create_mat(list,nnz,beta,idx1B,idx2B,valB,fmtB,nnzB,lastElem)
    else
       call psp_list_combine_listMat(m,n,1.0_dp,list,nnz,beta,idx1B,idx2B,valB,fmtB,nnzB,lastElem)
    end if
    ! convert format list to a sparse matrix C with nnz nonzero elements
    if (allocated(idx2C)) deallocate(idx2C)
    if (allocated(idx1C)) deallocate(idx1C)
    if (allocated(valC)) deallocate(valC)
    if (fmtC=='csc') then
       allocate(idx2C(n+1))
    else
       allocate(idx2C(nnz))
    end if
    allocate(idx1C(nnz))
    allocate(valC(nnz))
    call psp_list2spm(m,n,idx1C,idx2C,valC,fmtC,list,nnz,.false.)

    deallocate(valCLoc)
    deallocate(idx1CLoc)
    deallocate(idx2CLoc)
    if (nnz>0) call list_destroy(list)

  end subroutine psp_sst_dsum_spmspm

  subroutine psp_sst_zsum_spmspm(m,n,alpha,idx1A,idx2A,valA,fmtA,beta,&
       idx1B,idx2B,valB,fmtB,idx1C,idx2C,valC,fmtC,nnzA,nnzB)
    ! C = alpha*A+beta*B

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmtA, fmtB, fmtC ! storage format of sparse matrices, 'coo' or 'csc'
    integer, intent(in)  :: m, n
    complex(dp), intent(in) :: alpha, beta
    complex(dp), intent(in) :: valA(:), valB(:)
    integer, intent(in) :: idx1A(:), idx1B(:)
    integer, intent(in) :: idx2A(:), idx2B(:)
    integer, intent(in) :: nnzA, nnzB

    !**** INOUT ***********************************!
    complex(dp), allocatable, intent(inout) :: valC(:)
    integer, allocatable, intent(inout) :: idx1C(:)
    integer, allocatable, intent(inout) :: idx2C(:)

    !**** INTERNAL ********************************!
    type(zList), pointer :: list, lastElem
    integer :: cnt, cnt2, nnz
    complex(dp), allocatable :: valCLoc(:)
    integer, allocatable :: idx1CLoc(:)
    integer, allocatable :: idx2CLoc(:)
    !**********************************************!

    ! let CLoc=alpha*A
    if (fmtA=='csc') then
       nnz=n+1
       allocate(idx2CLoc(nnz))
       do cnt=1,nnz
          idx2CLoc(cnt)=idx2A(cnt)
       end do
    else
       nnz=nnzA
       allocate(idx2CLoc(nnz))
       do cnt=1,nnz
          idx2CLoc(cnt)=idx2A(cnt)
       end do
    end if
    nnz=nnzA!size(idx1A)
    allocate(idx1CLoc(nnz))
    allocate(valCLoc(nnz))
    do cnt=1,nnz
       idx1CLoc(cnt)=idx1A(cnt)
       valCLoc(cnt)=alpha*valA(cnt)
    end do

    ! convert the format of CLoc to fmtC
    call psp_sst_fmtCnvt(m,n,nnz,idx2CLoc,fmtA,fmtC)
    ! create list
    call psp_spm2list(m,n,idx1CLoc,idx2CLoc,valCLoc,nnz,fmtC,list)
    ! add B to list and update nnz
    if (nnz==0) then
       call psp_list_create_mat(list,nnz,beta,idx1B,idx2B,valB,fmtB,size(valB),lastElem)
    else
       call psp_list_combine_listMat(m,n,cmplx_1,list,nnz,beta,idx1B,idx2B,valB,fmtB,size(valB),lastElem)
    end if

    ! convert format list to a sparse matrix with nnz nonzero elements
    if (allocated(idx2C)) deallocate(idx2C)
    if (allocated(idx1C)) deallocate(idx1C)
    if (allocated(valC)) deallocate(valC)
    if (fmtC=='csc') then
       allocate(idx2C(n+1))
    else
       allocate(idx2C(nnz))
    end if
    allocate(idx1C(nnz))
    allocate(valC(nnz))
    call psp_list2spm(m,n,idx1C,idx2C,valC,fmtC,list,nnz,.false.)

    deallocate(valCLoc)
    deallocate(idx1CLoc)
    deallocate(idx2CLoc)

    if (nnz>0) call list_destroy(list)

  end subroutine psp_sst_zsum_spmspm




  subroutine psp_sst_dgespmm(M,N,K,opA,opB,alpha,row_ind,col_ptr,val,B,IB,JB,C,IC,JC,beta)
    ! sequential sparse matrix (ST format) dense matrix multiplication
    ! C(IC:IC+M-1,JC:JC+N-1) = alpha*opA(A)(1:M,1:K)*opB(B)(IB:IB+K-1,JB:JB+N-1) + beta*C(IC:IC+M-1,JC:JC+N-1)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K, IB, JB, IC, JC
    integer, intent(in) :: row_ind(:), col_ptr(:)
    real(dp), intent(in) :: val(:), B(:,:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: C(:,:)

    !**** LOCAL ***********************************!
    integer :: i, idx_col, cnt, IB_1, JB_1, IC_1, JC_1
    integer :: trA, trB, ot
    real(dp), allocatable :: C_loc(:,:)

    !**********************************************!
    IB_1=IB-1
    JB_1=JB-1
    IC_1=IC-1
    JC_1=JC-1

    if (alpha/=0.0_dp) then
       allocate(C_loc(M,N))
       C_loc=0.0_dp

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
          call die('invalid implementation')
       end if
       if (IB_1==0 .and. JB_1==0) then
          select case (ot)
          case (1)
             if (alpha/=1.0_dp) then
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      C_loc(row_ind(i),1:N)= alpha*val(i)*B(idx_col,1:N) &
                           + C_loc(row_ind(i),1:N)
                   end do
                end do
             else
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      C_loc(row_ind(i),1:N)= val(i)*B(idx_col,1:N) &
                           + C_loc(row_ind(i),1:N)
                   end do
                end do
             end if

          case (2)
             if (alpha/=1.0_dp) then
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      do cnt = 1,N
                         C_loc(row_ind(i),cnt)= alpha*val(i)*B(cnt,idx_col) &
                              + C_loc(row_ind(i),cnt)
                      end do
                   end do
                end do
             else
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      do cnt = 1,N
                         C_loc(row_ind(i),cnt)= val(i)*B(cnt,idx_col) &
                              + C_loc(row_ind(i),cnt)
                      end do
                   end do
                end do
             end if

          case (3)
             if (alpha/=1.0_dp) then
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      C_loc(idx_col,1:N)= alpha*val(i)*B(row_ind(i),1:N) &
                           + C_loc(idx_col,1:N)
                   end do
                end do
             else
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      C_loc(idx_col,1:N)= val(i)*B(row_ind(i)+IB_1,1:N) &
                           + C_loc(idx_col,1:N)
                   end do
                end do
             end if
          case (4)
             if (alpha/=1.0_dp) then
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      do cnt = 1,N
                         C_loc(idx_col,cnt)= alpha*val(i)*B(cnt,row_ind(i)) &
                              + C_loc(idx_col,cnt)
                      end do
                   end do
                end do
             else
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      do cnt = 1,N
                         C_loc(idx_col,cnt)= val(i)*B(cnt,row_ind(i)) &
                              + C_loc(idx_col,cnt)
                      end do
                   end do
                end do
             end if
          end select
       else
          select case (ot)
          case (1)
             if (alpha/=1.0_dp) then
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      C_loc(row_ind(i),1:N)= alpha*val(i)*B(idx_col+IB_1,JB:JB_1+N) &
                           + C_loc(row_ind(i),1:N)
                   end do
                end do
             else
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      C_loc(row_ind(i),1:N)= val(i)*B(idx_col+IB_1,JB:JB_1+N) &
                           + C_loc(row_ind(i),1:N)
                   end do
                end do
             end if

          case (2)
             if (alpha/=1.0_dp) then
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      do cnt = JB,JB_1+N
                         C_loc(row_ind(i),cnt-JB_1)= alpha*val(i)*B(cnt,idx_col+IB_1) &
                              + C_loc(row_ind(i),cnt-JB_1)
                      end do
                   end do
                end do
             else
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      do cnt = JB,JB_1+N
                         C_loc(row_ind(i),cnt-JB_1)= val(i)*B(cnt,idx_col+IB_1) &
                              + C_loc(row_ind(i),cnt-JB_1)
                      end do
                   end do
                end do
             end if

          case (3)
             if (alpha/=1.0_dp) then
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      C_loc(idx_col,1:N)= alpha*val(i)*B(row_ind(i)+IB_1,JB:JB_1+N) &
                           + C_loc(idx_col,1:N)
                   end do
                end do
             else
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      C_loc(idx_col,1:N)= val(i)*B(row_ind(i)+IB_1,JB:JB_1+N) &
                           + C_loc(idx_col,1:N)
                   end do
                end do
             end if
          case (4)
             if (alpha/=1.0_dp) then
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      do cnt = JB,JB_1+N
                         C_loc(idx_col,cnt-JB_1)= alpha*val(i)*B(cnt,row_ind(i)+IB_1) &
                              + C_loc(idx_col,cnt-JB_1)
                      end do
                   end do
                end do
             else
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      do cnt = JB,JB_1+N
                         C_loc(idx_col,cnt-JB_1)= val(i)*B(cnt,row_ind(i)+IB_1) &
                              + C_loc(idx_col,cnt-JB_1)
                      end do
                   end do
                end do
             end if
          end select
       end if
       if (IC==1 .and. JC==1) then
          if (beta==1.0_dp) then
             do cnt=1,N
                do i=1,M
                   C(i,cnt)=C_loc(i,cnt)+C(i,cnt)
                end do
             end do
          else if (beta==0.0_dp) then
             do cnt=1,N
                do i=1,M
                   C(i,cnt)=C_loc(i,cnt)
                end do
             end do
          else
             do cnt=1,N
                do i=1,M
                   C(i,cnt)=C_loc(i,cnt)+beta*C(i,cnt)
                end do
             end do
          end if
       else
          if (beta==1.0_dp) then
             do cnt=1,N
                do i=1,M
                   C(i+IC_1,cnt+JC_1)=C_loc(i,cnt)+C(i+IC_1,cnt+JC_1)
                end do
             end do
          else if (beta==0.0_dp) then
             do cnt=1,N
                do i=1,M
                   C(i+IC_1,cnt+JC_1)=C_loc(i,cnt)
                end do
             end do
          else
             do cnt=1,N
                do i=1,M
                   C(i+IC_1,cnt+JC_1)=C_loc(i,cnt)+beta*C(i+IC_1,cnt+JC_1)
                end do
             end do
          end if
       end if

       deallocate(C_loc)

    else
       if (IC==1 .and. JC==1) then
          if (beta/=1.0_dp) then
             if (beta==0.0_dp) then
                do cnt=1,N
                   do i=1,M
                      C(i,cnt)=0.0_dp
                   end do
                end do
             else
                do cnt=1,N
                   do i=1,M
                      C(i,cnt)=beta*C(i,cnt)
                   end do
                end do
             end if
          end if
       else
          if (beta/=1.0_dp) then
             if (beta==0.0_dp) then
                do cnt=1,N
                   do i=1,M
                      C(i+IC_1,cnt+JC_1)=0.0_dp
                   end do
                end do
             else
                do cnt=1,N
                   do i=1,M
                      C(i+IC_1,cnt+JC_1)=beta*C(i+IC_1,cnt+JC_1)
                   end do
                end do
             end if
          end if
       end if

    end if



  end subroutine psp_sst_dgespmm

  subroutine psp_sst_zgespmm(M,N,K,opA,opB,alpha,row_ind,col_ptr,val,B,IB,JB,C,IC,JC,beta)
    ! sequential sparse matrix (ST format) dense matrix multiplication
    ! C(IC:IC+M-1,JC:JC+N-1) = alpha*opA(A)(1:M,1:K)*opB(B)(IB:IB+K-1,JB:JB+N-1) + beta*C(IC:IC+M-1,JC:JC+N-1)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K, IB, JB, IC, JC
    integer, intent(in) :: row_ind(:), col_ptr(:)
    complex(dp), intent(in) :: val(:), B(:,:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: C(:,:)

    !**** LOCAL ***********************************!
    integer :: i, idx_col, cnt, IB_1, JB_1, IC_1, JC_1
    integer :: trA, trB, ot
    complex(dp), allocatable :: C_loc(:,:)

    !**********************************************!
    IB_1=IB-1
    JB_1=JB-1
    IC_1=IC-1
    JC_1=JC-1
    if (alpha/=cmplx_0) then
       allocate(C_loc(M,N))
       C_loc=0.0_dp

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
          call die('invalid implementation')
       end if
       if (IB_1==0 .and. JB_1==0) then
          select case (ot)
          case (1)
             if (alpha/=cmplx_1) then
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=M) then
                         C_loc(row_ind(i),1:N)= alpha*val(i)*B(idx_col,1:N) &
                              + C_loc(row_ind(i),1:N)
                      end if
                   end do
                end do
             else
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=M) then
                         C_loc(row_ind(i),1:N)= val(i)*B(idx_col,1:N) &
                              + C_loc(row_ind(i),1:N)
                      end if
                   end do
                end do
             end if
          case (2)
             if (alpha/=cmplx_1) then
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=M) then
                         do cnt = 1,N
                            if (trB==1) then
                               C_loc(row_ind(i),cnt)= alpha*val(i)*CONJG(B(cnt,idx_col)) &
                                    + C_loc(row_ind(i),cnt)
                            else
                               C_loc(row_ind(i),cnt)= alpha*val(i)*B(cnt,idx_col) &
                                    + C_loc(row_ind(i),cnt)
                            end if
                         end do
                      end if
                   end do
                end do
             else
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=M) then
                         do cnt = 1,N
                            if (trB==1) then
                               C_loc(row_ind(i),cnt)= val(i)*CONJG(B(cnt,idx_col)) &
                                    + C_loc(row_ind(i),cnt)
                            else
                               C_loc(row_ind(i),cnt)= val(i)*B(cnt,idx_col) &
                                    + C_loc(row_ind(i),cnt)
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          case (3)
             if (alpha/=cmplx_1) then
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         if (trA==1) then
                            C_loc(idx_col,1:N)= alpha*CONJG(val(i))*B(row_ind(i),1:N) &
                                 + C_loc(idx_col,1:N)
                         else
                            C_loc(idx_col,1:N)= alpha*val(i)*B(row_ind(i),1:N) &
                                 + C_loc(idx_col,1:N)
                         end if
                      end if
                   end do
                end do
             else
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         if (trA==1) then
                            C_loc(idx_col,1:N)= CONJG(val(i))*B(row_ind(i),1:N) &
                                 + C_loc(idx_col,1:N)
                         else
                            C_loc(idx_col,1:N)= val(i)*B(row_ind(i),1:N) &
                                 + C_loc(idx_col,1:N)
                         end if
                      end if
                   end do
                end do
             end if
          case (4)
             if (alpha/=cmplx_1) then
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         do cnt = 1,N
                            if (trA==2 .and. trB==1) then!if (trA==1 .and. trB==2) then
                               C_loc(idx_col,cnt)= alpha*val(i)*CONJG(B(cnt,row_ind(i))) &
                                    + C_loc(idx_col,cnt)
                            else if (trA==1 .and. trB==2) then
                               C_loc(idx_col,cnt)= alpha*CONJG(val(i))*B(cnt,row_ind(i)) &
                                    + C_loc(idx_col,cnt)
                            else if (trA==2 .and. trB==2) then
                               C_loc(idx_col,cnt)= alpha*val(i)*B(cnt,row_ind(i)) &
                                    + C_loc(idx_col,cnt)
                            else ! trA==1 .and. trB==1
                               C_loc(idx_col,cnt)= alpha*CONJG(val(i)*B(cnt,row_ind(i))) &
                                    + C_loc(idx_col,cnt)
                            end if
                         end do
                      end if
                   end do
                end do
             else
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         do cnt = 1,N
                            if (trA==2 .and. trB==1) then!if (trA==1 .and. trB==2) then
                               C_loc(idx_col,cnt)= val(i)*CONJG(B(cnt,row_ind(i))) &
                                    + C_loc(idx_col,cnt)
                            else if (trA==1 .and. trB==2) then
                               C_loc(idx_col,cnt)= CONJG(val(i))*B(cnt,row_ind(i)) &
                                    + C_loc(idx_col,cnt)
                            else if (trA==2 .and. trB==2) then
                               C_loc(idx_col,cnt)= val(i)*B(cnt,row_ind(i)) &
                                    + C_loc(idx_col,cnt)
                            else ! trA==1 .and. trB==1
                               C_loc(idx_col,cnt)= CONJG(val(i)*B(cnt,row_ind(i))) &
                                    + C_loc(idx_col,cnt)
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          end select
       else
          select case (ot)
          case (1)
             if (alpha/=cmplx_1) then
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=M) then
                         C_loc(row_ind(i),1:N)= alpha*val(i)*B(idx_col+IB_1,JB:JB_1+N) &
                              + C_loc(row_ind(i),1:N)
                      end if
                   end do
                end do
             else
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=M) then
                         C_loc(row_ind(i),1:N)= val(i)*B(idx_col+IB_1,JB:JB_1+N) &
                              + C_loc(row_ind(i),1:N)
                      end if
                   end do
                end do
             end if
          case (2)
             if (alpha/=cmplx_1) then
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=M) then
                         do cnt = JB,JB_1+N
                            if (trB==1) then
                               C_loc(row_ind(i),cnt-JB_1)= alpha*val(i)*CONJG(B(cnt,idx_col+IB_1)) &
                                    + C_loc(row_ind(i),cnt-JB_1)
                            else
                               C_loc(row_ind(i),cnt-JB_1)= alpha*val(i)*B(cnt,idx_col+IB_1) &
                                    + C_loc(row_ind(i),cnt-JB_1)
                            end if
                         end do
                      end if
                   end do
                end do
             else
                do idx_col=1,K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=M) then
                         do cnt = JB,JB_1+N
                            if (trB==1) then
                               C_loc(row_ind(i),cnt-JB_1)= val(i)*CONJG(B(cnt,idx_col+IB_1)) &
                                    + C_loc(row_ind(i),cnt-JB_1)
                            else
                               C_loc(row_ind(i),cnt-JB_1)= val(i)*B(cnt,idx_col+IB_1) &
                                    + C_loc(row_ind(i),cnt-JB_1)
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          case (3)
             if (alpha/=cmplx_1) then
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         if (trA==1) then
                            C_loc(idx_col,1:N)= alpha*CONJG(val(i))*B(row_ind(i)+IB_1,JB:JB_1+N) &
                                 + C_loc(idx_col,1:N)
                         else
                            C_loc(idx_col,1:N)= alpha*val(i)*B(row_ind(i)+IB_1,JB:JB_1+N) &
                                 + C_loc(idx_col,1:N)
                         end if
                      end if
                   end do
                end do
             else
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         if (trA==1) then
                            C_loc(idx_col,1:N)= CONJG(val(i))*B(row_ind(i)+IB_1,JB:JB_1+N) &
                                 + C_loc(idx_col,1:N)
                         else
                            C_loc(idx_col,1:N)= val(i)*B(row_ind(i)+IB_1,JB:JB_1+N) &
                                 + C_loc(idx_col,1:N)
                         end if
                      end if
                   end do
                end do
             end if
          case (4)
             if (alpha/=cmplx_1) then
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         do cnt = JB,JB_1+N
                            if (trA==2 .and. trB==1) then!if (trA==1 .and. trB==2) then
                               C_loc(idx_col,cnt-JB_1)= alpha*val(i)*CONJG(B(cnt,row_ind(i)+IB_1)) &
                                    + C_loc(idx_col,cnt-JB_1)
                            else if (trA==1 .and. trB==2) then
                               C_loc(idx_col,cnt-JB_1)= alpha*CONJG(val(i))*B(cnt,row_ind(i)+IB_1) &
                                    + C_loc(idx_col,cnt-JB_1)
                            else if (trA==2 .and. trB==2) then
                               C_loc(idx_col,cnt-JB_1)= alpha*val(i)*B(cnt,row_ind(i)+IB_1) &
                                    + C_loc(idx_col,cnt-JB_1)
                            else ! trA==1 .and. trB==1
                               C_loc(idx_col,cnt-JB_1)= alpha*CONJG(val(i)*B(cnt,row_ind(i)+IB_1)) &
                                    + C_loc(idx_col,cnt-JB_1)
                            end if
                         end do
                      end if
                   end do
                end do
             else
                do idx_col=1,M!K
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         do cnt = JB,JB_1+N
                            if (trA==2 .and. trB==1) then!if (trA==1 .and. trB==2) then
                               C_loc(idx_col,cnt-JB_1)= val(i)*CONJG(B(cnt,row_ind(i)+IB_1)) &
                                    + C_loc(idx_col,cnt-JB_1)
                            else if (trA==1 .and. trB==2) then
                               C_loc(idx_col,cnt-JB_1)= CONJG(val(i))*B(cnt,row_ind(i)+IB_1) &
                                    + C_loc(idx_col,cnt-JB_1)
                            else if (trA==2 .and. trB==2) then
                               C_loc(idx_col,cnt-JB_1)= val(i)*B(cnt,row_ind(i)+IB_1) &
                                    + C_loc(idx_col,cnt-JB_1)
                            else ! trA==1 .and. trB==1
                               C_loc(idx_col,cnt-JB_1)= CONJG(val(i)*B(cnt,row_ind(i)+IB_1)) &
                                    + C_loc(idx_col,cnt-JB_1)
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          end select

       end if
       if (IC==1 .and. JC==1) then
          if (beta==cmplx_1) then
             do cnt=1,N
                do i=1,M
                   C(i,cnt)=C_loc(i,cnt)+C(i,cnt)
                end do
             end do
          else if (beta==cmplx_0) then
             do cnt=1,N
                do i=1,M
                   C(i,cnt)=C_loc(i,cnt)
                end do
             end do
          else
             do cnt=1,N
                do i=1,M
                   C(i,cnt)=C_loc(i,cnt)+beta*C(i,cnt)
                end do
             end do
          end if
       else
          if (beta==cmplx_1) then
             do cnt=1,N
                do i=1,M
                   C(i+IC_1,cnt+JC_1)=C_loc(i,cnt)+C(i+IC_1,cnt+JC_1)
                end do
             end do
          else if (beta==cmplx_0) then
             do cnt=1,N
                do i=1,M
                   C(i+IC_1,cnt+JC_1)=C_loc(i,cnt)
                end do
             end do
          else
             do cnt=1,N
                do i=1,M
                   C(i+IC_1,cnt+JC_1)=C_loc(i,cnt)+beta*C(i+IC_1,cnt+JC_1)
                end do
             end do
          end if
       end if

       deallocate(C_loc)

    else
       if (IC==1 .and. JC==1) then
          if (beta/=cmplx_1) then
             if (beta==cmplx_0) then
                do cnt=1,N
                   do i=1,M
                      C(i,cnt)=cmplx_0
                   end do
                end do
             else
                do cnt=1,N
                   do i=1,M
                      C(i,cnt)=beta*C(i,cnt)
                   end do
                end do
             end if
          end if
       else
          if (beta/=cmplx_1) then
             if (beta==cmplx_0) then
                do cnt=1,N
                   do i=1,M
                      C(i+IC_1,cnt+JC_1)=cmplx_0
                   end do
                end do
             else
                do cnt=1,N
                   do i=1,M
                      C(i+IC_1,cnt+JC_1)=beta*C(i+IC_1,cnt+JC_1)
                   end do
                end do
             end if
          end if
       end if
    end if

  end subroutine psp_sst_zgespmm

  subroutine psp_sst_dgemspm(M,N,K,opA,opB,alpha,A,IA,JA,row_ind,col_ptr,val,C,IC,JC,beta)
    ! sequential sparse matrix (ST format) dense matrix multiplication
    ! C(IC:IC+M-1,JC:JC+N-1) = alpha*opA(A)(IA:IA+M-1,JA:JA+K-1)*opB(B)(1:K,1:N) + beta*C(IC:IC+M-1,JC:JC+N-1)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K, IA, JA, IC, JC
    integer, intent(in) :: row_ind(:), col_ptr(:)
    real(dp), intent(in) :: val(:), A(:,:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: C(:,:)

    !**** LOCAL ***********************************!

    integer :: i, idx_col, cnt, IA_1, JA_1, IC_1,JC_1
    integer :: trA, trB, ot
    real(dp), allocatable :: C_loc(:,:)

    !**********************************************!
    IA_1=IA-1
    JA_1=JA-1
    IC_1=IC-1
    JC_1=JC-1

    if (alpha/=0.0_dp) then
       allocate(C_loc(M,N))
       C_loc=0.0_dp

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
          call die('invalid implementation')
       end if
       if (IA==1 .and. JA==1) then
          select case (ot)
          case (1)
             if (alpha==1.0_dp) then
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         C_loc(1:M,idx_col)= val(i)*A(1:M,row_ind(i)) &
                              + C_loc(1:M,idx_col)
                      end if
                   end do
                end do
             else
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         C_loc(1:M,idx_col)= alpha*val(i)*A(1:M,row_ind(i)) &
                              + C_loc(1:M,idx_col)
                      end if
                   end do
                end do
             end if

          case (2)
             if (alpha==1.0_dp) then
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         C_loc(1:M,row_ind(i))= val(i)*A(1:M,idx_col) &
                              + C_loc(1:M,row_ind(i))
                      end if
                   end do
                end do
             else
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         C_loc(1:M,row_ind(i))= alpha*val(i)*A(1:M,idx_col) &
                              + C_loc(1:M,row_ind(i))
                      end if
                   end do
                end do
             end if

          case (3)
             if (alpha==1.0_dp) then
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         do cnt=1,M
                            C_loc(cnt,idx_col)= val(i)*A(row_ind(i),cnt) &
                                 + C_loc(cnt,idx_col)
                         end do
                      end if
                   end do
                end do
             else
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         do cnt=1,M
                            C_loc(cnt,idx_col)= alpha*val(i)*A(row_ind(i),cnt) &
                                 + C_loc(cnt,idx_col)
                         end do
                      end if
                   end do
                end do
             end if

          case (4)
             if (alpha==1.0_dp) then
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         do cnt=1,M
                            C_loc(cnt,row_ind(i))= alpha*val(i)*A(idx_col,cnt) &
                                 + C_loc(cnt,row_ind(i))
                         end do
                      end if
                   end do
                end do
             else
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         do cnt=1,M
                            C_loc(cnt,row_ind(i))= val(i)*A(idx_col,cnt) &
                                 + C_loc(cnt,row_ind(i))
                         end do
                      end if
                   end do
                end do
             end if
          end select
       else
          select case (ot)
          case (1)
             if (alpha==1.0_dp) then
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         C_loc(1:M,idx_col)= val(i)*A(IA:IA_1+M,row_ind(i)+JA_1) &
                              + C_loc(1:M,idx_col)
                      end if
                   end do
                end do
             else
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         C_loc(1:M,idx_col)= alpha*val(i)*A(IA:IA_1+M,row_ind(i)+JA_1) &
                              + C_loc(1:M,idx_col)
                      end if
                   end do
                end do
             end if

          case (2)
             if (alpha==1.0_dp) then
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         C_loc(1:M,row_ind(i))= val(i)*A(IA:IA_1+M,idx_col+JA_1) &
                              + C_loc(1:M,row_ind(i))
                      end if
                   end do
                end do
             else
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         C_loc(1:M,row_ind(i))= alpha*val(i)*A(IA:IA_1+M,idx_col+JA_1) &
                              + C_loc(1:M,row_ind(i))
                      end if
                   end do
                end do
             end if

          case (3)
             if (alpha==1.0_dp) then
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         do cnt=IA,IA_1+M
                            C_loc(cnt-IA_1,idx_col)= val(i)*A(row_ind(i)+JA_1,cnt) &
                                 + C_loc(cnt-IA_1,idx_col)
                         end do
                      end if
                   end do
                end do
             else
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         do cnt=IA,IA_1+M
                            C_loc(cnt-IA_1,idx_col)= alpha*val(i)*A(row_ind(i)+JA_1,cnt) &
                                 + C_loc(cnt-IA_1,idx_col)
                         end do
                      end if
                   end do
                end do
             end if

          case (4)
             if (alpha==1.0_dp) then
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         do cnt=IA,IA_1+M
                            C_loc(cnt-IA_1,row_ind(i))= alpha*val(i)*A(idx_col+JA_1,cnt) &
                                 + C_loc(cnt-IA_1,row_ind(i))
                         end do
                      end if
                   end do
                end do
             else
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         do cnt=IA,IA_1+M
                            C_loc(cnt-IA_1,row_ind(i))= val(i)*A(idx_col+JA_1,cnt) &
                                 + C_loc(cnt-IA_1,row_ind(i))
                         end do
                      end if
                   end do
                end do
             end if
          end select

       end if
       if (IC==1 .and. JC==1) then
          if (beta==1.0_dp) then
             do cnt=1,N
                do i=1,M
                   C(i,cnt)=C_loc(i,cnt)+C(i,cnt)
                end do
             end do
          else if (beta==0.0_dp) then
             do cnt=1,N
                do i=1,M
                   C(i,cnt)=C_loc(i,cnt)
                end do
             end do
          else
             do cnt=1,N
                do i=1,M
                   C(i,cnt)=C_loc(i,cnt)+beta*C(i,cnt)
                end do
             end do
          end if
          deallocate(C_loc)
       else
          if (beta==1.0_dp) then
             do cnt=1,N
                do i=1,M
                   C(i+IC_1,cnt+JC_1)=C_loc(i,cnt)+C(i+IC_1,cnt+JC_1)
                end do
             end do
          else if (beta==0.0_dp) then
             do cnt=1,N
                do i=1,M
                   C(i+IC_1,cnt+JC_1)=C_loc(i,cnt)
                end do
             end do
          else
             do cnt=1,N
                do i=1,M
                   C(i+IC_1,cnt+JC_1)=C_loc(i,cnt)+beta*C(i+IC_1,cnt+JC_1)
                end do
             end do
          end if
          deallocate(C_loc)
       end if

    else
       if (IC==1 .and. JC==1) then
          if (beta/=1.0_dp) then
             if (beta==0.0_dp) then
                do cnt=1,N
                   do i=1,M
                      C(i,cnt)=0.0_dp
                   end do
                end do
             else
                do cnt=1,N
                   do i=1,M
                      C(i,cnt)=beta*C(i,cnt)
                   end do
                end do
             end if
          end if
       else
          if (beta/=1.0_dp) then
             if (beta==0.0_dp) then
                do cnt=1,N
                   do i=1,M
                      C(i+IC_1,cnt+JC_1)=0.0_dp
                   end do
                end do
             else
                do cnt=1,N
                   do i=1,M
                      C(i+IC_1,cnt+JC_1)=beta*C(i+IC_1,cnt+JC_1)
                   end do
                end do
             end if
          end if
       end if

    end if

  end subroutine psp_sst_dgemspm


  subroutine psp_sst_zgemspm(M,N,K,opA,opB,alpha,A,IA,JA,row_ind,col_ptr,val,C,IC,JC,beta)
    ! sequential sparse matrix (ST format) dense matrix multiplication
    ! C(IC:IC+M-1,JC:JC+N-1) = alpha*opA(A)(IA:IA+M-1,JA:JA+K-1)*opB(B)(1:K,1:N) + beta*C(IC:IC+M-1,JC:JC+N-1)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K, IA, JA, IC, JC
    integer, intent(in) :: row_ind(:), col_ptr(:)
    complex(dp), intent(in) :: val(:), A(:,:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: C(:,:)

    !**** LOCAL ***********************************!

    integer :: i, idx_col, cnt, IA_1, JA_1, IC_1,JC_1

    integer :: trA, trB, ot
    complex(dp), allocatable :: C_loc(:,:)

    !**********************************************!
    IA_1=IA-1
    JA_1=JA-1
    IC_1=IC-1
    JC_1=JC-1

    if (alpha/=cmplx_0) then
       allocate(C_loc(M,N))
       C_loc=cmplx_0

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
          call die('invalid implementation')
       end if
       if (IA==1 .and. JA==1) then
          select case (ot)
          case (1)
             if (alpha/=cmplx_1) then
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         C_loc(1:M,idx_col)= alpha*val(i)*A(1:M,row_ind(i)) &
                              + C_loc(1:M,idx_col)
                      end if
                   end do
                end do
             else
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         C_loc(1:M,idx_col)= val(i)*A(1:M,row_ind(i)) &
                              + C_loc(1:M,idx_col)
                      end if
                   end do
                end do
             end if
          case (2)
             if (alpha/=cmplx_1) then
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         if (trB==1) then
                            C_loc(1:M,row_ind(i))= alpha*CONJG(val(i))*A(1:M,idx_col) &
                                 + C_loc(1:M,row_ind(i))
                         else
                            C_loc(1:M,row_ind(i))= alpha*val(i)*A(1:M,idx_col) &
                                 + C_loc(1:M,row_ind(i))
                         end if
                      end if
                   end do
                end do
             else
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         if (trB==1) then
                            C_loc(1:M,row_ind(i))= CONJG(val(i))*A(1:M,idx_col) &
                                 + C_loc(1:M,row_ind(i))
                         else
                            C_loc(1:M,row_ind(i))= val(i)*A(1:M,idx_col) &
                                 + C_loc(1:M,row_ind(i))
                         end if
                      end if
                   end do
                end do
             end if
          case (3)
             if (alpha/=cmplx_1) then
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         do cnt=1,M
                            if (trA==1) then
                               C_loc(cnt,idx_col)= alpha*val(i)*CONJG(A(row_ind(i),cnt)) &
                                    + C_loc(cnt,idx_col)
                            else
                               C_loc(cnt,idx_col)= alpha*val(i)*A(row_ind(i),cnt) &
                                    + C_loc(cnt,idx_col)
                            end if
                         end do
                      end if
                   end do
                end do
             else
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         do cnt=1,M
                            if (trA==1) then
                               C_loc(cnt,idx_col)= val(i)*CONJG(A(row_ind(i),cnt)) &
                                    + C_loc(cnt,idx_col)
                            else
                               C_loc(cnt,idx_col)= val(i)*A(row_ind(i),cnt) &
                                    + C_loc(cnt,idx_col)
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          case (4)
             if (alpha/=cmplx_1) then
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         do cnt=1,M
                            if (trA==1 .and. trB==2) then
                               C_loc(cnt,row_ind(i))= alpha*val(i)*CONJG(A(idx_col,cnt)) &
                                    + C_loc(cnt,row_ind(i))
                            else if (trA==2 .and. trB==1) then
                               C_loc(cnt,row_ind(i))= alpha*CONJG(val(i))*A(idx_col,cnt) &
                                    + C_loc(cnt,row_ind(i))
                            else if (trA==1 .and. trB==1) then
                               C_loc(cnt,row_ind(i))= alpha*CONJG(val(i)*A(idx_col,cnt)) &
                                    + C_loc(cnt,row_ind(i))
                            else
                               C_loc(cnt,row_ind(i))= alpha*val(i)*A(idx_col,cnt) &
                                    + C_loc(cnt,row_ind(i))
                            end if
                         end do
                      end if
                   end do
                end do
             else
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         do cnt=1,M
                            if (trA==1 .and. trB==2) then
                               C_loc(cnt,row_ind(i))= val(i)*CONJG(A(idx_col,cnt)) &
                                    + C_loc(cnt,row_ind(i))
                            else if (trA==2 .and. trB==1) then
                               C_loc(cnt,row_ind(i))= CONJG(val(i))*A(idx_col,cnt) &
                                    + C_loc(cnt,row_ind(i))
                            else if (trA==1 .and. trB==1) then
                               C_loc(cnt,row_ind(i))= CONJG(val(i)*A(idx_col,cnt)) &
                                    + C_loc(cnt,row_ind(i))
                            else
                               C_loc(cnt,row_ind(i))= val(i)*A(idx_col,cnt) &
                                    + C_loc(cnt,row_ind(i))
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          end select
       else
          select case (ot)
          case (1)
             if (alpha/=cmplx_1) then
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         C_loc(1:M,idx_col)= alpha*val(i)*A(IA:IA_1+M,row_ind(i)+JA_1) &
                              + C_loc(1:M,idx_col)
                      end if
                   end do
                end do
             else
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         C_loc(1:M,idx_col)= val(i)*A(IA:IA_1+M,row_ind(i)+JA_1) &
                              + C_loc(1:M,idx_col)
                      end if
                   end do
                end do
             end if
          case (2)
             if (alpha/=cmplx_1) then
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         if (trB==1) then
                            C_loc(1:M,row_ind(i))= alpha*CONJG(val(i))*A(IA:IA_1+M,idx_col+JA_1) &
                                 + C_loc(1:M,row_ind(i))
                         else
                            C_loc(1:M,row_ind(i))= alpha*val(i)*A(IA:IA_1+M,idx_col+JA_1) &
                                 + C_loc(1:M,row_ind(i))
                         end if
                      end if
                   end do
                end do
             else
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         if (trB==1) then
                            C_loc(1:M,row_ind(i))= CONJG(val(i))*A(IA:IA_1+M,idx_col+JA_1) &
                                 + C_loc(1:M,row_ind(i))
                         else
                            C_loc(1:M,row_ind(i))= val(i)*A(IA:IA_1+M,idx_col+JA_1) &
                                 + C_loc(1:M,row_ind(i))
                         end if
                      end if
                   end do
                end do
             end if
          case (3)
             if (alpha/=cmplx_1) then
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         do cnt=IA,IA_1+M
                            if (trA==1) then
                               C_loc(cnt-IA_1,idx_col)= alpha*val(i)*CONJG(A(row_ind(i)+JA_1,cnt)) &
                                    + C_loc(cnt-IA_1,idx_col)
                            else
                               C_loc(cnt-IA_1,idx_col)= alpha*val(i)*A(row_ind(i)+JA_1,cnt) &
                                    + C_loc(cnt-IA_1,idx_col)
                            end if
                         end do
                      end if
                   end do
                end do
             else
                do idx_col=1,N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=K) then
                         do cnt=IA,IA_1+M
                            if (trA==1) then
                               C_loc(cnt-IA_1,idx_col)= val(i)*CONJG(A(row_ind(i)+JA_1,cnt)) &
                                    + C_loc(cnt-IA_1,idx_col)
                            else
                               C_loc(cnt-IA_1,idx_col)= val(i)*A(row_ind(i)+JA_1,cnt) &
                                    + C_loc(cnt-IA_1,idx_col)
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          case (4)
             if (alpha/=cmplx_1) then
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         do cnt=IA,IA_1+M
                            if (trA==1 .and. trB==2) then
                               C_loc(cnt-IA_1,row_ind(i))= alpha*val(i)*CONJG(A(idx_col+JA_1,cnt)) &
                                    + C_loc(cnt-IA_1,row_ind(i))
                            else if (trA==2 .and. trB==1) then
                               C_loc(cnt-IA_1,row_ind(i))= alpha*CONJG(val(i))*A(idx_col+JA_1,cnt) &
                                    + C_loc(cnt-IA_1,row_ind(i))
                            else if (trA==1 .and. trB==1) then
                               C_loc(cnt-IA_1,row_ind(i))= alpha*CONJG(val(i)*A(idx_col+JA_1,cnt)) &
                                    + C_loc(cnt-IA_1,row_ind(i))
                            else
                               C_loc(cnt-IA_1,row_ind(i))= alpha*val(i)*A(idx_col+JA_1,cnt) &
                                    + C_loc(cnt-IA_1,row_ind(i))
                            end if
                         end do
                      end if
                   end do
                end do
             else
                do idx_col=1,K!N
                   do i=col_ptr(idx_col),col_ptr(idx_col+1)-1
                      if (row_ind(i)<=N) then
                         do cnt=IA,IA_1+M
                            if (trA==1 .and. trB==2) then
                               C_loc(cnt-IA_1,row_ind(i))= val(i)*CONJG(A(idx_col+JA_1,cnt)) &
                                    + C_loc(cnt-IA_1,row_ind(i))
                            else if (trA==2 .and. trB==1) then
                               C_loc(cnt-IA_1,row_ind(i))= CONJG(val(i))*A(idx_col+JA_1,cnt) &
                                    + C_loc(cnt-IA_1,row_ind(i))
                            else if (trA==1 .and. trB==1) then
                               C_loc(cnt-IA_1,row_ind(i))= CONJG(val(i)*A(idx_col+JA_1,cnt)) &
                                    + C_loc(cnt-IA_1,row_ind(i))
                            else
                               C_loc(cnt-IA_1,row_ind(i))= val(i)*A(idx_col+JA_1,cnt) &
                                    + C_loc(cnt-IA_1,row_ind(i))
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          end select
       end if
       if (IC==1 .and. JC==1) then
          if (beta==cmplx_1) then
             do cnt=1,N
                do i=1,M
                   C(i+IC_1,cnt+JC_1)=C_loc(i,cnt)+C(i+IC_1,cnt+JC_1)
                end do
             end do
          else if (beta==cmplx_0) then
             do cnt=1,N
                do i=1,M
                   C(i+IC_1,cnt+JC_1)=C_loc(i,cnt)
                end do
             end do
          else
             do cnt=1,N
                do i=1,M
                   C(i+IC_1,cnt+JC_1)=C_loc(i,cnt)+beta*C(i+IC_1,cnt+JC_1)
                end do
             end do
          end if
       else
          if (beta==cmplx_1) then
             do cnt=1,N
                do i=1,M
                   C(i,cnt)=C_loc(i,cnt)+C(i,cnt)
                end do
             end do
          else if (beta==cmplx_0) then
             do cnt=1,N
                do i=1,M
                   C(i,cnt)=C_loc(i,cnt)
                end do
             end do
          else
             do cnt=1,N
                do i=1,M
                   C(i,cnt)=C_loc(i,cnt)+beta*C(i,cnt)
                end do
             end do
          end if
       end if
       deallocate(C_loc)

    else
       if (IC==1 .and. JC==1) then
          if (beta/=cmplx_1) then
             if (beta==cmplx_0) then
                do cnt=1,N
                   do i=1,M
                      C(i+IC_1,cnt+JC_1)=cmplx_0
                   end do
                end do
             else
                do cnt=1,N
                   do i=1,M
                      C(i+IC_1,cnt+JC_1)=beta*C(i+IC_1,cnt+JC_1)
                   end do
                end do
             end if
          end if
       else
          if (beta/=cmplx_1) then
             if (beta==cmplx_0) then
                do cnt=1,N
                   do i=1,M
                      C(i,cnt)=cmplx_0
                   end do
                end do
             else
                do cnt=1,N
                   do i=1,M
                      C(i,cnt)=beta*C(i,cnt)
                   end do
                end do
             end if
          end if
       end if

    end if

  end subroutine psp_sst_zgemspm

  subroutine psp_sst_dgespmspm(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
       row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    ! sequential sparse matrix sparse matrix multiplication ('csc' format)
    ! C(1:M,1:N) = alpha*opA(A)(1:M,1:K)*opB(B)(1:K,1:N)+beta*C(1:M,1:N)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K
    integer, intent(in) :: row_ind_A(:), col_ptr_A(:)
    integer, intent(in) :: row_ind_B(:), col_ptr_B(:)
    real(dp), intent(in) :: val_A(:), val_B(:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: row_ind_C(:), col_ptr_C(:)
    real(dp), allocatable , intent(inout) :: val_C(:)

    !**** LOCAL ***********************************!
    integer :: trA, trB, ot
    !**********************************************!


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
       call die('invalid implementation')
    end if

    select case (ot)
    case (1)
       call psp_sst_dgespmspm_nn(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
            row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    case (2)
       call psp_sst_dgespmspm_nt(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
            row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    case (3)
       call psp_sst_dgespmspm_tn(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
            row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    case (4)
       call psp_sst_dgespmspm_tt(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
            row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    end select

  end subroutine psp_sst_dgespmspm

  subroutine psp_sst_zgespmspm(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
       row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    ! sequential sparse matrix sparse matrix multiplication ('csc' format)
    ! C(1:M,1:N) = alpha*opA(A)(1:M,1:K)*opB(B)(1:K,1:N)+beta*C(1:M,1:N)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K
    integer, intent(in) :: row_ind_A(:), col_ptr_A(:)
    integer, intent(in) :: row_ind_B(:), col_ptr_B(:)
    complex(dp), intent(in) :: val_A(:), val_B(:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: row_ind_C(:), col_ptr_C(:)
    complex(dp), allocatable , intent(inout) :: val_C(:)

    !**** LOCAL ***********************************!
    integer :: trA, trB, ot
    !**********************************************!


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
       call die('invalid implementation')
    end if

    select case (ot)
    case (1)
       call psp_sst_zgespmspm_nn(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
            row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    case (2)
       call psp_sst_zgespmspm_nt(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
            row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    case (3)
       call psp_sst_zgespmspm_tn(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
            row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    case (4)
       call psp_sst_zgespmspm_tt(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
            row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    end select

  end subroutine psp_sst_zgespmspm


  subroutine psp_sst_dgespmspm_nn(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
       row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    ! sequential sparse matrix sparse matrix multiplication ('csc' format)
    ! C(1:M,1:N) = alpha*A(1:M,1:K)*B(1:K,1:N)+beta*C(1:M,1:N)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K
    integer, intent(in) :: row_ind_A(:), col_ptr_A(:)
    integer, intent(in) :: row_ind_B(:), col_ptr_B(:)
    real(dp), intent(in) :: val_A(:), val_B(:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: row_ind_C(:), col_ptr_C(:)
    real(dp), allocatable , intent(inout) :: val_C(:)

    !**** LOCAL ***********************************!
    integer :: iB, jB, idx_col, cnt, cntA, cntAdded, crtGlb, crtLoc
    integer :: lenLoc, lenGlb
    real(dp) :: vB
    integer, allocatable :: row_ind_loc(:), col_ind_loc(:)
    real(dp), allocatable :: val_loc(:)
    type(dList), pointer :: listGlb, listLoc, lastGlb, lastLoc
    !**********************************************!

    allocate(val_loc(M))
    allocate(row_ind_loc(M))
    allocate(col_ind_loc(M))
    crtGlb=0
    lenGlb=0
    do jB=1,N ! loop over each column in B, correspondinly each column in C
       crtLoc=0
       lenLoc=0
       ! create a linked list
       do cnt=col_ptr_B(jB),col_ptr_B(jB+1)-1
          ! get the nonzero entry (iB,jB,vB) in B
          iB=row_ind_B(cnt)
          vB=val_B(cnt)
          ! get the iB'th column in A, the index range is col_ptr_A(iB):col_ptr_A(iB+1)-1
          cntAdded=0
          do cntA=col_ptr_A(iB),col_ptr_A(iB+1)-1
             cntAdded = cntAdded+1
             val_loc(cntAdded) = val_A(cntA)*vB*alpha
             row_ind_loc(cntAdded) = row_ind_A(cntA)
             col_ind_loc(cntAdded) = jB ! add to the jB'th column in C
          end do
          if (cntAdded>0) then
             crtLoc=crtLoc+1
             ! add the sparse matrix (row_ind_loc,col_ind_loc,val_loc) in linked list
             if (crtLoc==1) then
                call psp_list_create_mat(listLoc,lenLoc,1.0_dp,row_ind_loc,col_ind_loc,&
                     val_loc,'coo',cntAdded,lastLoc)
             else
                call psp_list_combine_listMat(M,N,1.0_dp,listLoc,lenLoc,1.0_dp,row_ind_loc,&
                     col_ind_loc,val_loc,'coo',cntAdded,lastLoc)
             end if
          end if
       end do
       if (crtLoc>0) then
          crtGlb=crtGlb+1
          if (crtGlb==1) then
             ! add the matrix stored in the linked list listLoc into listGlb, correspondingly to C
             listGlb=>listLoc ! TODO: any problem?
             lastGlb=>lastLoc
             lenGlb=lenLoc
          else
             ! add the matrix stored in the linked list listLoc into listGlb, correspondingly to C
             lastGlb%next => listLoc
             lastGlb => lastLoc
             lenGlb=lenGlb+lenLoc
          end if
       end if
       ! delete linked list listLoc
       listLoc => null()
       lastLoc => null()
    end do ! end do jB=1,N
    if (beta==0.0_dp) then
       ! convert format list to a sparse matrix with nnz nonzero elements
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    else
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (lenGlb==0) then
          call psp_list_create_mat(listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       else
          call psp_list_combine_listMat(M,N,1.0_dp,listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       end if
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (allocated(row_ind_C)) deallocate(row_ind_C)
       if (allocated(col_ptr_C)) deallocate(col_ptr_C)
       if (allocated(val_C)) deallocate(val_C)
       allocate(col_ptr_C(N+1))
       allocate(row_ind_C(lenGlb))
       allocate(val_C(lenGlb))
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    end if

    deallocate(val_loc)
    deallocate(row_ind_loc)
    deallocate(col_ind_loc)
    !call list_destroy(listGlb) ! why double free?

  end subroutine psp_sst_dgespmspm_nn

  subroutine psp_sst_zgespmspm_nn(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
       row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    ! sequential sparse matrix sparse matrix multiplication ('csc' format)
    ! C(1:M,1:N) = alpha*opA(A)(1:M,1:K)*opB(B)(1:K,1:N)+beta*C(1:M,1:N)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K
    integer, intent(in) :: row_ind_A(:), col_ptr_A(:)
    integer, intent(in) :: row_ind_B(:), col_ptr_B(:)
    complex(dp), intent(in) :: val_A(:), val_B(:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: row_ind_C(:), col_ptr_C(:)
    complex(dp), allocatable , intent(inout) :: val_C(:)

    !**** LOCAL ***********************************!
    integer :: iB, jB, idx_col, cnt, cntA, cntAdded, crtGlb, crtLoc
    integer :: trA, trB, ot, lenLoc, lenGlb
    complex(dp) :: vB
    integer, allocatable :: row_ind_loc(:), col_ind_loc(:)
    complex(dp), allocatable :: val_loc(:)
    type(zList), pointer :: listGlb, listLoc, lastGlb, lastLoc

    !**********************************************!

    allocate(val_loc(M))
    allocate(row_ind_loc(M))
    allocate(col_ind_loc(M))
    crtGlb=0
    lenGlb=0
    do jB=1,N ! loop over each column in B, correspondinly each column in C
       crtLoc=0
       lenLoc=0
       ! create a linked list
       do cnt=col_ptr_B(jB),col_ptr_B(jB+1)-1
          ! get the nonzero entry (iB,jB,vB) in B
          iB=row_ind_B(cnt)
          vB=val_B(cnt)
          ! get the iB'th column in A, the index range is col_ptr_A(iB):col_ptr_A(iB+1)-1
          cntAdded=0
          do cntA=col_ptr_A(iB),col_ptr_A(iB+1)-1
             cntAdded = cntAdded+1
             val_loc(cntAdded) = val_A(cntA)*vB*alpha
             row_ind_loc(cntAdded) = row_ind_A(cntA)
             col_ind_loc(cntAdded) = jB ! add to the jB'th column in C
          end do
          if (cntAdded>0) then
             crtLoc=crtLoc+1
             ! add the sparse matrix (row_ind_loc,col_ind_loc,val_loc) in linked list
             if (crtLoc==1) then
                call psp_list_create_mat(listLoc,lenLoc,cmplx_1,row_ind_loc,col_ind_loc,&
                     val_loc,'coo',cntAdded,lastLoc)
             else
                call psp_list_combine_listMat(M,N,cmplx_1,listLoc,lenLoc,cmplx_1,row_ind_loc,&
                     col_ind_loc,val_loc,'coo',cntAdded,lastLoc)
             end if
          end if
       end do
       if (crtLoc>0) then
          crtGlb=crtGlb+1
          if (crtGlb==1) then
             ! add the matrix stored in the linked list listLoc into listGlb, correspondingly to C
             listGlb=>listLoc ! TODO: any problem?
             lastGlb=>lastLoc
             lenGlb=lenLoc
          else
             ! add the matrix stored in the linked list listLoc into listGlb, correspondingly to C
             lastGlb%next => listLoc
             lastGlb => lastLoc
             lenGlb=lenGlb+lenLoc
          end if
       end if
       ! delete linked list listLoc
       listLoc => null()
       lastLoc => null()
    end do ! end do jB=1,N
    if (beta==cmplx_0) then
       ! convert format list to a sparse matrix with nnz nonzero elements
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    else
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (lenGlb==0) then
          call psp_list_create_mat(listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       else
          call psp_list_combine_listMat(M,N,cmplx_1,listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       end if
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (allocated(row_ind_C)) deallocate(row_ind_C)
       if (allocated(col_ptr_C)) deallocate(col_ptr_C)
       if (allocated(val_C)) deallocate(val_C)
       allocate(col_ptr_C(N+1))
       allocate(row_ind_C(lenGlb))
       allocate(val_C(lenGlb))
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    end if

    deallocate(val_loc)
    deallocate(row_ind_loc)
    deallocate(col_ind_loc)
    !call list_destroy(listGlb) ! why double free?

  end subroutine psp_sst_zgespmspm_nn

  subroutine psp_sst_dgespmspm_nt(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
       row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    ! sequential sparse matrix sparse matrix multiplication ('csc' format)
    ! C(1:M,1:N) = alpha*A(1:M,1:K)*B(1:K,1:N)+beta*C(1:M,1:N)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K
    integer, intent(in) :: row_ind_A(:), col_ptr_A(:)
    integer, intent(in) :: row_ind_B(:), col_ptr_B(:)
    real(dp), intent(in) :: val_A(:), val_B(:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: row_ind_C(:), col_ptr_C(:)
    real(dp), allocatable , intent(inout) :: val_C(:)

    !**** LOCAL ***********************************!
    integer :: iB, jB, cnt, cntA, cntAdded, crtGlb, trB, ot
    integer, allocatable :: crtLoc(:), lenLoc(:)
    integer :: lenGlb
    real(dp) :: vB
    integer, allocatable :: row_ind_loc(:), col_ind_loc(:)
    real(dp), allocatable :: val_loc(:)
    type(dList), pointer :: listGlb, listLoc, lastGlb
    type(dListPtrArray), dimension(N) :: lastLoc
    type(dListPtrArray), dimension(N) :: listArrayC
    !**********************************************!

    call psp_process_opM(opB,trB)
    ! operation table
    if (trB==1) then
       ot=1 ! 'c'
    else if (trB==2) then
       ot=2 ! 't'
    else
       call die('invalid implementation')
    end if

    allocate(val_loc(M))
    allocate(row_ind_loc(M))
    allocate(col_ind_loc(M))
    allocate(crtLoc(N))
    crtLoc=0
    allocate(lenLoc(N))
    lenLoc=0


    do jB=1,K ! loop over each column in B
       ! create a linked list
       do cnt=col_ptr_B(jB),col_ptr_B(jB+1)-1
          ! get the nonzero entry (iB,jB,vB) in B
          iB=row_ind_B(cnt)
          vB=val_B(cnt)
          ! get the jB'th column in A, the index range is col_ptr_A(jB):col_ptr_A(jB+1)-1
          cntAdded=0
          do cntA=col_ptr_A(jB),col_ptr_A(jB+1)-1
             cntAdded = cntAdded+1
             val_loc(cntAdded) = val_A(cntA)*vB*alpha
             row_ind_loc(cntAdded) = row_ind_A(cntA)
             col_ind_loc(cntAdded) = iB ! add to the iB'th column in C
          end do
          if (cntAdded>0) then ! add to the iB'th column in C
             crtLoc(iB)=crtLoc(iB)+1
             ! add the sparse matrix (row_ind_loc,col_ind_loc,val_loc) in the iB'th linked list
             if (crtLoc(iB)==1) then
                call psp_list_create_mat(listArrayC(iB)%ptr,lenLoc(iB),1.0_dp,row_ind_loc,col_ind_loc,&
                     val_loc,'coo',cntAdded,lastLoc(iB)%ptr)
             else
                call psp_list_combine_listMat(M,N,1.0_dp,listArrayC(iB)%ptr,lenLoc(iB),1.0_dp,row_ind_loc,&
                     col_ind_loc,val_loc,'coo',cntAdded,lastLoc(iB)%ptr)
             end if
          end if
       end do
    end do ! end do jB=1,K

    if (.false.) then
       ! bug in the first method
       crtGlb=0
       do cnt=1,N
          if (crtLoc(cnt)>0) then
             !call psp_list_print('ch1',listArrayC(cnt)%ptr, cnt)
             crtGlb=crtGlb+1
             if (crtGlb==1) then
                ! add the matrix stored in the linked list listLoc(cnt) into listGlb, correspondingly to C
                listGlb=>listArrayC(cnt)%ptr
                lastGlb=>lastLoc(cnt)%ptr
                lenGlb=lenLoc(cnt)
             else
                ! add the matrix stored in the linked list listLoc(cnt) into listGlb, correspondingly to C
                lastGlb%next => listArrayC(cnt)%ptr
                lastGlb => lastLoc(cnt)%ptr
                lenGlb=lenGlb+lenLoc(cnt)
             end if
          end if
       end do
    else
       !second method
       lenGlb=0
       cnt=1
       do while (cnt<=N)
          if (crtLoc(cnt)>0) then
             listGlb=>listArrayC(cnt)%ptr
             lastGlb=>lastLoc(cnt)%ptr
             lenGlb=lenLoc(cnt)
             exit
          else
             cnt=cnt+1
          end if
       end do
       cntA=cnt
       if (N>1) then
          do cnt=cntA,N-1
             lenGlb=lenGlb+lenLoc(cnt+1)
             if (crtLoc(cnt)>0) then
                if (crtLoc(cnt+1)>0) then
                   lastLoc(cnt)%ptr%next=>listArrayC(cnt+1)%ptr
                else
                   lastLoc(cnt+1)%ptr=>lastLoc(cnt)%ptr
                   listArrayC(cnt+1)%ptr=>listArrayC(cnt)%ptr
                end if
             else
                if (crtLoc(cnt+1)>0) then
                   lastLoc(cnt)%ptr%next=>listArrayC(cnt+1)%ptr
                else
                   lastLoc(cnt+1)%ptr=>lastLoc(cnt)%ptr
                   listArrayC(cnt+1)%ptr=>listArrayC(cnt)%ptr
                end if
             end if
          end do
       end if
    end if

    if (beta==0.0_dp) then
       ! convert format list to a sparse matrix with nnz nonzero elements
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    else
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (lenGlb==0) then
          call psp_list_create_mat(listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       else
          call psp_list_combine_listMat(M,N,1.0_dp,listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       end if
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (allocated(row_ind_C)) deallocate(row_ind_C)
       if (allocated(col_ptr_C)) deallocate(col_ptr_C)
       if (allocated(val_C)) deallocate(val_C)
       allocate(col_ptr_C(N+1))
       allocate(row_ind_C(lenGlb))
       allocate(val_C(lenGlb))
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    end if

    deallocate(val_loc)
    deallocate(row_ind_loc)
    deallocate(col_ind_loc)
    deallocate(lenLoc)
    deallocate(crtLoc)
    !call list_destroy(listGlb) ! why double free?

  end subroutine psp_sst_dgespmspm_nt

  subroutine psp_sst_zgespmspm_nt(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
       row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    ! sequential sparse matrix sparse matrix multiplication ('csc' format)
    ! C(1:M,1:N) = alpha*opA(A)(1:M,1:K)*opB(B)(1:K,1:N)+beta*C(1:M,1:N)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K
    integer, intent(in) :: row_ind_A(:), col_ptr_A(:)
    integer, intent(in) :: row_ind_B(:), col_ptr_B(:)
    complex(dp), intent(in) :: val_A(:), val_B(:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: row_ind_C(:), col_ptr_C(:)
    complex(dp), allocatable , intent(inout) :: val_C(:)

    !**** LOCAL ***********************************!
    integer :: iB, jB, cnt, cntA, cntAdded, crtGlb, trB, ot
    integer, allocatable :: crtLoc(:), lenLoc(:)
    integer :: lenGlb
    complex(dp) :: vB
    integer, allocatable :: row_ind_loc(:), col_ind_loc(:)
    complex(dp), allocatable :: val_loc(:)
    type(zList), pointer :: listGlb, listLoc, lastGlb
    type(zListPtrArray), dimension(N) :: lastLoc
    type(zListPtrArray), dimension(N) :: listArrayC
    !**********************************************!

    call psp_process_opM(opB,trB)
    ! operation table
    if (trB==1) then
       ot=1 ! 'c'
    else if (trB==2) then
       ot=2 ! 't'
    else
       call die('invalid implementation')
    end if

    allocate(val_loc(M))
    allocate(row_ind_loc(M))
    allocate(col_ind_loc(M))
    allocate(crtLoc(N))
    crtLoc=0
    allocate(lenLoc(N))
    lenLoc=0
    !allocate(listArrayC(n)) ! each element points to a list corresponding to a column in C
    !allocate(lastLoc(n))


    do jB=1,K ! loop over each column in B
       ! create a linked list
       do cnt=col_ptr_B(jB),col_ptr_B(jB+1)-1
          ! get the nonzero entry (iB,jB,vB) in B
          iB=row_ind_B(cnt)
          if (ot==1) then
             vB=CONJG(val_B(cnt))
          else
             vB=val_B(cnt)
          end if
          ! get the jB'th column in A, the index range is col_ptr_A(jB):col_ptr_A(jB+1)-1
          cntAdded=0
          do cntA=col_ptr_A(jB),col_ptr_A(jB+1)-1
             cntAdded = cntAdded+1
             val_loc(cntAdded) = val_A(cntA)*vB*alpha
             row_ind_loc(cntAdded) = row_ind_A(cntA)
             col_ind_loc(cntAdded) = iB ! add to the iB'th column in C
          end do
          if (cntAdded>0) then ! add to the iB'th column in C
             crtLoc(iB)=crtLoc(iB)+1
             ! add the sparse matrix (row_ind_loc,col_ind_loc,val_loc) in the iB'th linked list
             if (crtLoc(iB)==1) then
                call psp_list_create_mat(listArrayC(iB)%ptr,lenLoc(iB),cmplx_1,row_ind_loc,col_ind_loc,&
                     val_loc,'coo',cntAdded,lastLoc(iB)%ptr)
             else
                call psp_list_combine_listMat(M,N,cmplx_1,listArrayC(iB)%ptr,lenLoc(iB),cmplx_1,row_ind_loc,&
                     col_ind_loc,val_loc,'coo',cntAdded,lastLoc(iB)%ptr)
             end if
          end if
       end do
    end do ! end do jB=1,K

    if (.false.) then
       ! bug in the first method
       crtGlb=0
       do cnt=1,N
          if (crtLoc(cnt)>0) then
             !call psp_list_print('ch1',listArrayC(cnt)%ptr, cnt)
             crtGlb=crtGlb+1
             if (crtGlb==1) then
                ! add the matrix stored in the linked list listLoc(cnt) into listGlb, correspondingly to C
                listGlb=>listArrayC(cnt)%ptr
                lastGlb=>lastLoc(cnt)%ptr
                lenGlb=lenLoc(cnt)
             else
                ! add the matrix stored in the linked list listLoc(cnt) into listGlb, correspondingly to C
                lastGlb%next => listArrayC(cnt)%ptr
                lastGlb => lastLoc(cnt)%ptr
                lenGlb=lenGlb+lenLoc(cnt)
             end if
          end if
       end do
    else
       !second method
       lenGlb=0
       cnt=1
       do while (cnt<=N)
          if (crtLoc(cnt)>0) then
             listGlb=>listArrayC(cnt)%ptr
             lastGlb=>lastLoc(cnt)%ptr
             lenGlb=lenLoc(cnt)
             exit
          else
             cnt=cnt+1
          end if
       end do
       cntA=cnt
       if (N>1) then
          do cnt=cntA,N-1
             lenGlb=lenGlb+lenLoc(cnt+1)
             if (crtLoc(cnt)>0) then
                if (crtLoc(cnt+1)>0) then
                   lastLoc(cnt)%ptr%next=>listArrayC(cnt+1)%ptr
                else
                   lastLoc(cnt+1)%ptr=>lastLoc(cnt)%ptr
                   listArrayC(cnt+1)%ptr=>listArrayC(cnt)%ptr
                end if
             else
                if (crtLoc(cnt+1)>0) then
                   lastLoc(cnt)%ptr%next=>listArrayC(cnt+1)%ptr
                else
                   lastLoc(cnt+1)%ptr=>lastLoc(cnt)%ptr
                   listArrayC(cnt+1)%ptr=>listArrayC(cnt)%ptr
                end if
             end if
          end do
       end if
    end if

    if (beta==cmplx_0) then
       ! convert format list to a sparse matrix with nnz nonzero elements
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    else
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (lenGlb==0) then
          call psp_list_create_mat(listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       else
          call psp_list_combine_listMat(M,N,cmplx_1,listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       end if
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (allocated(row_ind_C)) deallocate(row_ind_C)
       if (allocated(col_ptr_C)) deallocate(col_ptr_C)
       if (allocated(val_C)) deallocate(val_C)
       allocate(col_ptr_C(N+1))
       allocate(row_ind_C(lenGlb))
       allocate(val_C(lenGlb))
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    end if

    deallocate(val_loc)
    deallocate(row_ind_loc)
    deallocate(col_ind_loc)
    deallocate(lenLoc)
    deallocate(crtLoc)
    !call list_destroy(listGlb) ! why double free?

  end subroutine psp_sst_zgespmspm_nt

  subroutine psp_sst_dgespmspm_tn(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
       row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    ! sequential sparse matrix sparse matrix multiplication ('csc' format)
    ! C(1:M,1:N) = alpha*opA(A)(1:M,1:K)*B(1:K,1:N)+beta*C(1:M,1:N)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K
    integer, intent(in) :: row_ind_A(:), col_ptr_A(:)
    integer, intent(in) :: row_ind_B(:), col_ptr_B(:)
    real(dp), intent(in) :: val_A(:), val_B(:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: row_ind_C(:), col_ptr_C(:)
    real(dp), allocatable , intent(inout) :: val_C(:)

    !**** LOCAL ***********************************!
    integer :: iB, jB, jA, idx_col, cntA, cntB, trA, ot
    integer :: lenGlb, numB, numA
    integer, allocatable :: idxA(:), idxB(:)
    real(dp), allocatable :: valA(:), valB(:)
    real(dp) :: sol
    type(dList), pointer :: listGlb, lastGlb
    type(dNodeData) :: elem

    !**********************************************!

    allocate(valA(K))
    allocate(valB(K))
    allocate(idxA(K))
    allocate(idxB(K))

    call psp_process_opM(opA,trA)
    ! operation table
    if (trA==1) then
       ot=1 ! 'c'
    else if (trA==2) then
       ot=2 ! 't'
    else
       call die('invalid implementation')
    end if

    lenGlb=0
    do jB=1,N
       if (col_ptr_B(jB+1)/=col_ptr_B(jB)) then ! the jB'th column of B is not empty
          ! get jB'th column in B
          numB=0
          do cntB=col_ptr_B(jB),col_ptr_B(jB+1)-1
             numB=numB+1
             valB(numB)=val_B(cntB)
             idxB(numB)=row_ind_B(cntB)
          end do
          do jA=1,M
             if (col_ptr_A(jA+1)/=col_ptr_A(jA)) then ! the jA'th column of A is not empty
                ! get jA'th column in A
                numA=0
                do cntA=col_ptr_A(jA),col_ptr_A(jA+1)-1
                   numA=numA+1
                   valA(numA)=val_A(cntA)
                   idxA(numA)=row_ind_A(cntA)
                end do
                ! compute inner product of two sparse vectors
                if (ot==1) then
                   call psp_sst_DOTC(K,idxA,valA,idxB,valB,sol,numA,numB)
                else
                   call psp_sst_DOT(K,idxA,valA,idxB,valB,sol,numA,numB)
                end if
                ! put (iC,jC,vC)=(jA,jB,sol) in C
                lenGlb=lenGlb+1
                elem%row_ind=jA
                elem%col_ind=jB
                elem%val=sol*alpha
                if (lenGlb==1) then
                   call list_create(listGlb,elem)
                   lastGlb=>listGlb
                else
                   call list_insert(lastGlb,elem)
                   lastGlb=>list_next(lastGlb)
                end if
             end if
          end do
       end if
    end do

    ! list format to ST format for the sparse matrix C
    if (beta==0.0_dp) then
       ! convert format list to a sparse matrix with nnz nonzero elements
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    else
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (lenGlb==0) then
          call psp_list_create_mat(listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       else
          call psp_list_combine_listMat(M,N,1.0_dp,listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       end if
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (allocated(row_ind_C)) deallocate(row_ind_C)
       if (allocated(col_ptr_C)) deallocate(col_ptr_C)
       if (allocated(val_C)) deallocate(val_C)
       allocate(col_ptr_C(N+1))
       allocate(row_ind_C(lenGlb))
       allocate(val_C(lenGlb))
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    end if

    deallocate(valA)
    deallocate(valB)
    deallocate(idxA)
    deallocate(idxB)
    !call list_destroy(listGlb) ! why double free?

  end subroutine psp_sst_dgespmspm_tn

  subroutine psp_sst_zgespmspm_tn(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
       row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    ! sequential sparse matrix sparse matrix multiplication ('csc' format)
    ! C(1:M,1:N) = alpha*opA(A)(1:M,1:K)*B(1:K,1:N)+beta*C(1:M,1:N)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K
    integer, intent(in) :: row_ind_A(:), col_ptr_A(:)
    integer, intent(in) :: row_ind_B(:), col_ptr_B(:)
    complex(dp), intent(in) :: val_A(:), val_B(:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: row_ind_C(:), col_ptr_C(:)
    complex(dp), allocatable , intent(inout) :: val_C(:)

    !**** LOCAL ***********************************!
    integer :: iB, jB, jA, idx_col, cntA, cntB, trA, ot
    integer :: lenGlb, numB, numA
    integer, allocatable :: idxA(:), idxB(:)
    complex(dp), allocatable :: valA(:), valB(:)
    complex(dp) :: sol
    type(zList), pointer :: listGlb, lastGlb
    type(zNodeData) :: elem

    !**********************************************!

    allocate(valA(K))
    allocate(valB(K))
    allocate(idxA(K))
    allocate(idxB(K))

    call psp_process_opM(opA,trA)
    ! operation table
    if (trA==1) then
       ot=1 ! 'c'
    else if (trA==2) then
       ot=2 ! 't'
    else
       call die('invalid implementation')
    end if

    lenGlb=0
    do jB=1,N
       if (col_ptr_B(jB+1)/=col_ptr_B(jB)) then ! the jB'th column of B is not empty
          ! get jB'th column in B
          numB=0
          do cntB=col_ptr_B(jB),col_ptr_B(jB+1)-1
             numB=numB+1
             valB(numB)=val_B(cntB)
             idxB(numB)=row_ind_B(cntB)
          end do
          do jA=1,M
             if (col_ptr_A(jA+1)/=col_ptr_A(jA)) then ! the jA'th column of A is not empty
                ! get jA'th column in A
                numA=0
                do cntA=col_ptr_A(jA),col_ptr_A(jA+1)-1
                   numA=numA+1
                   valA(numA)=val_A(cntA)
                   idxA(numA)=row_ind_A(cntA)
                end do
                ! compute inner product of two sparse vectors
                if (ot==1) then
                   call psp_sst_DOTC(K,idxA,valA,idxB,valB,sol,numA,numB)
                else
                   call psp_sst_DOT(K,idxA,valA,idxB,valB,sol,numA,numB)
                end if
                ! put (iC,jC,vC)=(jA,jB,sol) in C
                lenGlb=lenGlb+1
                elem%row_ind=jA
                elem%col_ind=jB
                elem%val=sol*alpha
                if (lenGlb==1) then
                   call list_create(listGlb,elem)
                   lastGlb=>listGlb
                else
                   call list_insert(lastGlb,elem)
                   lastGlb=>list_next(lastGlb)
                end if
             end if
          end do
       end if
    end do

    ! list format to ST format for the sparse matrix C
    if (beta==cmplx_0) then
       ! convert format list to a sparse matrix with nnz nonzero elements
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    else
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (lenGlb==0) then
          call psp_list_create_mat(listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       else
          call psp_list_combine_listMat(M,N,cmplx_1,listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       end if
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (allocated(row_ind_C)) deallocate(row_ind_C)
       if (allocated(col_ptr_C)) deallocate(col_ptr_C)
       if (allocated(val_C)) deallocate(val_C)
       allocate(col_ptr_C(N+1))
       allocate(row_ind_C(lenGlb))
       allocate(val_C(lenGlb))
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    end if

    deallocate(valA)
    deallocate(valB)
    deallocate(idxA)
    deallocate(idxB)
    !call list_destroy(listGlb) ! why double free?

  end subroutine psp_sst_zgespmspm_tn

  subroutine psp_sst_dgespmspm_tt(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
       row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    ! sequential sparse matrix sparse matrix multiplication ('csc' format)
    ! C(1:M,1:N) = alpha*opA(A)(1:M,1:K)*opB(B)(1:K,1:N)+beta*C(1:M,1:N)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K
    integer, intent(in) :: row_ind_A(:), col_ptr_A(:)
    integer, intent(in) :: row_ind_B(:), col_ptr_B(:)
    real(dp), intent(in) :: val_A(:), val_B(:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: row_ind_C(:), col_ptr_C(:)
    real(dp), allocatable , intent(inout) :: val_C(:)

    !**** LOCAL ***********************************!
    integer :: iB, jB, jA, idx_col, cntA, cntB, trA, trB, otA, otB
    integer :: lenGlb, numB, numA
    integer, allocatable :: idxA(:), idxB(:), lenLoc(:)
    real(dp), allocatable :: valA(:), valB(:)
    real(dp) :: sol
    type(dList), pointer :: listGlb, lastGlb
    type(dNodeData) :: elem
    type(dListPtrArray), dimension(N) :: crtLoc
    type(dListPtrArray), dimension(N) :: listArrayB

    !**********************************************!

    allocate(valA(K))
    allocate(valB(K))
    allocate(idxA(K))
    allocate(idxB(K))
    allocate(lenLoc(N))

    call psp_process_opM(opA,trA)
    ! operation table
    if (trA==1) then
       otA=1 ! 'c'
    else if (trA==2) then
       otA=2 ! 't'
    else
       call die('invalid implementation')
    end if

    call psp_process_opM(opB,trB)
    ! operation table
    if (trB==1) then
       otB=1 ! 'c'
    else if (trB==2) then
       otB=2 ! 't'
    else
       call die('invalid implementation')
    end if

    ! compute opB(B) in a format of an array of lists
    ! each list corresponds to a column in opB(B)
    lenLoc=0

    do jB=1,K
       do cntB=col_ptr_B(jB),col_ptr_B(jB+1)-1
          ! the entry (row_ind_B(cntB),jB,val_B(cntB)) in B is 
          ! the entry (jB,row_ind_B(cntB),val_B(cntB)) in opB(B)
          lenLoc(row_ind_B(cntB))=lenLoc(row_ind_B(cntB))+1
          elem%row_ind=jB
          elem%col_ind=row_ind_B(cntB)
          elem%val=val_B(cntB)
          if (lenLoc(row_ind_B(cntB))==1) then
             call list_create(listArrayB(row_ind_B(cntB))%ptr,elem)
             crtLoc(row_ind_B(cntB))%ptr=>listArrayB(row_ind_B(cntB))%ptr
          else
             call list_insert(crtLoc(row_ind_B(cntB))%ptr,elem)
             crtLoc(row_ind_B(cntB))%ptr=>list_next(crtLoc(row_ind_B(cntB))%ptr)
          end if
       end do
    end do

    lenGlb=0
    do jB=1,N ! loop over columns of opB(B)
       numB=lenLoc(jB)
       if (numB/=0) then ! the jB'th column of opB(B) is not empty
          ! get jB'th column in opB(B)
          crtLoc(jB)%ptr=>listArrayB(jB)%ptr
          do cntB=1,lenLoc(jB)
             elem=list_get_data(crtLoc(jB)%ptr)
             crtLoc(jB)%ptr=>list_next(crtLoc(jB)%ptr)
             valB(cntB)=elem%val
             idxB(cntB)=elem%row_ind
          end do
          do jA=1,M
             if (col_ptr_A(jA+1)/=col_ptr_A(jA)) then ! the jA'th column of A is not empty
                ! get jA'th column in A
                numA=0
                do cntA=col_ptr_A(jA),col_ptr_A(jA+1)-1
                   numA=numA+1
                   valA(numA)=val_A(cntA)
                   idxA(numA)=row_ind_A(cntA)
                end do
                ! compute inner product of two sparse vectors
                if (otA==1) then
                   call psp_sst_DOTC(K,idxA,valA,idxB,valB,sol,numA,numB)
                else
                   call psp_sst_DOT(K,idxA,valA,idxB,valB,sol,numA,numB)
                end if
                ! put (iC,jC,vC)=(jA,jB,sol) in C
                lenGlb=lenGlb+1
                elem%row_ind=jA
                elem%col_ind=jB
                elem%val=sol*alpha
                if (lenGlb==1) then
                   call list_create(listGlb,elem)
                   lastGlb=>listGlb
                else
                   call list_insert(lastGlb,elem)
                   lastGlb=>list_next(lastGlb)
                end if
             end if
          end do
       end if
    end do

    ! list format to ST format for the sparse matrix C
    if (beta==0.0_dp) then
       ! convert format list to a sparse matrix with nnz nonzero elements
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    else
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (lenGlb==0) then
          call psp_list_create_mat(listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       else
          call psp_list_combine_listMat(M,N,1.0_dp,listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       end if
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (allocated(row_ind_C)) deallocate(row_ind_C)
       if (allocated(col_ptr_C)) deallocate(col_ptr_C)
       if (allocated(val_C)) deallocate(val_C)
       allocate(col_ptr_C(N+1))
       allocate(row_ind_C(lenGlb))
       allocate(val_C(lenGlb))
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    end if

    deallocate(valA)
    deallocate(valB)
    deallocate(idxA)
    deallocate(idxB)
    deallocate(lenLoc)
    !call list_destroy(listGlb) ! why double free?

  end subroutine psp_sst_dgespmspm_tt

  subroutine psp_sst_zgespmspm_tt(M,N,K,opA,opB,alpha,row_ind_A,col_ptr_A,val_A, &
       row_ind_B,col_ptr_B,val_B,row_ind_C,col_ptr_C,val_C,beta)
    ! sequential sparse matrix sparse matrix multiplication ('csc' format)
    ! C(1:M,1:N) = alpha*opA(A)(1:M,1:K)*opB(B)(1:K,1:N)+beta*C(1:M,1:N)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: M, N, K
    integer, intent(in) :: row_ind_A(:), col_ptr_A(:)
    integer, intent(in) :: row_ind_B(:), col_ptr_B(:)
    complex(dp), intent(in) :: val_A(:), val_B(:), alpha, beta
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)

    !**** INOUT ***********************************!

    integer, allocatable , intent(inout) :: row_ind_C(:), col_ptr_C(:)
    complex(dp), allocatable , intent(inout) :: val_C(:)

    !**** LOCAL ***********************************!
    integer :: iB, jB, jA, idx_col, cntA, cntB, trA, trB, otA, otB
    integer :: lenGlb, numB, numA
    integer, allocatable :: idxA(:), idxB(:), lenLoc(:)
    complex(dp), allocatable :: valA(:), valB(:)
    complex(dp) :: sol
    type(zList), pointer :: listGlb, lastGlb
    type(zNodeData) :: elem
    type(zListPtrArray), dimension(N) :: crtLoc
    type(zListPtrArray), dimension(N) :: listArrayB

    !**********************************************!

    allocate(valA(K))
    allocate(valB(K))
    allocate(idxA(K))
    allocate(idxB(K))
    allocate(lenLoc(N))

    call psp_process_opM(opA,trA)
    ! operation table
    if (trA==1) then
       otA=1 ! 'c'
    else if (trA==2) then
       otA=2 ! 't'
    else
       call die('invalid implementation')
    end if

    call psp_process_opM(opB,trB)
    ! operation table
    if (trB==1) then
       otB=1 ! 'c'
    else if (trB==2) then
       otB=2 ! 't'
    else
       call die('invalid implementation')
    end if

    ! compute opB(B) in a format of an array of lists
    ! each list corresponds to a column in opB(B)
    lenLoc=0
    do jB=1,K
       do cntB=col_ptr_B(jB),col_ptr_B(jB+1)-1
          ! the entry (row_ind_B(cntB),jB,val_B(cntB)) in B is
          ! the entry (jB,row_ind_B(cntB),val_B(cntB)) in opB(B)
          lenLoc(row_ind_B(cntB))=lenLoc(row_ind_B(cntB))+1
          elem%row_ind=jB
          elem%col_ind=row_ind_B(cntB)
          if (otB==1) then
             elem%val=CONJG(val_B(cntB))
          else
             elem%val=val_B(cntB)
          end if
          if (lenLoc(row_ind_B(cntB))==1) then
             call list_create(listArrayB(row_ind_B(cntB))%ptr,elem)
             crtLoc(row_ind_B(cntB))%ptr=>listArrayB(row_ind_B(cntB))%ptr
          else
             call list_insert(crtLoc(row_ind_B(cntB))%ptr,elem)
             crtLoc(row_ind_B(cntB))%ptr=>list_next(crtLoc(row_ind_B(cntB))%ptr)
          end if
       end do
    end do

    lenGlb=0
    do jB=1,N ! loop over columns of opB(B)
       numB=lenLoc(jB)
       if (numB/=0) then ! the jB'th column of opB(B) is not empty
          ! get jB'th column in opB(B)
          crtLoc(jB)%ptr=>listArrayB(jB)%ptr
          do cntB=1,lenLoc(jB)
             elem=list_get_data(crtLoc(jB)%ptr)
             crtLoc(jB)%ptr=>list_next(crtLoc(jB)%ptr)
             valB(cntB)=elem%val
             idxB(cntB)=elem%row_ind
          end do
          do jA=1,M
             if (col_ptr_A(jA+1)/=col_ptr_A(jA)) then ! the jA'th column of A is not empty
                ! get jA'th column in A
                numA=0
                do cntA=col_ptr_A(jA),col_ptr_A(jA+1)-1
                   numA=numA+1
                   valA(numA)=val_A(cntA)
                   idxA(numA)=row_ind_A(cntA)
                end do
                ! compute inner product of two sparse vectors
                if (otA==1) then
                   call psp_sst_DOTC(K,idxA,valA,idxB,valB,sol,numA,numB)
                else
                   call psp_sst_DOT(K,idxA,valA,idxB,valB,sol,numA,numB)
                end if
                ! put (iC,jC,vC)=(jA,jB,sol) in C
                lenGlb=lenGlb+1
                elem%row_ind=jA
                elem%col_ind=jB
                elem%val=sol*alpha
                if (lenGlb==1) then
                   call list_create(listGlb,elem)
                   lastGlb=>listGlb
                else
                   call list_insert(lastGlb,elem)
                   lastGlb=>list_next(lastGlb)
                end if
             end if
          end do
       end if
    end do

    ! list format to ST format for the sparse matrix C
    if (beta==cmplx_0) then
       ! convert format list to a sparse matrix with nnz nonzero elements
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    else
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (lenGlb==0) then
          call psp_list_create_mat(listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       else
          call psp_list_combine_listMat(M,N,cmplx_1,listGlb,lenGlb,beta,row_ind_C,col_ptr_C,val_C,'csc',size(val_C),lastGlb)
       end if
       ! convert format list to a sparse matrix with nnz nonzero elements
       if (allocated(row_ind_C)) deallocate(row_ind_C)
       if (allocated(col_ptr_C)) deallocate(col_ptr_C)
       if (allocated(val_C)) deallocate(val_C)
       allocate(col_ptr_C(N+1))
       allocate(row_ind_C(lenGlb))
       allocate(val_C(lenGlb))
       call psp_list2spm(M,N,row_ind_C,col_ptr_C,val_C,'csc',listGlb,lenGlb,.false.)
    end if

    deallocate(valA)
    deallocate(valB)
    deallocate(idxA)
    deallocate(idxB)
    deallocate(lenLoc)
    !call list_destroy(listGlb) ! why double free?

  end subroutine psp_sst_zgespmspm_tt



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


END MODULE psp_spBLAS_Level3
