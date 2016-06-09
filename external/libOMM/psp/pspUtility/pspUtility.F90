MODULE pspUtility
  use pspVariable

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

  interface psp_sst_gespmm
     module procedure psp_sst_dgespmm
     module procedure psp_sst_zgespmm
  end interface psp_sst_gespmm

  interface psp_sst_gemspm
     module procedure psp_sst_dgemspm
     module procedure psp_sst_zgemspm
  end interface psp_sst_gemspm

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

  public :: psp_coo2csc
  public :: psp_csc2coo
  public :: psp_copy_spm2st
  public :: psp_copy_m
  public :: psp_copy_v
  public :: psp_idx_glb2loc
  public :: psp_sst_gespmm
  public :: psp_sst_gemspm
  public :: psp_process_opM
  public :: init_random_seed

contains

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

  end subroutine psp_copy_dspm2st

  subroutine psp_copy_zspm2st(M,N,A,IA,JA,B_idx1,B_idx2,B_val,B_dim1,B_dim2,IB,JB,beta)
    ! B(IB:IB+M-1,JB:JB+N-1) = A(IA:IA+M-1,JA:JA+N-1)
    ! other entries of B are set to be zero
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

  end subroutine psp_copy_zspm2st

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


END MODULE pspUtility
