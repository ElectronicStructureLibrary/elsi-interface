MODULE pspMatSum
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

  interface psp_sum_spmm
     ! compute the sum of a sparse matrix and a dense matrix
     module procedure psp_dsum_spmm
     module procedure psp_zsum_spmm
  end interface psp_sum_spmm

  interface psp_sum_spmspm
     ! compute the sum of two sparse matrices
     module procedure psp_dsum_spmspm
     module procedure psp_zsum_spmspm
  end interface psp_sum_spmspm

  interface die
     module procedure die
  end interface die

  public :: psp_sum_spmm
  public :: psp_sum_spmspm

contains

  subroutine psp_dsum_spmm(A,B,C,alpha,beta)
    ! compute the sum of a sparse matrix and a dense matrix, C=alpha*A+beta*B
    ! A is a sparse matri, B is a dense matrix, C is a dense matrix
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in)  ::   alpha, beta ! scalar
    type(psp_matrix_spm), intent(in) :: A ! sparse matrix A
    real(dp), intent(in) :: B(:,:) ! dense matrix B

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: C(:,:) ! dense matrix C

    !**********************************************!

    if (A%str_type=='coo') then
       call psp_sst_sum_spmm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ind,A%dval,'coo',beta,B,C)
    else
       call psp_sst_sum_spmm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ptr,A%dval,'csc',beta,B,C)
    end if

  end subroutine psp_dsum_spmm

  subroutine psp_zsum_spmm(A,B,C,alpha,beta)
    ! compute the sum of a sparse matrix and a dense matrix, C=alpha*A+beta*B
    ! A is a sparse matri, B is a dense matrix, C is a dense matrix
    implicit none

    !**** INPUT ***********************************!

    complex(dp), intent(in)  ::   alpha, beta ! scalar
    type(psp_matrix_spm), intent(in) :: A ! sparse matrix A
    complex(dp), intent(in) :: B(:,:) ! dense matrix B

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: C(:,:) ! dense matrix C

    !**********************************************!

    if (A%str_type=='coo') then
       call psp_sst_sum_spmm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ind,A%zval,'coo',beta,B,C)
    else
       call psp_sst_sum_spmm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ptr,A%zval,'csc',beta,B,C)
    end if

  end subroutine psp_zsum_spmm

  subroutine psp_dsum_spmspm(A,B,C,alpha,beta)
    ! compute the sum of a sparse matrix and a dense matrix, C=alpha*A+beta*B
    ! A is a sparse matri, B is a dense matrix, C is a dense matrix
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in)  ::   alpha, beta ! scalar
    type(psp_matrix_spm), intent(in) :: A, B ! sparse matrices A and B
    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout) :: C ! sparse matrix C

    !**********************************************!

    select case (A%str_type)
    case ('coo')
       select case (B%str_type)
       case('coo')
          select case (C%str_type)
          case ('coo')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ind,A%dval,'coo',beta,&
                  B%row_ind,B%col_ind,B%dval,'coo',C%row_ind,C%col_ind,C%dval,'coo',size(A%row_ind),size(B%row_ind))
          case ('csc')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ind,A%dval,'coo',beta,&
                  B%row_ind,B%col_ind,B%dval,'coo',C%row_ind,C%col_ptr,C%dval,'csc',size(A%row_ind),size(B%row_ind))
          end select
       case('csc')
          select case (C%str_type)
          case ('coo')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ind,A%dval,'coo',beta,&
                  B%row_ind,B%col_ptr,B%dval,'csc',C%row_ind,C%col_ind,C%dval,'coo',size(A%row_ind),size(B%row_ind))
          case ('csc')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ind,A%dval,'coo',beta,&
                  B%row_ind,B%col_ptr,B%dval,'csc',C%row_ind,C%col_ptr,C%dval,'csc',size(A%row_ind),size(B%row_ind))
          end select
       end select
    case('csc')
       select case (B%str_type)
       case('coo')
          select case (C%str_type)
          case ('coo')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ptr,A%dval,'csc',beta,&
                  B%row_ind,B%col_ind,B%dval,'coo',C%row_ind,C%col_ind,C%dval,'coo',size(A%row_ind),size(B%row_ind))
          case ('csc')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ptr,A%dval,'csc',beta,&
                  B%row_ind,B%col_ind,B%dval,'coo',C%row_ind,C%col_ptr,C%dval,'csc',size(A%row_ind),size(B%row_ind))
          end select
       case('csc')
          select case (C%str_type)
          case ('coo')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ptr,A%dval,'csc',beta,&
                  B%row_ind,B%col_ptr,B%dval,'csc',C%row_ind,C%col_ind,C%dval,'coo',size(A%row_ind),size(B%row_ind))
          case ('csc')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ptr,A%dval,'csc',beta,&
                  B%row_ind,B%col_ptr,B%dval,'csc',C%row_ind,C%col_ptr,C%dval,'csc',size(A%row_ind),size(B%row_ind))
          end select
       end select
    end select
    C%nnz=size(C%dval)

  end subroutine psp_dsum_spmspm

  subroutine psp_zsum_spmspm(A,B,C,alpha,beta)
    ! compute the sum of a sparse matrix and a dense matrix, C=alpha*A+beta*B
    ! A is a sparse matri, B is a dense matrix, C is a dense matrix
    implicit none

    !**** INPUT ***********************************!

    complex(dp), intent(in)  ::   alpha, beta ! scalar
    type(psp_matrix_spm), intent(in) :: A, B ! sparse matrices A and B
    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout) :: C ! sparse matrix C

    !**********************************************!

    select case (A%str_type)
    case ('coo')
       select case (B%str_type)
       case('coo')
          select case (C%str_type)
          case ('coo')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ind,A%zval,'coo',beta,&
                  B%row_ind,B%col_ind,B%zval,'coo',C%row_ind,C%col_ind,C%zval,'coo',size(A%row_ind),size(B%row_ind))
          case ('csc')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ind,A%zval,'coo',beta,&
                  B%row_ind,B%col_ind,B%zval,'coo',C%row_ind,C%col_ptr,C%zval,'csc',size(A%row_ind),size(B%row_ind))
          end select
       case('csc')
          select case (C%str_type)
          case ('coo')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ind,A%zval,'coo',beta,&
                  B%row_ind,B%col_ptr,B%zval,'csc',C%row_ind,C%col_ind,C%zval,'coo',size(A%row_ind),size(B%row_ind))
          case ('csc')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ind,A%zval,'coo',beta,&
                  B%row_ind,B%col_ptr,B%zval,'csc',C%row_ind,C%col_ptr,C%zval,'csc',size(A%row_ind),size(B%row_ind))
          end select
       end select
    case('csc')
       select case (B%str_type)
       case('coo')
          select case (C%str_type)
          case ('coo')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ptr,A%zval,'csc',beta,&
                  B%row_ind,B%col_ind,B%zval,'coo',C%row_ind,C%col_ind,C%zval,'coo',size(A%row_ind),size(B%row_ind))
          case ('csc')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ptr,A%zval,'csc',beta,&
                  B%row_ind,B%col_ind,B%zval,'coo',C%row_ind,C%col_ptr,C%zval,'csc',size(A%row_ind),size(B%row_ind))
          end select
       case('csc')
          select case (C%str_type)
          case ('coo')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ptr,A%zval,'csc',beta,&
                  B%row_ind,B%col_ptr,B%zval,'csc',C%row_ind,C%col_ind,C%zval,'coo',size(A%row_ind),size(B%row_ind))
          case ('csc')
             call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,alpha,A%row_ind,A%col_ptr,A%zval,'csc',beta,&
                  B%row_ind,B%col_ptr,B%zval,'csc',C%row_ind,C%col_ptr,C%zval,'csc',size(A%row_ind),size(B%row_ind))
          end select
       end select
    end select
    C%nnz=size(C%zval)

  end subroutine psp_zsum_spmspm

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


END MODULE pspMatSum
