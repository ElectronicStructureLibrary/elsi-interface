MODULE pspVariable

  private

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** TYPES *************************************!

  ! This is the derived type that encapsulates all matrix storage possibilities and hides the details from the user.
  type psp_matrix_spm
     character(3) :: str_type ! label identifying the storage format, 'csc', 'coo'

     logical :: is_initialized=.false. ! has the matrix been initialized?
     logical :: is_serial ! is the matrix serial or parallel distributed?
     logical :: is_real ! is the matrix real or complex (both kind dp)?
     logical :: is_square ! is the matrix square?

     integer :: nnz      ! number of nonzero entries
     integer :: glb_dim1 ! (global) row dimension size of the matrix
     integer :: glb_dim2 ! (global) column dimension size of the matrix
     integer :: loc_dim1 ! (local) row dimension size of the local matrix
     integer :: loc_dim2 ! (local) column dimension size of the local matrix
     integer, allocatable  :: row_ind(:) ! array for row indices in the coo format
     integer, allocatable  :: col_ind(:) ! array for column indices in the coo format
     integer, allocatable  :: col_ptr(:) ! array for column pointers in the csc format
     integer :: desc(9)    ! (global and local input) integer array
     ! On entry, DESCB  is an integer array of dimension 9. This
     ! is the array descriptor for the sparse matrix. The entries
     ! of desc are similar to those in Scalapack

     real(dp), allocatable  :: dval(:)  ! matrix elements for a real matrix

     complex(dp), allocatable  :: zval(:)  ! matrix elements for a complex matrix
  end type psp_matrix_spm

  !**********************************************!

  integer, external :: numroc ! it is a function to compute local size

  !**** INTERFACES ********************************!

  interface psp_register_spm
     module procedure psp_register_dspm
     module procedure psp_register_zspm
  end interface psp_register_spm

  interface psp_deallocate_spm
     module procedure psp_deallocate_spm
  end interface psp_deallocate_spm

  public :: psp_matrix_spm
  public :: psp_register_spm
  public :: psp_deallocate_spm

contains

  subroutine psp_register_dspm(m_name,idx1,idx2,val,desc,fmt,loc_dims,nprow,npcol)
    implicit none

    !**** INPUT ***********************************!
    character(3), intent(in), target :: fmt ! 'coo' or 'csc'
    integer, intent(in), target :: loc_dims(2) ! local dimensions of the sparse matrix
    integer, intent(in), target :: desc(9) ! BLACS array descriptor
    integer, intent(in), target :: idx1(:) ! one-dimensional array for row indices, local
    integer, intent(in), target :: idx2(:) ! one-dimensional array for column indices in 'coo' format or column pointers in 'csc' format, local
    integer :: nprow, npcol
    real(dp), intent(in), target :: val(:) ! one-dimensional array containing the nonzero local matrix elements

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout) :: m_name ! sparse matrix to be allocated

    !**********************************************!

    m_name%desc=desc
    m_name%nnz=size(val)
    m_name%glb_dim1=desc(3)
    m_name%glb_dim2=desc(4)
    m_name%loc_dim1=loc_dims(1)
    m_name%loc_dim2=loc_dims(2)
    if (m_name%glb_dim1==m_name%glb_dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type=fmt
    m_name%is_serial=.false.
    m_name%is_real=.true.

    if (allocated(m_name%dval)) deallocate(m_name%dval)
    allocate(m_name%dval(m_name%nnz))
    m_name%dval = val

    if (allocated(m_name%row_ind)) deallocate(m_name%row_ind)
    allocate(m_name%row_ind(m_name%nnz))
    m_name%row_ind = idx1

    if (fmt.EQ.'coo') then
       if (allocated(m_name%col_ind)) deallocate(m_name%col_ind)
       allocate(m_name%col_ind(m_name%nnz))
       m_name%col_ind = idx2
    else
       if (allocated(m_name%col_ind)) deallocate(m_name%col_ind)
       allocate(m_name%col_ind(m_name%loc_dim2+1))
       m_name%col_ptr = idx2
    end if

    m_name%is_initialized=.true.

  end subroutine psp_register_dspm


  subroutine psp_register_zspm(m_name,idx1,idx2,val,desc,fmt,loc_dims,nprow,npcol)
    implicit none

    !**** INPUT ***********************************!
    character(3), intent(in), target :: fmt ! 'coo' or 'csc'
    integer, intent(in), target :: loc_dims(2) ! local dimensions of the sparse matrix
    integer, intent(in), target :: desc(9) ! BLACS array descriptor
    integer, intent(in), target :: idx1(:) ! one-dimensional array for row indices, local
    integer, intent(in), target :: idx2(:) ! one-dimensional array for column indices in 'coo' format or column pointers in 'csc' format, local
    integer :: nprow, npcol

    complex(dp), intent(in), target :: val(:) ! one-dimensional array containing the nonzero local matrix elements

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout) :: m_name ! sparse matrix to be allocated

    !**********************************************!

    m_name%desc=desc
    m_name%nnz=size(val)
    m_name%glb_dim1=desc(3)
    m_name%glb_dim2=desc(4)
    m_name%loc_dim1=loc_dims(1)
    m_name%loc_dim2=loc_dims(2)
    if (m_name%glb_dim1==m_name%glb_dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type=fmt
    m_name%is_serial=.false.
    m_name%is_real=.true.


    if (allocated(m_name%zval)) deallocate(m_name%zval)
    allocate(m_name%zval(m_name%nnz))
    m_name%zval = val

    if (allocated(m_name%row_ind)) deallocate(m_name%row_ind)
    allocate(m_name%row_ind(m_name%nnz))
    m_name%row_ind = idx1

    if (fmt.EQ.'coo') then
       if (allocated(m_name%col_ind)) deallocate(m_name%col_ind)
       allocate(m_name%col_ind(m_name%nnz))
       m_name%col_ind = idx2
    else
       if (allocated(m_name%col_ind)) deallocate(m_name%col_ind)
       allocate(m_name%col_ind(m_name%loc_dim2+1))
       m_name%col_ptr = idx2
    end if

    m_name%is_initialized=.true.

  end subroutine psp_register_zspm

  !================================================!
  ! deallocate matrix                              !
  !================================================!

  subroutine psp_deallocate_spm(m_name)
    implicit none

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout) :: m_name ! matrix to be deallocated

    !**********************************************!

    if (allocated(m_name%row_ind)) deallocate(m_name%row_ind)
    if (allocated(m_name%col_ind)) deallocate(m_name%col_ind)
    if (allocated(m_name%col_ptr)) deallocate(m_name%col_ptr)
    if (allocated(m_name%dval)) deallocate(m_name%dval)
    if (allocated(m_name%zval)) deallocate(m_name%zval)

    m_name%is_initialized=.false.

  end subroutine psp_deallocate_spm



END MODULE pspVariable
