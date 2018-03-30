module ms_m_register
  use MatrixSwitch_ops

  implicit none

  !**** INTERFACES ********************************!

  interface m_register_sden
     module procedure m_register_sdden
     module procedure m_register_szden
  end interface m_register_sden

  interface m_register_pdbc
     module procedure m_register_pddbc
     module procedure m_register_pzdbc
  end interface m_register_pdbc

  !************************************************!

contains

  !================================================!
  ! register matrix: simple dense serial           !
  !================================================!
  subroutine m_register_sdden(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in), target :: A(:,:) ! two-dimensional array containing the matrix elements

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name ! matrix to be allocated

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!

    dim=shape(A)
    m_name%dim1=dim(1)
    m_name%dim2=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='den'
    m_name%is_serial=.true.
    m_name%is_real=.true.
    m_name%is_sparse=.false.

    m_name%dval => A

    m_name%is_initialized=.true.

  end subroutine m_register_sdden

  subroutine m_register_szden(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    complex(dp), intent(in), target :: A(:,:) ! two-dimensional array containing the matrix elements

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name ! matrix to be allocated

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!

    dim=shape(A)
    m_name%dim1=dim(1)
    m_name%dim2=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='den'
    m_name%is_serial=.true.
    m_name%is_real=.false.
    m_name%is_sparse=.false.

    m_name%zval => A

    m_name%is_initialized=.true.

  end subroutine m_register_szden

  !================================================!
  ! register matrix: dense block cyclic parallel   !
  !================================================!
  subroutine m_register_pddbc(m_name,A,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9) ! BLACS array descriptor

    real(dp), intent(in), target :: A(:,:) ! two-dimensional array containing the local matrix elements

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name ! matrix to be allocated

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!

    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    dim=shape(A)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='dbc'
    m_name%is_serial=.false.
    m_name%is_real=.true.
    m_name%is_sparse=.false.

    m_name%dval => A

    m_name%is_initialized=.true.

  end subroutine m_register_pddbc

  subroutine m_register_pzdbc(m_name,A,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: desc(9) ! BLACS array descriptor

    complex(dp), intent(in), target :: A(:,:) ! two-dimensional array containing the local matrix elements

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name ! matrix to be allocated

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!

    allocate(m_name%iaux1(9))
    m_name%iaux1_is_allocated=.true.
    m_name%iaux1=desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
    dim=shape(A)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type='dbc'
    m_name%is_serial=.false.
    m_name%is_real=.false.
    m_name%is_sparse=.false.

    m_name%zval => A

    m_name%is_initialized=.true.

  end subroutine m_register_pzdbc

end module ms_m_register
