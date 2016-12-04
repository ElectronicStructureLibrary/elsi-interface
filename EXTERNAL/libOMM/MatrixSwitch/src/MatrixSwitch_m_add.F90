module MatrixSwitch_m_add
  use MatrixSwitch_ops

  implicit none

contains

  !================================================!
  ! implementation: reference                      !
  !================================================!

#ifdef PSP
  subroutine m_add_pdcscpddbcref(A,C,alpha,beta)
    implicit none
    include 'mpif.h'

    !**** INPUT ***********************************!

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, l

    !**********************************************!

    C%dval=beta*C%dval

    do i=1,A%spm%loc_dim2
       do j=0,A%spm%col_ptr(i+1)-A%spm%col_ptr(i)-1
          l=A%spm%col_ptr(i)+j
          C%dval(A%spm%row_ind(l),i)=C%dval(A%spm%row_ind(l),i)+alpha*A%spm%dval(l)
       end do
    end do

  end subroutine m_add_pdcscpddbcref

  subroutine m_add_pzcscpzdbcref(A,C,alpha,beta)
    implicit none
    include 'mpif.h'

    !**** INPUT ***********************************!

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, l

    !**********************************************!

    C%zval=beta*C%zval

    do i=1,A%spm%loc_dim2
       do j=0,A%spm%col_ptr(i+1)-A%spm%col_ptr(i)-1
          l=A%spm%col_ptr(i)+j
          C%zval(A%spm%row_ind(l),i)=C%zval(A%spm%row_ind(l),i)+alpha*A%spm%zval(l)
       end do
    end do

  end subroutine m_add_pzcscpzdbcref
#endif

  subroutine m_add_sddenref(A,trA,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    logical, intent(in) :: trA

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j

    !**********************************************!

    C%dval=beta*C%dval

    if (.not. trA) then
       C%dval=C%dval+alpha*A%dval
    else if (trA) then
       do i=1,C%dim1
          do j=1,C%dim2
             C%dval(i,j)=C%dval(i,j)+alpha*A%dval(j,i)
          end do
       end do
    end if

  end subroutine m_add_sddenref

  subroutine m_add_szdenref(A,tcA,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j

    !**********************************************!

    C%zval=beta*C%zval

    if (tcA==0) then
       C%zval=C%zval+alpha*A%zval
    else if (tcA==1) then
       do i=1,C%dim1
          do j=1,C%dim2
             C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(j,i))
          end do
       end do
    else if (tcA==2) then
       do i=1,C%dim1
          do j=1,C%dim2
             C%zval(i,j)=C%zval(i,j)+alpha*A%zval(j,i)
          end do
       end do
    end if

  end subroutine m_add_szdenref

  subroutine m_add_sdcscsddenref(A,trA,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    logical, intent(in) :: trA

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, l

    !**********************************************!

    C%dval=beta*C%dval

    if (.not. trA) then
       do j=1,C%dim2
          do l=A%iaux4(j)+1,A%iaux4(j+1)
             i=A%iaux3(l)
             C%dval(i,j)=C%dval(i,j)+alpha*A%dval(l,1)
          end do
       end do
    else if (trA) then
       do i=1,C%dim1
          do l=A%iaux4(i)+1,A%iaux4(i+1)
             j=A%iaux3(l)
             C%dval(i,j)=C%dval(i,j)+alpha*A%dval(l,1)
          end do
       end do
    end if

  end subroutine m_add_sdcscsddenref

  subroutine m_add_sdcsrsddenref(A,trA,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    logical, intent(in) :: trA

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, l

    !**********************************************!

    C%dval=beta*C%dval

    if (.not. trA) then
       do i=1,C%dim1
          do l=A%iaux3(i)+1,A%iaux3(i+1)
             j=A%iaux4(l)
             C%dval(i,j)=C%dval(i,j)+alpha*A%dval(l,1)
          end do
       end do
    else if (trA) then
       do j=1,C%dim2
          do l=A%iaux3(j)+1,A%iaux3(j+1)
             i=A%iaux4(l)
             C%dval(i,j)=C%dval(i,j)+alpha*A%dval(l,1)
          end do
       end do
    end if

  end subroutine m_add_sdcsrsddenref

  subroutine m_add_szcscszdenref(A,tcA,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, l

    !**********************************************!

    C%zval=beta*C%zval

    if (tcA==0) then
       do j=1,C%dim2
          do l=A%iaux4(j)+1,A%iaux4(j+1)
             i=A%iaux3(l)
             C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)
          end do
       end do
    else if (tcA==1) then
       do i=1,C%dim1
          do l=A%iaux4(i)+1,A%iaux4(i+1)
             j=A%iaux3(l)
             C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(l,1))
          end do
       end do
    else if (tcA==2) then
       do i=1,C%dim1
          do l=A%iaux4(i)+1,A%iaux4(i+1)
             j=A%iaux3(l)
             C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)
          end do
       end do
    end if

  end subroutine m_add_szcscszdenref

  subroutine m_add_szcsrszdenref(A,tcA,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, l

    !**********************************************!

    C%zval=beta*C%zval

    if (tcA==0) then
       do i=1,C%dim1
          do l=A%iaux3(i)+1,A%iaux3(i+1)
             j=A%iaux4(l)
             C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)
          end do
       end do
    else if (tcA==1) then
       do j=1,C%dim2
          do l=A%iaux3(j)+1,A%iaux3(j+1)
             i=A%iaux4(l)
             C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(l,1))
          end do
       end do
    else if (tcA==2) then
       do j=1,C%dim2
          do l=A%iaux3(j)+1,A%iaux3(j+1)
             i=A%iaux4(l)
             C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)
          end do
       end do
    end if

  end subroutine m_add_szcsrszdenref

end module MatrixSwitch_m_add
