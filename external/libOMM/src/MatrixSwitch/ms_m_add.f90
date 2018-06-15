module ms_m_add
  use MatrixSwitch_ops

  implicit none

contains

  !================================================!
  ! implementation: reference                      !
  !================================================!

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

end module ms_m_add
