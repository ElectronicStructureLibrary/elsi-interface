module MatrixSwitch_mm_multiply
  use MatrixSwitch_ops

  implicit none

contains

  !================================================!
  ! implementation: reference                      !
  !================================================!

  subroutine mm_multiply_sddenref(A,trA,B,trB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    logical, intent(in) :: trA
    logical, intent(in) :: trB

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

    C%dval=beta*C%dval

    if ((.not. trA) .and. (.not. trB)) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim2
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(i,k)*B%dval(k,j)
             end do
          end do
       end do
    else if (trA .and. (.not. trB)) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim1
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(k,i)*B%dval(k,j)
             end do
          end do
       end do
    else if ((.not. trA) .and. trB) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim2
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(i,k)*B%dval(j,k)
             end do
          end do
       end do
    else if (trA .and. trB) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim1
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(k,i)*B%dval(j,k)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_sddenref

  subroutine mm_multiply_szdenref(A,tcA,B,tcB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA
    integer, intent(in) :: tcB

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

    C%zval=beta*C%zval

    if ((tcA==0) .and. (tcB==0)) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(i,k)*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==0)) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(k,i))*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==0)) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(k,i)*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==1)) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(i,k)*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==1)) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(k,i))*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==1)) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(k,i)*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==2)) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(i,k)*B%zval(j,k)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==2)) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(k,i))*B%zval(j,k)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==2)) then
       do i=1,C%dim1
          do j=1,C%dim2
             do k=1,A%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(k,i)*B%zval(j,k)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_szdenref

  subroutine mm_multiply_sdcscsddensddenref(A,trA,B,trB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    logical, intent(in) :: trA
    logical, intent(in) :: trB

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l

    !**********************************************!

    C%dval=beta*C%dval

    if ((.not. trA) .and. (.not. trB)) then
       do k=1,A%dim2
          do l=A%iaux4(k)+1,A%iaux4(k+1)
             i=A%iaux3(l)
             do j=1,C%dim2
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(l,1)*B%dval(k,j)
             end do
          end do
       end do
    else if (trA .and. (.not. trB)) then
       do i=1,C%dim1
          do l=A%iaux4(i)+1,A%iaux4(i+1)
             k=A%iaux3(l)
             do j=1,C%dim2
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(l,1)*B%dval(k,j)
             end do
          end do
       end do
    else if ((.not. trA) .and. trB) then
       do k=1,A%dim2
          do l=A%iaux4(k)+1,A%iaux4(k+1)
             i=A%iaux3(l)
             do j=1,C%dim2
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(l,1)*B%dval(j,k)
             end do
          end do
       end do
    else if (trA .and. trB) then
       do i=1,C%dim1
          do l=A%iaux4(i)+1,A%iaux4(i+1)
             k=A%iaux3(l)
             do j=1,C%dim2
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(l,1)*B%dval(j,k)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_sdcscsddensddenref

  subroutine mm_multiply_sddensdcscsddenref(A,trA,B,trB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    logical, intent(in) :: trA
    logical, intent(in) :: trB

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l

    !**********************************************!

    C%dval=beta*C%dval

    if ((.not. trA) .and. (.not. trB)) then
       do j=1,C%dim2
          do l=B%iaux4(j)+1,B%iaux4(j+1)
             k=B%iaux3(l)
             do i=1,C%dim1
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(i,k)*B%dval(l,1)
             end do
          end do
       end do
    else if (trA .and. (.not. trB)) then
       do j=1,C%dim2
          do l=B%iaux4(j)+1,B%iaux4(j+1)
             k=B%iaux3(l)
             do i=1,C%dim1
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(k,i)*B%dval(l,1)
             end do
          end do
       end do
    else if ((.not. trA) .and. trB) then
       do k=1,A%dim2
          do l=B%iaux4(k)+1,B%iaux4(k+1)
             j=B%iaux3(l)
             do i=1,C%dim1
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(i,k)*B%dval(l,1)
             end do
          end do
       end do
    else if (trA .and. trB) then
       do k=1,A%dim1
          do l=B%iaux4(k)+1,B%iaux4(k+1)
             j=B%iaux3(l)
             do i=1,C%dim1
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(k,i)*B%dval(l,1)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_sddensdcscsddenref

  subroutine mm_multiply_sddensddensdcscref(A,trA,B,trB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    logical, intent(in) :: trA
    logical, intent(in) :: trB

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l

    !**********************************************!

    C%dval=beta*C%dval

    if ((.not. trA) .and. (.not. trB)) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim2
                C%dval(l,1)=C%dval(l,1)+alpha*A%dval(i,k)*B%dval(k,j)
             end do
          end do
       end do
    else if (trA .and. (.not. trB)) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim1
                C%dval(l,1)=C%dval(l,1)+alpha*A%dval(k,i)*B%dval(k,j)
             end do
          end do
       end do
    else if ((.not. trA) .and. trB) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim2
                C%dval(l,1)=C%dval(l,1)+alpha*A%dval(i,k)*B%dval(j,k)
             end do
          end do
       end do
    else if (trA .and. trB) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim1
                C%dval(l,1)=C%dval(l,1)+alpha*A%dval(k,i)*B%dval(j,k)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_sddensddensdcscref

  subroutine mm_multiply_sdcsrsddensddenref(A,trA,B,trB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    logical, intent(in) :: trA
    logical, intent(in) :: trB

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l

    !**********************************************!

    C%dval=beta*C%dval

    if ((.not. trA) .and. (.not. trB)) then
       do i=1,C%dim1
          do l=A%iaux3(i)+1,A%iaux3(i+1)
             k=A%iaux4(l)
             do j=1,C%dim2
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(l,1)*B%dval(k,j)
             end do
          end do
       end do
    else if (trA .and. (.not. trB)) then
       do k=1,A%dim1
          do l=A%iaux3(k)+1,A%iaux3(k+1)
             i=A%iaux4(l)
             do j=1,C%dim2
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(l,1)*B%dval(k,j)
             end do
          end do
       end do
    else if ((.not. trA) .and. trB) then
       do i=1,C%dim1
          do l=A%iaux3(i)+1,A%iaux3(i+1)
             k=A%iaux4(l)
             do j=1,C%dim2
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(l,1)*B%dval(j,k)
             end do
          end do
       end do
    else if (trA .and. trB) then
       do k=1,A%dim1
          do l=A%iaux3(k)+1,A%iaux3(k+1)
             i=A%iaux4(l)
             do j=1,C%dim2
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(l,1)*B%dval(j,k)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_sdcsrsddensddenref

  subroutine mm_multiply_sddensdcsrsddenref(A,trA,B,trB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    logical, intent(in) :: trA
    logical, intent(in) :: trB

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l

    !**********************************************!

    C%dval=beta*C%dval

    if ((.not. trA) .and. (.not. trB)) then
       do k=1,A%dim2
          do l=B%iaux3(k)+1,B%iaux3(k+1)
             j=B%iaux4(l)
             do i=1,C%dim1
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(i,k)*B%dval(l,1)
             end do
          end do
       end do
    else if (trA .and. (.not. trB)) then
       do k=1,A%dim1
          do l=B%iaux3(k)+1,B%iaux3(k+1)
             j=B%iaux4(l)
             do i=1,C%dim1
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(k,i)*B%dval(l,1)
             end do
          end do
       end do
    else if ((.not. trA) .and. trB) then
       do j=1,C%dim2
          do l=B%iaux3(j)+1,B%iaux3(j+1)
             k=B%iaux4(l)
             do i=1,C%dim1
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(i,k)*B%dval(l,1)
             end do
          end do
       end do
    else if (trA .and. trB) then
       do j=1,C%dim2
          do l=B%iaux3(j)+1,B%iaux3(j+1)
             k=B%iaux4(l)
             do i=1,C%dim1
                C%dval(i,j)=C%dval(i,j)+alpha*A%dval(k,i)*B%dval(l,1)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_sddensdcsrsddenref

  subroutine mm_multiply_sddensddensdcsrref(A,trA,B,trB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    logical, intent(in) :: trA
    logical, intent(in) :: trB

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l

    !**********************************************!

    C%dval=beta*C%dval

    if ((.not. trA) .and. (.not. trB)) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim2
                C%dval(l,1)=C%dval(l,1)+alpha*A%dval(i,k)*B%dval(k,j)
             end do
          end do
       end do
    else if (trA .and. (.not. trB)) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim1
                C%dval(l,1)=C%dval(l,1)+alpha*A%dval(k,i)*B%dval(k,j)
             end do
          end do
       end do
    else if ((.not. trA) .and. trB) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim2
                C%dval(l,1)=C%dval(l,1)+alpha*A%dval(i,k)*B%dval(j,k)
             end do
          end do
       end do
    else if (trA .and. trB) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim1
                C%dval(l,1)=C%dval(l,1)+alpha*A%dval(k,i)*B%dval(j,k)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_sddensddensdcsrref

  subroutine mm_multiply_szcscszdenszdenref(A,tcA,B,tcB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA
    integer, intent(in) :: tcB

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l

    !**********************************************!

    C%zval=beta*C%zval

    if ((tcA==0) .and. (tcB==0)) then
       do k=1,A%dim2
          do l=A%iaux4(k)+1,A%iaux4(k+1)
             i=A%iaux3(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==0)) then
       do i=1,C%dim1
          do l=A%iaux4(i)+1,A%iaux4(i+1)
             k=A%iaux3(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(l,1))*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==0)) then
       do i=1,C%dim1
          do l=A%iaux4(i)+1,A%iaux4(i+1)
             k=A%iaux3(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==1)) then
       do k=1,A%dim2
          do l=A%iaux4(k)+1,A%iaux4(k+1)
             i=A%iaux3(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==1)) then
       do i=1,C%dim1
          do l=A%iaux4(i)+1,A%iaux4(i+1)
             k=A%iaux3(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(l,1))*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==1)) then
       do i=1,C%dim1
          do l=A%iaux4(i)+1,A%iaux4(i+1)
             k=A%iaux3(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==2)) then
       do k=1,A%dim2
          do l=A%iaux4(k)+1,A%iaux4(k+1)
             i=A%iaux3(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)*B%zval(j,k)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==2)) then
       do i=1,C%dim1
          do l=A%iaux4(i)+1,A%iaux4(i+1)
             k=A%iaux3(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(l,1))*B%zval(j,k)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==2)) then
       do i=1,C%dim1
          do l=A%iaux4(i)+1,A%iaux4(i+1)
             k=A%iaux3(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)*B%zval(j,k)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_szcscszdenszdenref

  subroutine mm_multiply_szdenszcscszdenref(A,tcA,B,tcB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA
    integer, intent(in) :: tcB

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l

    !**********************************************!

    C%zval=beta*C%zval

    if ((tcA==0) .and. (tcB==0)) then
       do j=1,C%dim2
          do l=B%iaux4(j)+1,B%iaux4(j+1)
             k=B%iaux3(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(i,k)*B%zval(l,1)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==0)) then
       do j=1,C%dim2
          do l=B%iaux4(j)+1,B%iaux4(j+1)
             k=B%iaux3(l)
             do i=1,C%dim1
               C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(k,i))*B%zval(l,1)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==0)) then
       do j=1,C%dim2
          do l=B%iaux4(j)+1,B%iaux4(j+1)
             k=B%iaux3(l)
             do i=1,C%dim1
               C%zval(i,j)=C%zval(i,j)+alpha*A%zval(k,i)*B%zval(l,1)
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==1)) then
       do k=1,A%dim2
          do l=B%iaux4(k)+1,B%iaux4(k+1)
             j=B%iaux3(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(i,k)*conjg(B%zval(l,1))
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==1)) then
       do k=1,A%dim1
          do l=B%iaux4(k)+1,B%iaux4(k+1)
             j=B%iaux3(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(k,i))*conjg(B%zval(l,1))
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==1)) then
       do k=1,A%dim1
          do l=B%iaux4(k)+1,B%iaux4(k+1)
             j=B%iaux3(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(k,i)*conjg(B%zval(l,1))
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==2)) then
       do k=1,A%dim2
          do l=B%iaux4(k)+1,B%iaux4(k+1)
             j=B%iaux3(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(i,k)*B%zval(l,1)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==2)) then
       do k=1,A%dim1
          do l=B%iaux4(k)+1,B%iaux4(k+1)
             j=B%iaux3(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(k,i))*B%zval(l,1)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==2)) then
       do k=1,A%dim1
          do l=B%iaux4(k)+1,B%iaux4(k+1)
             j=B%iaux3(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(k,i)*B%zval(l,1)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_szdenszcscszdenref

  subroutine mm_multiply_szdenszdenszcscref(A,tcA,B,tcB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA
    integer, intent(in) :: tcB

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l

    !**********************************************!

    C%zval=beta*C%zval

    if ((tcA==0) .and. (tcB==0)) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim2
                C%zval(l,1)=C%zval(l,1)+alpha*A%zval(i,k)*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==0)) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim1
                C%zval(l,1)=C%zval(l,1)+alpha*conjg(A%zval(k,i))*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==0)) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim1
                C%zval(l,1)=C%zval(l,1)+alpha*A%zval(k,i)*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==1)) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim2
                C%zval(l,1)=C%zval(l,1)+alpha*A%zval(i,k)*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==1)) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim1
                C%zval(l,1)=C%zval(l,1)+alpha*conjg(A%zval(k,i))*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==1)) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim1
                C%zval(l,1)=C%zval(l,1)+alpha*A%zval(k,i)*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==2)) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim2
                C%zval(l,1)=C%zval(l,1)+alpha*A%zval(i,k)*B%zval(j,k)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==2)) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim1
                C%zval(l,1)=C%zval(l,1)+alpha*conjg(A%zval(k,i))*B%zval(j,k)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==2)) then
       do j=1,C%dim2
          do l=C%iaux4(j)+1,C%iaux4(j+1)
             i=C%iaux3(l)
             do k=1,A%dim1
                C%zval(l,1)=C%zval(l,1)+alpha*A%zval(k,i)*B%zval(j,k)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_szdenszdenszcscref

  subroutine mm_multiply_szcsrszdenszdenref(A,tcA,B,tcB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA
    integer, intent(in) :: tcB

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l

    !**********************************************!

    C%zval=beta*C%zval

    if ((tcA==0) .and. (tcB==0)) then
       do i=1,C%dim1
          do l=A%iaux3(i)+1,A%iaux3(i+1)
             k=A%iaux4(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==0)) then
       do k=1,A%dim1
          do l=A%iaux3(k)+1,A%iaux3(k+1)
             i=A%iaux4(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(l,1))*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==0)) then
       do k=1,A%dim1
          do l=A%iaux3(k)+1,A%iaux3(k+1)
             i=A%iaux4(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==1)) then
       do i=1,C%dim1
          do l=A%iaux3(i)+1,A%iaux3(i+1)
             k=A%iaux4(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==1)) then
       do k=1,A%dim1
          do l=A%iaux3(k)+1,A%iaux3(k+1)
             i=A%iaux4(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(l,1))*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==1)) then
       do k=1,A%dim1
          do l=A%iaux3(k)+1,A%iaux3(k+1)
             i=A%iaux4(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==2)) then
       do i=1,C%dim1
          do l=A%iaux3(i)+1,A%iaux3(i+1)
             k=A%iaux4(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)*B%zval(j,k)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==2)) then
       do k=1,A%dim1
          do l=A%iaux3(k)+1,A%iaux3(k+1)
             i=A%iaux4(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(l,1))*B%zval(j,k)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==2)) then
       do k=1,A%dim1
          do l=A%iaux3(k)+1,A%iaux3(k+1)
             i=A%iaux4(l)
             do j=1,C%dim2
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(l,1)*B%zval(j,k)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_szcsrszdenszdenref

  subroutine mm_multiply_szdenszcsrszdenref(A,tcA,B,tcB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA
    integer, intent(in) :: tcB

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l

    !**********************************************!

    C%zval=beta*C%zval

    if ((tcA==0) .and. (tcB==0)) then
       do k=1,A%dim2
          do l=B%iaux3(k)+1,B%iaux3(k+1)
             j=B%iaux4(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(i,k)*B%zval(l,1)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==0)) then
       do k=1,A%dim1
          do l=B%iaux3(k)+1,B%iaux3(k+1)
             j=B%iaux4(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(k,i))*B%zval(l,1)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==0)) then
       do k=1,A%dim1
          do l=B%iaux3(k)+1,B%iaux3(k+1)
             j=B%iaux4(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(k,i)*B%zval(l,1)
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==1)) then
       do j=1,C%dim2
          do l=B%iaux3(j)+1,B%iaux3(j+1)
             k=B%iaux4(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(i,k)*conjg(B%zval(l,1))
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==1)) then
       do j=1,C%dim2
          do l=B%iaux3(j)+1,B%iaux3(j+1)
             k=B%iaux4(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(k,i))*conjg(B%zval(l,1))
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==1)) then
       do j=1,C%dim2
          do l=B%iaux3(j)+1,B%iaux3(j+1)
             k=B%iaux4(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(k,i)*conjg(B%zval(l,1))
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==2)) then
       do j=1,C%dim2
          do l=B%iaux3(j)+1,B%iaux3(j+1)
             k=B%iaux4(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(i,k)*B%zval(l,1)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==2)) then
       do j=1,C%dim2
          do l=B%iaux3(j)+1,B%iaux3(j+1)
             k=B%iaux4(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*conjg(A%zval(k,i))*B%zval(l,1)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==2)) then
       do j=1,C%dim2
          do l=B%iaux3(j)+1,B%iaux3(j+1)
             k=B%iaux4(l)
             do i=1,C%dim1
                C%zval(i,j)=C%zval(i,j)+alpha*A%zval(k,i)*B%zval(l,1)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_szdenszcsrszdenref

  subroutine mm_multiply_szdenszdenszcsrref(A,tcA,B,tcB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA
    integer, intent(in) :: tcB

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l

    !**********************************************!

    C%zval=beta*C%zval

    if ((tcA==0) .and. (tcB==0)) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim2
                C%zval(l,1)=C%zval(l,1)+alpha*A%zval(i,k)*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==0)) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim1
                C%zval(l,1)=C%zval(l,1)+alpha*conjg(A%zval(k,i))*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==0)) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim1
                C%zval(l,1)=C%zval(l,1)+alpha*A%zval(k,i)*B%zval(k,j)
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==1)) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim2
                C%zval(l,1)=C%zval(l,1)+alpha*A%zval(i,k)*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==1)) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim1
                C%zval(l,1)=C%zval(l,1)+alpha*conjg(A%zval(k,i))*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==1)) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim1
                C%zval(l,1)=C%zval(l,1)+alpha*A%zval(k,i)*conjg(B%zval(j,k))
             end do
          end do
       end do
    else if ((tcA==0) .and. (tcB==2)) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim2
                C%zval(l,1)=C%zval(l,1)+alpha*A%zval(i,k)*B%zval(j,k)
             end do
          end do
       end do
    else if ((tcA==1) .and. (tcB==2)) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim1
                C%zval(l,1)=C%zval(l,1)+alpha*conjg(A%zval(k,i))*B%zval(j,k)
             end do
          end do
       end do
    else if ((tcA==2) .and. (tcB==2)) then
       do i=1,C%dim1
          do l=C%iaux3(i)+1,C%iaux3(i+1)
             j=C%iaux4(l)
             do k=1,A%dim1
                C%zval(l,1)=C%zval(l,1)+alpha*A%zval(k,i)*B%zval(j,k)
             end do
          end do
       end do
    end if

  end subroutine mm_multiply_szdenszdenszcsrref

  !================================================!
  ! implementation: sparse-dense 1D distributed    !
  !================================================!

#ifdef PSP
  subroutine mm_multiply_pddbcpdcscpddbct1D(A,B,trB,C,alpha,beta)
    implicit none
    include 'mpif.h'

    !**** INPUT ***********************************!

    logical, intent(in) :: trB

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l, m
    integer :: n_comm, nnz_recv, loc_dim_recv, info
    integer, allocatable :: col_ptr_recv(:), row_ind_recv(:)

    real(dp), allocatable :: dval_recv(:)

    !**** EXTERNAL ********************************!

    integer, external :: indxl2g

    !**********************************************!

    C%dval=beta*C%dval

    do n_comm=0,ms_mpi_size-1
       if (n_comm==ms_mpi_rank) then
          loc_dim_recv=B%spm%loc_dim2
          nnz_recv=B%spm%nnz
       end if
       call mpi_bcast(loc_dim_recv,1,mpi_integer,n_comm,ms_mpi_comm,info)
       call mpi_bcast(nnz_recv,    1,mpi_integer,n_comm,ms_mpi_comm,info)
       allocate(col_ptr_recv(loc_dim_recv+1))
       allocate(row_ind_recv(nnz_recv))
       allocate(dval_recv(nnz_recv))
       if (n_comm==ms_mpi_rank) then
          col_ptr_recv(1:loc_dim_recv+1)=B%spm%col_ptr(1:loc_dim_recv+1)
          row_ind_recv(1:nnz_recv)=B%spm%row_ind(1:nnz_recv)
          dval_recv(1:nnz_recv)=B%spm%dval(1:nnz_recv)
       end if
       call mpi_bcast(col_ptr_recv(1),loc_dim_recv+1,mpi_integer,         n_comm,ms_mpi_comm,info)
       call mpi_bcast(row_ind_recv(1),nnz_recv,      mpi_integer,         n_comm,ms_mpi_comm,info)
       call mpi_bcast(dval_recv(1),   nnz_recv,      mpi_double_precision,n_comm,ms_mpi_comm,info)
       if (.not. trB) then
          do j=1,loc_dim_recv
             do m=col_ptr_recv(j),col_ptr_recv(j+1)-1
                k=indxl2g(row_ind_recv(m),B%spm%desc(6),n_comm,B%spm%desc(8),ms_mpi_size)
                C%dval(:,j)=C%dval(:,j)+alpha*A%dval(:,k)*dval_recv(m)
             end do
          end do
       else
          do k=1,loc_dim_recv
             do m=col_ptr_recv(k),col_ptr_recv(k+1)-1
                j=indxl2g(row_ind_recv(m),B%spm%desc(6),n_comm,B%spm%desc(8),ms_mpi_size)
                C%dval(:,j)=C%dval(:,j)+alpha*A%dval(:,k)*dval_recv(m)
             end do
          end do
       end if
       deallocate(dval_recv)
       deallocate(row_ind_recv)
       deallocate(col_ptr_recv)
    end do

  end subroutine mm_multiply_pddbcpdcscpddbct1D

  subroutine mm_multiply_pdcscpddbcpddbct1D(A,trA,B,C,alpha,beta)
    implicit none
    include 'mpif.h'

    !**** INPUT ***********************************!

    logical, intent(in) :: trA

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l, m
    integer :: n_comm, nnz_recv, loc_dim_recv, info
    integer, allocatable :: col_ptr_recv(:), row_ind_recv(:)

    real(dp), allocatable :: dval_recv(:)

    !**** EXTERNAL ********************************!

    integer, external :: indxl2g

    !**********************************************!

    C%dval=beta*C%dval

    do n_comm=0,ms_mpi_size-1
       if (n_comm==ms_mpi_rank) then
          loc_dim_recv=A%spm%loc_dim2
          nnz_recv=A%spm%nnz
       end if
       call mpi_bcast(loc_dim_recv,1,mpi_integer,n_comm,ms_mpi_comm,info)
       call mpi_bcast(nnz_recv,    1,mpi_integer,n_comm,ms_mpi_comm,info)
       allocate(col_ptr_recv(loc_dim_recv+1))
       allocate(row_ind_recv(nnz_recv))
       allocate(dval_recv(nnz_recv))
       if (n_comm==ms_mpi_rank) then
          col_ptr_recv(1:loc_dim_recv+1)=A%spm%col_ptr(1:loc_dim_recv+1)
          row_ind_recv(1:nnz_recv)=A%spm%row_ind(1:nnz_recv)
          dval_recv(1:nnz_recv)=A%spm%dval(1:nnz_recv)
       end if
       call mpi_bcast(col_ptr_recv(1),loc_dim_recv+1,mpi_integer,         n_comm,ms_mpi_comm,info)
       call mpi_bcast(row_ind_recv(1),nnz_recv,      mpi_integer,         n_comm,ms_mpi_comm,info)
       call mpi_bcast(dval_recv(1),   nnz_recv,      mpi_double_precision,n_comm,ms_mpi_comm,info)
       if (.not. trA) then
          do l=1,loc_dim_recv
             k=indxl2g(l,A%spm%desc(5),n_comm,A%spm%desc(7),ms_mpi_size)
             do m=col_ptr_recv(l),col_ptr_recv(l+1)-1
                i=row_ind_recv(m)
                C%dval(i,:)=C%dval(i,:)+alpha*dval_recv(m)*B%dval(k,:)
             end do
          end do
       else
          do l=1,loc_dim_recv
             i=indxl2g(l,A%spm%desc(5),n_comm,A%spm%desc(7),ms_mpi_size)
             do m=col_ptr_recv(l),col_ptr_recv(l+1)-1
                k=row_ind_recv(m)
                C%dval(i,:)=C%dval(i,:)+alpha*dval_recv(m)*B%dval(k,:)
             end do
          end do
       end if
       deallocate(dval_recv)
       deallocate(row_ind_recv)
       deallocate(col_ptr_recv)
    end do

  end subroutine mm_multiply_pdcscpddbcpddbct1D

  subroutine mm_multiply_pddbcpddbcpdcsct1D(A,trA,B,trB,C,alpha,beta)
    implicit none
    include 'mpif.h'

    !**** INPUT ***********************************!

    logical, intent(in) :: trA
    logical, intent(in) :: trB

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l, m
    integer :: n_comm, nnz_recv, loc_dim_recv, info
    integer, allocatable :: col_ptr_recv(:), row_ind_recv(:)

    real(dp), allocatable :: dval_recv(:)

    !**** EXTERNAL ********************************!

    integer, external :: indxl2g

    !**********************************************!

    do n_comm=0,ms_mpi_size-1
       if (n_comm==ms_mpi_rank) then
          loc_dim_recv=C%spm%loc_dim2
          nnz_recv=C%spm%nnz
       end if
       call mpi_bcast(loc_dim_recv,1,mpi_integer,n_comm,ms_mpi_comm,info)
       call mpi_bcast(nnz_recv,    1,mpi_integer,n_comm,ms_mpi_comm,info)
       allocate(col_ptr_recv(loc_dim_recv+1))
       allocate(row_ind_recv(nnz_recv))
       allocate(dval_recv(nnz_recv))
       if (n_comm==ms_mpi_rank) then
          col_ptr_recv(1:loc_dim_recv+1)=C%spm%col_ptr(1:loc_dim_recv+1)
          row_ind_recv(1:nnz_recv)=C%spm%row_ind(1:nnz_recv)
       end if
       call mpi_bcast(col_ptr_recv(1),loc_dim_recv+1,mpi_integer,n_comm,ms_mpi_comm,info)
       call mpi_bcast(row_ind_recv(1),nnz_recv,      mpi_integer,n_comm,ms_mpi_comm,info)
       dval_recv=0.0_dp
       if (.not. trA) then
          do l=1,loc_dim_recv
             j=indxl2g(l,C%spm%desc(5),n_comm,C%spm%desc(7),ms_mpi_size)
             do m=col_ptr_recv(l),col_ptr_recv(l+1)-1
                i=row_ind_recv(m)
                do k=1,A%iaux2(2)
                   dval_recv(m)=dval_recv(m)+A%dval(i,k)*B%dval(j,k)
                end do
             end do
          end do
       else
          do j=1,loc_dim_recv
             do m=col_ptr_recv(j),col_ptr_recv(j+1)-1
                i=indxl2g(row_ind_recv(m),C%spm%desc(6),n_comm,C%spm%desc(8),ms_mpi_size)
                do k=1,A%iaux2(1)
                   dval_recv(m)=dval_recv(m)+A%dval(k,i)*B%dval(k,j)
                end do
             end do
          end do
       end if
       if (n_comm==ms_mpi_rank) then
          call mpi_reduce(mpi_in_place,dval_recv(1),nnz_recv,mpi_double_precision,mpi_sum,n_comm,ms_mpi_comm,info)
          C%spm%dval=beta*C%spm%dval+alpha*dval_recv
       else
          call mpi_reduce(dval_recv(1),mpi_proc_null,nnz_recv,mpi_double_precision,mpi_sum,n_comm,ms_mpi_comm,info)
       end if
       deallocate(dval_recv)
       deallocate(row_ind_recv)
       deallocate(col_ptr_recv)
    end do

  end subroutine mm_multiply_pddbcpddbcpdcsct1D

  subroutine mm_multiply_pzdbcpzcscpzdbct1D(A,B,tcB,C,alpha,beta)
    implicit none
    include 'mpif.h'

    !**** INPUT ***********************************!

    integer, intent(in) :: tcB

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l, m
    integer :: n_comm, nnz_recv, loc_dim_recv, info
    integer, allocatable :: col_ptr_recv(:), row_ind_recv(:)

    complex(dp), allocatable :: zval_recv(:)

    !**** EXTERNAL ********************************!

    integer, external :: indxl2g

    !**********************************************!

    C%zval=beta*C%zval

    do n_comm=0,ms_mpi_size-1
       if (n_comm==ms_mpi_rank) then
          loc_dim_recv=B%spm%loc_dim2
          nnz_recv=B%spm%nnz
       end if
       call mpi_bcast(loc_dim_recv,1,mpi_integer,n_comm,ms_mpi_comm,info)
       call mpi_bcast(nnz_recv,    1,mpi_integer,n_comm,ms_mpi_comm,info)
       allocate(col_ptr_recv(loc_dim_recv+1))
       allocate(row_ind_recv(nnz_recv))
       allocate(zval_recv(nnz_recv))
       if (n_comm==ms_mpi_rank) then
          col_ptr_recv(1:loc_dim_recv+1)=B%spm%col_ptr(1:loc_dim_recv+1)
          row_ind_recv(1:nnz_recv)=B%spm%row_ind(1:nnz_recv)
          zval_recv(1:nnz_recv)=B%spm%zval(1:nnz_recv)
       end if
       call mpi_bcast(col_ptr_recv(1),loc_dim_recv+1,mpi_integer,       n_comm,ms_mpi_comm,info)
       call mpi_bcast(row_ind_recv(1),nnz_recv,      mpi_integer,       n_comm,ms_mpi_comm,info)
       call mpi_bcast(zval_recv(1),   nnz_recv,      mpi_double_complex,n_comm,ms_mpi_comm,info)
       if (tcB==0) then
          do j=1,loc_dim_recv
             do m=col_ptr_recv(j),col_ptr_recv(j+1)-1
                k=indxl2g(row_ind_recv(m),B%spm%desc(6),n_comm,B%spm%desc(8),ms_mpi_size)
                C%zval(:,j)=C%zval(:,j)+alpha*A%zval(:,k)*zval_recv(m)
             end do
          end do
       else if (tcB==1) then
          do k=1,loc_dim_recv
             do m=col_ptr_recv(k),col_ptr_recv(k+1)-1
                j=indxl2g(row_ind_recv(m),B%spm%desc(6),n_comm,B%spm%desc(8),ms_mpi_size)
                C%zval(:,j)=C%zval(:,j)+alpha*A%zval(:,k)*conjg(zval_recv(m))
             end do
          end do
       else if (tcB==2) then
          do k=1,loc_dim_recv
             do m=col_ptr_recv(k),col_ptr_recv(k+1)-1
                j=indxl2g(row_ind_recv(m),B%spm%desc(6),n_comm,B%spm%desc(8),ms_mpi_size)
                C%zval(:,j)=C%zval(:,j)+alpha*A%zval(:,k)*zval_recv(m)
             end do
          end do
       end if
       deallocate(zval_recv)
       deallocate(row_ind_recv)
       deallocate(col_ptr_recv)
    end do

  end subroutine mm_multiply_pzdbcpzcscpzdbct1D

  subroutine mm_multiply_pzcscpzdbcpzdbct1D(A,tcA,B,C,alpha,beta)
    implicit none
    include 'mpif.h'

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l, m
    integer :: n_comm, nnz_recv, loc_dim_recv, info
    integer, allocatable :: col_ptr_recv(:), row_ind_recv(:)

    complex(dp), allocatable :: zval_recv(:)

    !**** EXTERNAL ********************************!

    integer, external :: indxl2g

    !**********************************************!

    C%zval=beta*C%zval

    do n_comm=0,ms_mpi_size-1
       if (n_comm==ms_mpi_rank) then
          loc_dim_recv=A%spm%loc_dim2
          nnz_recv=A%spm%nnz
       end if
       call mpi_bcast(loc_dim_recv,1,mpi_integer,n_comm,ms_mpi_comm,info)
       call mpi_bcast(nnz_recv,    1,mpi_integer,n_comm,ms_mpi_comm,info)
       allocate(col_ptr_recv(loc_dim_recv+1))
       allocate(row_ind_recv(nnz_recv))
       allocate(zval_recv(nnz_recv))
       if (n_comm==ms_mpi_rank) then
          col_ptr_recv(1:loc_dim_recv+1)=A%spm%col_ptr(1:loc_dim_recv+1)
          row_ind_recv(1:nnz_recv)=A%spm%row_ind(1:nnz_recv)
          zval_recv(1:nnz_recv)=A%spm%zval(1:nnz_recv)
       end if
       call mpi_bcast(col_ptr_recv(1),loc_dim_recv+1,mpi_integer,       n_comm,ms_mpi_comm,info)
       call mpi_bcast(row_ind_recv(1),nnz_recv,      mpi_integer,       n_comm,ms_mpi_comm,info)
       call mpi_bcast(zval_recv(1),   nnz_recv,      mpi_double_complex,n_comm,ms_mpi_comm,info)
       if (tcA==0) then
          do l=1,loc_dim_recv
             k=indxl2g(l,A%spm%desc(5),n_comm,A%spm%desc(7),ms_mpi_size)
             do m=col_ptr_recv(l),col_ptr_recv(l+1)-1
                i=row_ind_recv(m)
                C%zval(i,:)=C%zval(i,:)+alpha*zval_recv(m)*B%zval(k,:)
             end do
          end do
       else if (tcA==1) then
          do l=1,loc_dim_recv
             i=indxl2g(l,A%spm%desc(5),n_comm,A%spm%desc(7),ms_mpi_size)
             do m=col_ptr_recv(l),col_ptr_recv(l+1)-1
                k=row_ind_recv(m)
                C%zval(i,:)=C%zval(i,:)+alpha*conjg(zval_recv(m))*B%zval(k,:)
             end do
          end do
       else if (tcA==2) then
          do l=1,loc_dim_recv
             i=indxl2g(l,A%spm%desc(5),n_comm,A%spm%desc(7),ms_mpi_size)
             do m=col_ptr_recv(l),col_ptr_recv(l+1)-1
                k=row_ind_recv(m)
                C%zval(i,:)=C%zval(i,:)+alpha*zval_recv(m)*B%zval(k,:)
             end do
          end do
       end if
       deallocate(zval_recv)
       deallocate(row_ind_recv)
       deallocate(col_ptr_recv)
    end do

  end subroutine mm_multiply_pzcscpzdbcpzdbct1D

  subroutine mm_multiply_pzdbcpzdbcpzcsct1D(A,tcA,B,tcB,C,alpha,beta)
    implicit none
    include 'mpif.h'

    !**** INPUT ***********************************!

    integer, intent(in) :: tcA
    integer, intent(in) :: tcB

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    type(matrix), intent(in) :: A
    type(matrix), intent(in) :: B

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: i, j, k, l, m
    integer :: n_comm, nnz_recv, loc_dim_recv, info
    integer, allocatable :: col_ptr_recv(:), row_ind_recv(:)

    complex(dp), allocatable :: zval_recv(:)

    !**** EXTERNAL ********************************!

    integer, external :: indxl2g

    !**********************************************!

    do n_comm=0,ms_mpi_size-1
       if (n_comm==ms_mpi_rank) then
          loc_dim_recv=C%spm%loc_dim2
          nnz_recv=C%spm%nnz
       end if
       call mpi_bcast(loc_dim_recv,1,mpi_integer,n_comm,ms_mpi_comm,info)
       call mpi_bcast(nnz_recv,    1,mpi_integer,n_comm,ms_mpi_comm,info)
       allocate(col_ptr_recv(loc_dim_recv+1))
       allocate(row_ind_recv(nnz_recv))
       allocate(zval_recv(nnz_recv))
       if (n_comm==ms_mpi_rank) then
          col_ptr_recv(1:loc_dim_recv+1)=C%spm%col_ptr(1:loc_dim_recv+1)
          row_ind_recv(1:nnz_recv)=C%spm%row_ind(1:nnz_recv)
       end if
       call mpi_bcast(col_ptr_recv(1),loc_dim_recv+1,mpi_integer,n_comm,ms_mpi_comm,info)
       call mpi_bcast(row_ind_recv(1),nnz_recv,      mpi_integer,n_comm,ms_mpi_comm,info)
       zval_recv=cmplx_0
       if (tcA==0) then
          if (tcB==1) then
             do l=1,loc_dim_recv
                j=indxl2g(l,C%spm%desc(5),n_comm,C%spm%desc(7),ms_mpi_size)
                do m=col_ptr_recv(l),col_ptr_recv(l+1)-1
                   i=row_ind_recv(m)
                   do k=1,A%iaux2(2)
                      zval_recv(m)=zval_recv(m)+A%zval(i,k)*conjg(B%zval(j,k))
                   end do
                end do
             end do
          else if (tcB==2) then
             do l=1,loc_dim_recv
                j=indxl2g(l,C%spm%desc(5),n_comm,C%spm%desc(7),ms_mpi_size)
                do m=col_ptr_recv(l),col_ptr_recv(l+1)-1
                   i=row_ind_recv(m)
                   do k=1,A%iaux2(2)
                      zval_recv(m)=zval_recv(m)+A%zval(i,k)*B%zval(j,k)
                   end do
                end do
             end do
          end if
       else if (tcB==0) then
          if (tcA==1) then
             do j=1,loc_dim_recv
                do m=col_ptr_recv(j),col_ptr_recv(j+1)-1
                   i=indxl2g(row_ind_recv(m),C%spm%desc(6),n_comm,C%spm%desc(8),ms_mpi_size)
                   do k=1,A%iaux2(1)
                      zval_recv(m)=zval_recv(m)+conjg(A%zval(k,i))*B%zval(k,j)
                   end do
                end do
             end do
          else if (tcA==2) then
             do j=1,loc_dim_recv
                do m=col_ptr_recv(j),col_ptr_recv(j+1)-1
                   i=indxl2g(row_ind_recv(m),C%spm%desc(6),n_comm,C%spm%desc(8),ms_mpi_size)
                   do k=1,A%iaux2(1)
                      zval_recv(m)=zval_recv(m)+A%zval(k,i)*B%zval(k,j)
                   end do
                end do
             end do
          end if
       end if
       if (n_comm==ms_mpi_rank) then
          call mpi_reduce(mpi_in_place,zval_recv(1),nnz_recv,mpi_double_complex,mpi_sum,n_comm,ms_mpi_comm,info)
          C%spm%zval=beta*C%spm%zval+alpha*zval_recv
       else
          call mpi_reduce(zval_recv(1),mpi_proc_null,nnz_recv,mpi_double_complex,mpi_sum,n_comm,ms_mpi_comm,info)
       end if
       deallocate(zval_recv)
       deallocate(row_ind_recv)
       deallocate(col_ptr_recv)
    end do

  end subroutine mm_multiply_pzdbcpzdbcpzcsct1D
#endif

end module MatrixSwitch_mm_multiply
