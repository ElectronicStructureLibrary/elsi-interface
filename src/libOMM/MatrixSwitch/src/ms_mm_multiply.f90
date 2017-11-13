module ms_mm_multiply
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

end module ms_mm_multiply
