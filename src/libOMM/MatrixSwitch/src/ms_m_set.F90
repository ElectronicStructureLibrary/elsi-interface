module ms_m_set
  use MatrixSwitch_ops

  implicit none

contains

  !================================================!
  ! implementation: reference                      !
  !================================================!

  subroutine m_set_sddenref(C,seC,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: seC

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: luC

    integer :: i, j

    !**********************************************!

    call process_seM(seC,luC)

    if (luC==0) then
       do i=1,C%dim1
          do j=1,C%dim2
             if (i/=j) then
                C%dval(i,j)=alpha
             else
                C%dval(i,j)=beta
             end if
          end do
       end do
    else if (luC==1) then
       do i=1,C%dim1
          do j=i,C%dim2
             if (i/=j) then
                C%dval(i,j)=alpha
             else
                C%dval(i,j)=beta
             end if
          end do
       end do
    else if (luC==2) then
       do i=1,C%dim1
          do j=1,i
             if (i/=j) then
                C%dval(i,j)=alpha
             else
                C%dval(i,j)=beta
             end if
          end do
       end do
    end if

  end subroutine m_set_sddenref

  subroutine m_set_szdenref(C,seC,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: seC

    complex(dp), intent(in) :: alpha
    complex(dp), intent(in) :: beta

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C

    !**** INTERNAL ********************************!

    integer :: luC

    integer :: i, j

    !**********************************************!

    call process_seM(seC,luC)

    if (luC==0) then
       do i=1,C%dim1
          do j=1,C%dim2
             if (i/=j) then
                C%zval(i,j)=alpha
             else
                C%zval(i,j)=beta
             end if
          end do
       end do
    else if (luC==1) then
       do i=1,C%dim1
          do j=i,C%dim2
             if (i/=j) then
                C%zval(i,j)=alpha
             else
                C%zval(i,j)=beta
             end if
          end do
       end do
    else if (luC==2) then
       do i=1,C%dim1
          do j=1,i
             if (i/=j) then
                C%zval(i,j)=alpha
             else
                C%zval(i,j)=beta
             end if
          end do
       end do
    end if

  end subroutine m_set_szdenref

end module ms_m_set
