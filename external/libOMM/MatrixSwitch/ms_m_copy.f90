module ms_m_copy
  use MatrixSwitch_ops

  implicit none

contains

  !================================================!
  ! implementation: reference                      !
  !================================================!

  subroutine m_copy_sdensdenref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if

  end subroutine m_copy_sdensdenref

  subroutine m_copy_sdensdenref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          do i=1,m_name%dim1
             do j=1,m_name%dim2
                if (abs(A%dval(i,j))>abs_threshold) then
                   m_name%dval(i,j)=A%dval(i,j)-soft_threshold*A%dval(i,j)/abs(A%dval(i,j))
                 else
                   m_name%dval(i,j)=0.0_dp
                end if
             end do
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          do i=1,m_name%dim1
             do j=1,m_name%dim2
                if (abs(A%zval(i,j))>abs_threshold) then
                   m_name%zval(i,j)=A%zval(i,j)-soft_threshold*A%zval(i,j)/abs(A%zval(i,j))
                else
                   m_name%zval(i,j)=cmplx_0
                end if
             end do
          end do
       end if

  end subroutine m_copy_sdensdenref_thre

  subroutine m_copy_pdbcpdbcref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

       allocate(m_name%iaux1(9))
       m_name%iaux1_is_allocated=.true.
       allocate(m_name%iaux2(2))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux1=A%iaux1
       m_name%iaux2=A%iaux2
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),m_name%iaux2(2)))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),m_name%iaux2(2)))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if

  end subroutine m_copy_pdbcpdbcref

  subroutine m_copy_pdbcpdbcref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j

    !**********************************************!

       allocate(m_name%iaux1(9))
       m_name%iaux1_is_allocated=.true.
       allocate(m_name%iaux2(2))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux1=A%iaux1
       m_name%iaux2=A%iaux2
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),m_name%iaux2(2)))
          m_name%dval_is_allocated=.true.
          do i=1,m_name%iaux2(1)
             do j=1,m_name%iaux2(2)
                if (abs(A%dval(i,j))>abs_threshold) then
                   m_name%dval(i,j)=A%dval(i,j)-soft_threshold*A%dval(i,j)/abs(A%dval(i,j))
                else
                   m_name%dval(i,j)=0.0_dp
                end if
             end do
          end do
       else
          allocate(m_name%zval(m_name%iaux2(1),m_name%iaux2(2)))
          m_name%zval_is_allocated=.true.
          do i=1,m_name%iaux2(1)
             do j=1,m_name%iaux2(2)
                if (abs(A%zval(i,j))>abs_threshold) then
                   m_name%zval(i,j)=A%zval(i,j)-soft_threshold*A%zval(i,j)/abs(A%zval(i,j))
                else
                   m_name%zval(i,j)=cmplx_0
                end if
             end do
          end do
       end if

  end subroutine m_copy_pdbcpdbcref_thre

  subroutine m_copy_scooscooref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       m_name%iaux3=A%iaux3
       m_name%iaux4=A%iaux4
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if

  end subroutine m_copy_scooscooref

  subroutine m_copy_sdenscooref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=m_name%dim1*m_name%dim2
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          k=0
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                k=k+1
                m_name%iaux3(k)=j
                m_name%iaux4(k)=i
                m_name%dval(k,1)=A%dval(j,i)
             end do
          end do
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          k=0
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                k=k+1
                m_name%iaux3(k)=j
                m_name%iaux4(k)=i
                m_name%zval(k,1)=A%zval(j,i)
             end do
          end do
       end if

  end subroutine m_copy_sdenscooref

  subroutine m_copy_sdenscooref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=0
       if (m_name%is_real) then
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%dval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%iaux2(1)))
             m_name%iaux3_is_allocated=.true.
             allocate(m_name%iaux4(m_name%iaux2(1)))
             m_name%iaux4_is_allocated=.true.
             allocate(m_name%dval(m_name%iaux2(1),1))
             m_name%dval_is_allocated=.true.
             k=0
             do i=1,m_name%dim2
                do j=1,m_name%dim1
                   if (abs(A%dval(j,i))>abs_threshold) then
                      k=k+1
                      m_name%iaux3(k)=j
                      m_name%iaux4(k)=i
                      m_name%dval(k,1)=A%dval(j,i)-soft_threshold*A%dval(j,i)/abs(A%dval(j,i))
                   end if
                end do
             end do
          end if
       else
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%zval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%iaux2(1)))
             m_name%iaux3_is_allocated=.true.
             allocate(m_name%iaux4(m_name%iaux2(1)))
             m_name%iaux4_is_allocated=.true.
             allocate(m_name%zval(m_name%iaux2(1),1))
             m_name%zval_is_allocated=.true.
             k=0
             do i=1,m_name%dim2
                do j=1,m_name%dim1
                   if (abs(A%zval(j,i))>abs_threshold) then
                      k=k+1
                      m_name%iaux3(k)=j
                      m_name%iaux4(k)=i
                      m_name%zval(k,1)=A%zval(j,i)-soft_threshold*A%zval(j,i)/abs(A%zval(j,i))
                   end if
                end do
             end do
          end if
       end if

  end subroutine m_copy_sdenscooref_thre

  subroutine m_copy_scscscscref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%dim2+1))
       m_name%iaux4_is_allocated=.true.
       m_name%iaux3=A%iaux3
       m_name%iaux4=A%iaux4
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if

  end subroutine m_copy_scscscscref

  subroutine m_copy_scsrscsrref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%dim1+1))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       m_name%iaux3=A%iaux3
       m_name%iaux4=A%iaux4
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if

  end subroutine m_copy_scsrscsrref

  subroutine m_copy_sdenscscref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=m_name%dim1*m_name%dim2
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%dim2+1))
       m_name%iaux4_is_allocated=.true.
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          k=0
          m_name%iaux4(1)=0
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                k=k+1
                m_name%iaux3(k)=j
                m_name%dval(k,1)=A%dval(j,i)
             end do
             m_name%iaux4(i+1)=k
          end do
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          k=0
          m_name%iaux4(1)=0
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                k=k+1
                m_name%iaux3(k)=j
                m_name%zval(k,1)=A%zval(j,i)
             end do
             m_name%iaux4(i)=k+1
          end do
       end if

  end subroutine m_copy_sdenscscref

  subroutine m_copy_sdenscscref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=0
       if (m_name%is_real) then
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%dval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%iaux2(1)))
             m_name%iaux3_is_allocated=.true.
             allocate(m_name%iaux4(m_name%dim2+1))
             m_name%iaux4_is_allocated=.true.
             allocate(m_name%dval(m_name%iaux2(1),1))
             m_name%dval_is_allocated=.true.
             k=0
             m_name%iaux4(1)=0
             do i=1,m_name%dim2
                do j=1,m_name%dim1
                   if (abs(A%dval(j,i))>abs_threshold) then
                      k=k+1
                      m_name%iaux3(k)=j
                      m_name%dval(k,1)=A%dval(j,i)-soft_threshold*A%dval(j,i)/abs(A%dval(j,i))
                   end if
                end do
                m_name%iaux4(i+1)=k
             end do
          end if
       else
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%zval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%iaux2(1)))
             m_name%iaux3_is_allocated=.true.
             allocate(m_name%iaux4(m_name%dim2+1))
             m_name%iaux4_is_allocated=.true.
             allocate(m_name%zval(m_name%iaux2(1),1))
             m_name%zval_is_allocated=.true.
             k=0
             m_name%iaux4(1)=0
             do i=1,m_name%dim2
                do j=1,m_name%dim1
                   if (abs(A%zval(j,i))>abs_threshold) then
                      k=k+1
                      m_name%iaux3(k)=j
                      m_name%zval(k,1)=A%zval(j,i)-soft_threshold*A%zval(j,i)/abs(A%zval(j,i))
                   end if
                end do
                m_name%iaux4(i)=k+1
             end do
          end if
       end if

  end subroutine m_copy_sdenscscref_thre

  subroutine m_copy_sdenscsrref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=m_name%dim1*m_name%dim2
       allocate(m_name%iaux3(m_name%dim1+1))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          k=0
          m_name%iaux3(1)=0
          do i=1,m_name%dim1
             do j=1,m_name%dim2
                k=k+1
                m_name%iaux4(k)=j
                m_name%dval(k,1)=A%dval(i,j)
             end do
             m_name%iaux3(i+1)=k
          end do
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          k=0
          m_name%iaux3(1)=0
          do i=1,m_name%dim1
             do j=1,m_name%dim2
                k=k+1
                m_name%iaux4(k)=j
                m_name%zval(k,1)=A%zval(i,j)
             end do
             m_name%iaux3(i)=k+1
          end do
       end if

  end subroutine m_copy_sdenscsrref

  subroutine m_copy_sdenscsrref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=0
       if (m_name%is_real) then
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%dval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%dim1+1))
             m_name%iaux3_is_allocated=.true.
             allocate(m_name%iaux4(m_name%iaux2(1)))
             m_name%iaux4_is_allocated=.true.
             allocate(m_name%dval(m_name%iaux2(1),1))
             m_name%dval_is_allocated=.true.
             k=0
             m_name%iaux3(1)=0
             do i=1,m_name%dim1
                do j=1,m_name%dim2
                   if (abs(A%dval(i,j))>abs_threshold) then
                      k=k+1
                      m_name%iaux4(k)=j
                      m_name%dval(k,1)=A%dval(i,j)-soft_threshold*A%dval(i,j)/abs(A%dval(i,j))
                   end if
                end do
                m_name%iaux3(i+1)=k
             end do
          end if
       else
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%zval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%dim1+1))
             m_name%iaux3_is_allocated=.true.
             allocate(m_name%iaux4(m_name%iaux2(1)))
             m_name%iaux4_is_allocated=.true.
             allocate(m_name%zval(m_name%iaux2(1),1))
             m_name%zval_is_allocated=.true.
             k=0
             m_name%iaux3(1)=0
             do i=1,m_name%dim1
                do j=1,m_name%dim2
                   if (abs(A%zval(i,j))>abs_threshold) then
                      k=k+1
                      m_name%iaux4(k)=j
                      m_name%zval(k,1)=A%zval(i,j)-soft_threshold*A%zval(i,j)/abs(A%zval(i,j))
                   end if
                end do
                m_name%iaux3(i)=k+1
             end do
          end if
       end if

  end subroutine m_copy_sdenscsrref_thre

  subroutine m_copy_scoosdenref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
          do i=1,A%iaux2(1)
             m_name%dval(A%iaux3(i),A%iaux4(i))=A%dval(i,1)
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
          do i=1,A%iaux2(1)
             m_name%zval(A%iaux3(i),A%iaux4(i))=A%zval(i,1)
          end do
       end if

  end subroutine m_copy_scoosdenref

  subroutine m_copy_scoosdenref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
          do i=1,A%iaux2(1)
             if (abs(A%dval(i,1))>abs_threshold) then
                m_name%dval(A%iaux3(i),A%iaux4(i))=A%dval(i,1)-soft_threshold*A%dval(i,1)/abs(A%dval(i,1))
             else
                m_name%dval(A%iaux3(i),A%iaux4(i))=0.0_dp
             end if
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
          do i=1,A%iaux2(1)
             if (abs(A%zval(i,1))>abs_threshold) then
                m_name%zval(A%iaux3(i),A%iaux4(i))=A%zval(i,1)-soft_threshold*A%zval(i,1)/abs(A%zval(i,1))
             else
                m_name%zval(A%iaux3(i),A%iaux4(i))=cmplx_0
             end if
          end do
       end if

  end subroutine m_copy_scoosdenref_thre

  subroutine m_copy_scscsdenref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
          do i=1,A%dim2
             do j=1,A%iaux4(i+1)-A%iaux4(i)
                k=A%iaux4(i)+j
                m_name%dval(A%iaux3(k),i)=A%dval(k,1)
             end do
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
          do i=1,A%dim2
             do j=1,A%iaux4(i+1)-A%iaux4(i)
                k=A%iaux4(i)+j
                m_name%zval(A%iaux3(k),i)=A%zval(k,1)
             end do
          end do
       end if

  end subroutine m_copy_scscsdenref

  subroutine m_copy_scscsdenref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
          do i=1,A%dim2
             do j=1,A%iaux4(i+1)-A%iaux4(i)
                k=A%iaux4(i)+j
                if (abs(A%dval(k,1))>abs_threshold) then
                   m_name%dval(A%iaux3(k),i)=A%dval(k,1)-soft_threshold*A%dval(k,1)/abs(A%dval(k,1))
                else
                   m_name%dval(A%iaux3(k),i)=0.0_dp
                end if
             end do
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
          do i=1,A%dim2
             do j=1,A%iaux4(i+1)-A%iaux4(i)
                k=A%iaux4(i)+j
                if (abs(A%zval(k,1))>abs_threshold) then
                   m_name%zval(A%iaux3(k),i)=A%zval(k,1)-soft_threshold*A%zval(k,1)/abs(A%zval(k,1))
                else
                   m_name%zval(A%iaux3(k),i)=0.0_dp
                end if
             end do
          end do
       end if

  end subroutine m_copy_scscsdenref_thre

  subroutine m_copy_scsrsdenref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
          do i=1,A%dim1
             do j=1,A%iaux3(i+1)-A%iaux3(i)
                k=A%iaux3(i)+j
                m_name%dval(i,A%iaux4(k))=A%dval(k,1)
             end do
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
          do i=1,A%dim1
             do j=1,A%iaux3(i+1)-A%iaux3(i)
                k=A%iaux3(i)+j
                m_name%zval(i,A%iaux4(k))=A%zval(k,1)
             end do
          end do
       end if

  end subroutine m_copy_scsrsdenref

  subroutine m_copy_scsrsdenref_thre(m_name,A,abs_threshold,soft_threshold)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in) :: abs_threshold, soft_threshold

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
          do i=1,A%dim1
             do j=1,A%iaux3(i+1)-A%iaux3(i)
                k=A%iaux3(i)+j
                if (abs(A%dval(k,1))>abs_threshold) then
                   m_name%dval(i,A%iaux4(k))=A%dval(k,1)-soft_threshold*A%dval(k,1)/abs(A%dval(k,1))
                else
                   m_name%dval(i,A%iaux4(k))=0.0_dp
                end if
             end do
          end do
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
          do i=1,A%dim1
             do j=1,A%iaux3(i+1)-A%iaux3(i)
                k=A%iaux3(i)+j
                if (abs(A%zval(k,1))>abs_threshold) then
                   m_name%zval(i,A%iaux4(k))=A%zval(k,1)-soft_threshold*A%zval(k,1)/abs(A%zval(k,1))
                else
                   m_name%zval(i,A%iaux4(k))=0.0_dp
                end if
             end do
          end do
       end if

  end subroutine m_copy_scsrsdenref_thre

  subroutine m_copy_scscscooref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if
       do i=1,A%dim2
          do j=1,A%iaux4(i+1)-A%iaux4(i)
             k=A%iaux4(i)+j
             m_name%iaux3(k)=A%iaux3(k)
             m_name%iaux4(k)=i
          end do
       end do

  end subroutine m_copy_scscscooref

  subroutine m_copy_scsrscooref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          m_name%zval=A%zval
       end if
       do i=1,A%dim1
          do j=1,A%iaux3(i+1)-A%iaux3(i)
             k=A%iaux3(i)+j
             m_name%iaux3(k)=i
             m_name%iaux4(k)=A%iaux4(k)
          end do
       end do

  end subroutine m_copy_scsrscooref

  subroutine m_copy_scooscscref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k
    integer, allocatable :: sort_temp(:)

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%iaux2(1)))
       m_name%iaux3_is_allocated=.true.
       allocate(m_name%iaux4(m_name%dim2+1))
       m_name%iaux4_is_allocated=.true.
       m_name%iaux4=0
       do i=1,m_name%iaux2(1)
          m_name%iaux4(A%iaux4(i)+1)=m_name%iaux4(A%iaux4(i)+1)+1
       end do
       do i=1,m_name%dim2
          m_name%iaux4(i+1)=m_name%iaux4(i+1)+m_name%iaux4(i)
       end do
       m_name%iaux3=0
       allocate(sort_temp(m_name%dim2))
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          do i=1,m_name%iaux2(1)
             do j=m_name%iaux4(A%iaux4(i))+1,m_name%iaux4(A%iaux4(i)+1)
                if (m_name%iaux3(j)==0) then
                   sort_temp(A%iaux4(i))=j
                   m_name%iaux3(j)=A%iaux3(i)
                   m_name%dval(j,1)=A%dval(i,1)
                   exit
                else if (m_name%iaux3(j)>A%iaux3(i)) then
                   k=sort_temp(A%iaux4(i))
                   m_name%iaux3(j+1:k+1)=m_name%iaux3(j:k)
                   m_name%dval(j+1:k+1,1)=m_name%dval(j:k,1)
                   sort_temp(A%iaux4(i))=k+1
                   m_name%iaux3(j)=A%iaux3(i)
                   m_name%dval(j,1)=A%dval(i,1)
                   exit
                end if
             end do
          end do
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          do i=1,m_name%iaux2(1)
             do j=m_name%iaux4(A%iaux4(i))+1,m_name%iaux4(A%iaux4(i)+1)
                if (m_name%iaux3(j)==0) then
                   sort_temp(A%iaux4(i))=j
                   m_name%iaux3(j)=A%iaux3(i)
                   m_name%zval(j,1)=A%zval(i,1)
                   exit
                else if (m_name%iaux3(j)>A%iaux3(i)) then
                   k=sort_temp(A%iaux4(i))
                   m_name%iaux3(j+1:k+1)=m_name%iaux3(j:k)
                   m_name%zval(j+1:k+1,1)=m_name%zval(j:k,1)
                   sort_temp(A%iaux4(i))=k+1
                   m_name%iaux3(j)=A%iaux3(i)
                   m_name%zval(j,1)=A%zval(i,1)
                   exit
                end if
             end do
          end do
       end if
       deallocate(sort_temp)

  end subroutine m_copy_scooscscref

  subroutine m_copy_scooscsrref(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    type(matrix), intent(inout) :: A

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name

    !**** INTERNAL ********************************!

    integer :: i, j, k
    integer, allocatable :: sort_temp(:)

    !**********************************************!

       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux4_is_allocated=.true.
       allocate(m_name%iaux3(m_name%dim1+1))
       m_name%iaux3_is_allocated=.true.
       m_name%iaux3=0
       do i=1,m_name%iaux2(1)
          m_name%iaux3(A%iaux3(i)+1)=m_name%iaux3(A%iaux3(i)+1)+1
       end do
       do i=1,m_name%dim1
          m_name%iaux3(i+1)=m_name%iaux3(i+1)+m_name%iaux3(i)
       end do
       m_name%iaux4=0
       allocate(sort_temp(m_name%dim1))
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval_is_allocated=.true.
          do i=1,m_name%iaux2(1)
             do j=m_name%iaux3(A%iaux3(i))+1,m_name%iaux3(A%iaux3(i)+1)
                if (m_name%iaux4(j)==0) then
                   sort_temp(A%iaux3(i))=j
                   m_name%iaux4(j)=A%iaux4(i)
                   m_name%dval(j,1)=A%dval(i,1)
                   exit
                else if (m_name%iaux4(j)>A%iaux4(i)) then
                   k=sort_temp(A%iaux3(i))
                   m_name%iaux4(j+1:k+1)=m_name%iaux4(j:k)
                   m_name%dval(j+1:k+1,1)=m_name%dval(j:k,1)
                   sort_temp(A%iaux3(i))=k+1
                   m_name%iaux4(j)=A%iaux4(i)
                   m_name%dval(j,1)=A%dval(i,1)
                   exit
                end if
             end do
          end do
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval_is_allocated=.true.
          do i=1,m_name%iaux2(1)
             do j=m_name%iaux3(A%iaux3(i))+1,m_name%iaux3(A%iaux3(i)+1)
                if (m_name%iaux4(j)==0) then
                   sort_temp(A%iaux3(i))=j
                   m_name%iaux4(j)=A%iaux4(i)
                   m_name%zval(j,1)=A%zval(i,1)
                   exit
                else if (m_name%iaux4(j)>A%iaux4(i)) then
                   k=sort_temp(A%iaux3(i))
                   m_name%iaux4(j+1:k+1)=m_name%iaux4(j:k)
                   m_name%zval(j+1:k+1,1)=m_name%zval(j:k,1)
                   sort_temp(A%iaux3(i))=k+1
                   m_name%iaux4(j)=A%iaux4(i)
                   m_name%zval(j,1)=A%zval(i,1)
                   exit
                end if
             end do
          end do
       end if
       deallocate(sort_temp)

  end subroutine m_copy_scooscsrref

end module ms_m_copy
