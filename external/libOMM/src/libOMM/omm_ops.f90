module omm_ops
use MatrixSwitch
use omm_params
use MatrixSwitch_ops, only: ms_lap_icontxt,&
                            ms_mpi_rank,&
                            ms_lap_nprow,&
                            ms_lap_npcol,&
                            ms_lap_bs_def,&
                            ms_mpi_comm

implicit none

contains

!================================================!
! calculate the gradient of the Kim functional:  !
! G=2*(2*H*C-S*C*HW-H*C*SW)                      !
!================================================!
subroutine calc_G(HW,SW,G,HC,SC,m_operation)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in) :: m_operation

  type(matrix), intent(in) :: HW ! hamiltonian matrix in WF basis
  type(matrix), intent(in) :: SW ! overlap matrix in WF basis

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: G ! gradient of Kim functional
  type(matrix), intent(inout) :: HC ! work matrix
  type(matrix), intent(inout) :: SC ! work matrix

  !**********************************************!

  call m_add(HC,'n',G,4.0_dp,0.0_dp,m_operation)
  call mm_multiply(SW,'n',HC,'n',G,-2.0_dp,1.0_dp,m_operation)
  call mm_multiply(HW,'n',SC,'n',G,-2.0_dp,1.0_dp,m_operation)

end subroutine calc_G

!================================================!
! calculate operator matrix in WF basis          !
!================================================!
subroutine calc_AW(A,C,AW,AC,m_operation)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in) :: m_operation

  type(matrix), intent(in) :: A ! operator matrix in AO basis
  type(matrix), intent(in) :: C ! WF coeffs. matrix

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: AC ! work matrix
  type(matrix), intent(inout) :: AW ! operator matrix in WF basis

  !**********************************************!

  call mm_multiply(C,'n',A,'n',AC,1.0_dp,0.0_dp,m_operation)
  call mm_multiply(AC,'n',C,'c',AW,1.0_dp,0.0_dp,m_operation)

end subroutine calc_AW

!================================================!
! calculate operator matrix in WF basis          !
!================================================!
subroutine calc_HW_callback(H,C,HW,HC,m_operation)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in) :: m_operation

  type(matrix), intent(in) :: C ! WF coeffs. matrix

  !**** CALLBACK ********************************!

  interface
    subroutine H(C,HC)
      use MatrixSwitch
      implicit none
      type(matrix), intent(in) :: C ! WF coeffs. matrix
      type(matrix), intent(inout) :: HC ! work matrix
    end subroutine H
  end interface

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: HC ! work matrix
  type(matrix), intent(inout) :: HW ! operator matrix in WF basis

  !**********************************************!

  call H(C,HC)
  call mm_multiply(HC,'n',C,'c',HW,1.0_dp,0.0_dp,m_operation)

end subroutine calc_HW_callback

!================================================!
! calculate operator matrix in AO basis          !
!================================================!
subroutine calc_A(AW,C1,A,CAW,m_operation)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in) :: m_operation

  type(matrix), intent(in) :: AW ! Operator matrix in WF basis
  type(matrix), intent(in) :: C1 ! WF coeffs. matrix

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: A ! operator matrix in orbital basis
  type(matrix), intent(inout) :: CAW ! work matrix

  !**********************************************!

  call mm_multiply(AW,'n',C1,'n',CAW,1.0_dp,0.0_dp,m_operation)
  call mm_multiply(C1,'c',CAW,'n',A,1.0_dp,0.0_dp,m_operation)

end subroutine calc_A

!================================================!
! calculate operator matrix in AO basis          !
!================================================!
subroutine calc_A2(AW,C1,C2,A,CAW,m_operation)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in) :: m_operation

  type(matrix), intent(in) :: AW ! Operator matrix in WF basis
  type(matrix), intent(in) :: C1 ! WF coeffs. matrix
  type(matrix), intent(in) :: C2 ! pre-multiplied WF coeffs. matrix

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: A ! operator matrix in orbital basis
  type(matrix), intent(inout) :: CAW ! work matrix

  !**********************************************!

  call mm_multiply(AW,'n',C2,'n',CAW,1.0_dp,0.0_dp,m_operation)
  call mm_multiply(C1,'c',CAW,'n',A,1.0_dp,0.0_dp,m_operation)

end subroutine calc_A2

!================================================!
! calculate coeffs. of the quartic line search   !
! equation using analytical expressions          !
!================================================!
subroutine calc_coeff(HW,SW,HWd,SWd,HWdd,SWdd,SWdH,coeff,m_operation)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in) :: m_operation

  type(matrix), intent(in) :: HW ! hamiltonian matrix in WF basis
  type(matrix), intent(in) :: SW ! overlap matrix in WF basis
  type(matrix), intent(in) :: HWd ! g^T*h*c
  type(matrix), intent(in) :: SWd ! g^T*s*c
  type(matrix), intent(in) :: HWdd ! g^T*h*g
  type(matrix), intent(in) :: SWdd ! g^T*h*g

  !**** INOUT ***********************************!

  real(dp), intent(inout) :: coeff(0:4) ! coeffs. of the quartic equation

  type(matrix), intent(inout) :: SWdH ! work matrix

  !**** LOCAL ***********************************!

  real(dp) :: TrH
  real(dp) :: TrHS
  real(dp) :: TrHd
  real(dp) :: TrHdS
  real(dp) :: TrHSd
  real(dp) :: TrHdd
  real(dp) :: TrHddS
  real(dp) :: TrHSdd
  real(dp) :: TrHdSd
  real(dp) :: TrHdSdH
  real(dp) :: TrHddSd
  real(dp) :: TrHdSdd
  real(dp) :: TrHddSdd

  !**********************************************!

  call m_trace(HW,TrH,m_operation)
  call mm_trace(HW,SW,TrHS,m_operation)

  call m_trace(HWd,TrHd,m_operation)
  call mm_trace(SW,HWd,TrHdS,m_operation)
  call mm_trace(HW,SWd,TrHSd,m_operation)

  call m_trace(HWdd,TrHdd,m_operation)
  call mm_trace(HWdd,SW,TrHddS,m_operation)
  call mm_trace(HW,SWdd,TrHSdd,m_operation)
  call m_add(SWd,'c',SWdH,1.0_dp,0.0_dp,m_operation)
  call mm_trace(SWdH,HWd,TrHdSd,m_operation)
  call mm_trace(SWd,HWd,TrHdSdH,m_operation)

  call mm_trace(HWdd,SWd,TrHddSd,m_operation)
  call mm_trace(SWdd,HWd,TrHdSdd,m_operation)

  call mm_trace(HWdd,SWdd,TrHddSdd,m_operation)

  coeff(0)=2.0_dp*TrH-TrHS
  coeff(1)=2.0_dp*(2.0_dp*TrHd-TrHdS-TrHSd)
  coeff(2)=2.0_dp*(TrHdd-TrHdSd-TrHdSdH)-TrHddS-TrHSdd
  coeff(3)=-2.0_dp*(TrHddSd+TrHdSdd)
  coeff(4)=-TrHddSdd

end subroutine calc_coeff

!================================================!
! matrix factorization                           !
! C = U^T*U                                      !
!================================================!
subroutine m_factorize(C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**********************************************!

  if (C%is_real) then
    call m_dfactorize(C,label)
  else
    call m_zfactorize(C,label)
  end if

end subroutine m_factorize

subroutine m_dfactorize(C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot, info

  !**********************************************!

  if (.not. C%is_square) call die('m_dfactorize: matrix C is not square')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'lap') then
          ot=1
        else
          call die('m_dfactorize: invalid implementation')
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        else if (label .eq. 'psp') then
          ot=2
        else if (label .eq. 't1D') then
          ot=2
        else
          call die('m_dfactorize: invalid implementation')
        end if
      end if
  else
    call die('m_dfactorize: invalid implementation')
  end if

  select case (ot)
    case (1)
      call dpotrf('u',C%dim1,C%dval,C%dim1,info)
      if (info/=0) call die('m_dfactorize: error in dpotrf')
    case (2)
      call pdpotrf('u',C%dim1,C%dval,1,1,C%iaux1,info)
      if (info/=0) call die('m_dfactorize: error in pdpotrf')
  end select

end subroutine m_dfactorize

subroutine m_zfactorize(C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot, info

  !**********************************************!

  if (.not. C%is_square) call die('m_zfactorize: matrix C is not square')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'lap') then
          ot=1
        else
          call die('m_zfactorize: invalid implementation')
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        else if (label .eq. 'psp') then
          ot=2
        else if (label .eq. 't1D') then
          ot=2
        else
          call die('m_zfactorize: invalid implementation')
        end if
      end if
  else
    call die('m_zfactorize: invalid implementation')
  end if

  select case (ot)
    case (1)
      call zpotrf('u',C%dim1,C%zval,C%dim1,info)
      if (info/=0) call die('m_zfactorize: error in zpotrf')
    case (2)
      call pzpotrf('u',C%dim1,C%zval,1,1,C%iaux1,info)
      if (info/=0) call die('m_zfactorize: error in pzpotrf')
  end select

end subroutine m_zfactorize

!================================================!
! reduction of generalized eigenproblem to       !
! standard form                                  !
! C := A^(-T)*C*A^(-1)                           !
!================================================!
subroutine m_reduce(A,C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: A ! matrix A
  type(matrix), intent(inout) :: C ! matrix C

  !**********************************************!

  if ((C%is_real) .and. &
      (C%is_real)) then
    call m_dreduce(A,C,label)
  else if ((.not. C%is_real) .and. &
           (.not. C%is_real)) then
    call m_zreduce(A,C,label)
  else
    call die('m_reduce: invalid implementation')
  end if

end subroutine m_reduce

subroutine m_dreduce(A,C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: A ! matrix A
  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  character(1) :: c1, c2
  character(5) :: m_storage

  integer :: ot, i, j, info

  real(dp) :: scale_Cholesky

  type(matrix) :: work1

  !**********************************************!

  if (.not. A%is_square) call die('m_zreduce: matrix A is not square')
  if (.not. C%is_square) call die('m_dreduce: matrix C is not square')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'lap') then
          ot=1
        else
          call die('m_dreduce: invalid implementation')
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        else if (label .eq. 'psp') then
          ot=2
        else if (label .eq. 't1D') then
          ot=2
        else
          call die('m_dreduce: invalid implementation')
        end if
      end if
  else
    call die('m_dreduce: invalid implementation')
  end if

  select case (ot)
    case (1)
      call dsygst(1,'u',C%dim1,C%dval,C%dim1,A%dval,A%dim1,info)
      if (info/=0) call die('m_dreduce: error in dsygst')
      do i=1,C%dim1-1
        do j=i+1,C%dim1
          C%dval(j,i)=C%dval(i,j)
        end do
      end do
    case (2)
      call pdsygst(1,'u',C%dim1,C%dval,1,1,C%iaux1,A%dval,1,1,A%iaux1,scale_Cholesky,info)
      if (info/=0) call die('m_dreduce: error in pdsygst')
      if (C%is_serial) then
        c1='s'
      else
        c1='p'
      end if
      if (C%is_real) then
        c2='d'
      else
        c2='z'
      end if
      ! Fill in missing triangle of C.
      do i = 1, C%dim1
          do j = 1, i-1
              call m_set_element(C, i, j, 0.0_dp, 0.0_dp, label)
          end do
      end do
      write(m_storage,'(a1,a1,a3)') c1, c2, C%str_type
      call m_allocate(work1,C%dim1,C%dim2,m_storage)
      call m_add(C,'t',work1,1.0_dp,0.0_dp,label)
      do i = 1, work1%dim1
          call m_set_element(work1, i, i, 0.0_dp, 0.0_dp, label)
      end do
      call m_add(work1,'n',C,1.0_dp,1.0_dp)
      call m_deallocate(work1)
  end select

end subroutine m_dreduce

subroutine m_zreduce(A,C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: A ! matrix A
  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  character(1) :: c1, c2
  character(5) :: m_storage

  integer :: ot, i, j, info

  real(dp) :: scale_Cholesky

  type(matrix) :: work1

  !**********************************************!

  if (.not. A%is_square) call die('m_zreduce: matrix A is not square')
  if (.not. C%is_square) call die('m_zreduce: matrix C is not square')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'lap') then
          ot=1
        else
          call die('m_zreduce: invalid implementation')
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        else if (label .eq. 'psp') then
          ot=2
        else if (label .eq. 't1D') then
          ot=2
        else
          call die('m_zreduce: invalid implementation')
        end if
      end if
  else
    call die('m_zreduce: invalid implementation')
  end if

  select case (ot)
    case (1)
      call zhegst(1,'u',C%dim1,C%zval,C%dim1,A%zval,A%dim1,info)
      if (info/=0) call die('m_zreduce: error in zhegst')
      do i=1,C%dim1-1
        do j=i+1,C%dim1
          C%zval(j,i)=conjg(C%zval(i,j))
        end do
      end do
    case (2)
      call pzhegst(1,'u',C%dim1,C%zval,1,1,C%iaux1,A%zval,1,1,A%iaux1,scale_Cholesky,info)
      if (info/=0) call die('m_zreduce: error in pzhegst')
      if (C%is_serial) then
        c1='s'
      else
        c1='p'
      end if
      if (C%is_real) then
        c2='d'
      else
        c2='z'
      end if
      ! Fill in missing triangle of C.
      do i = 1, C%dim1
          do j = 1, i-1
              call m_set_element(C, i, j, 0.0_dp, 0.0_dp, label)
          end do
      end do
      write(m_storage,'(a1,a1,a3)') c1, c2, C%str_type
      call m_allocate(work1,C%dim1,C%dim2,m_storage)
      call m_add(C,'c',work1,1.0_dp,0.0_dp,label)
      do i = 1, work1%dim1
          call m_set_element(work1, i, i, 0.0_dp, 0.0_dp, label)
      end do
      call m_add(work1,'n',C,1.0_dp,1.0_dp)
      call m_deallocate(work1)
  end select

end subroutine m_zreduce

!================================================!
! transformation of coefficients                 !
!================================================!
subroutine m_transform(A,C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: A ! matrix A
  type(matrix), intent(inout) :: C ! matrix C

  !**********************************************!

  if ((C%is_real) .and. &
      (C%is_real)) then
    call m_dtransform(A,C,label)
  else if ((.not. C%is_real) .and. &
           (.not. C%is_real)) then
    call m_ztransform(A,C,label)
  else
    call die('m_transform: invalid implementation')
  end if

end subroutine m_transform

subroutine m_dtransform(A,C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: A ! matrix A
  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot

  !**********************************************!

  if (.not. A%is_square) call die('m_dtransform: matrix C is not square')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'lap') then
          ot=1
        else
          call die('m_dtransform: invalid implementation')
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        else if (label .eq. 'psp') then
          ot=2
        else if (label .eq. 't1D') then
          ot=2
        else
          call die('m_dtransform: invalid implementation')
        end if
      end if
  else
    call die('m_dtransform: invalid implementation')
  end if

  select case (ot)
    case (1)
      call dtrmm('r','u','t','n',C%dim1,C%dim2,1.0_dp,A%dval,A%dim1,C%dval,C%dim1)
    case (2)
      call pdtrmm('r','u','t','n',C%dim1,C%dim2,1.0_dp,A%dval,1,1,A%iaux1,C%dval,1,1,C%iaux1)
  end select

end subroutine m_dtransform

subroutine m_ztransform(A,C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: A ! matrix A
  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot

  !**********************************************!

  if (.not. A%is_square) call die('m_ztransform: matrix C is not square')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'lap') then
          ot=1
        else
          call die('m_ztransform: invalid implementation')
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        else if (label .eq. 'psp') then
          ot=2
        else if (label .eq. 't1D') then
          ot=2
        else
          call die('m_ztransform: invalid implementation')
        end if
      end if
  else
    call die('m_ztransform: invalid implementation')
  end if

  select case (ot)
    case (1)
      call ztrmm('r','u','c','n',C%dim1,C%dim2,cmplx_1,A%zval,A%dim1,C%zval,C%dim1)
    case (2)
      call pztrmm('r','u','c','n',C%dim1,C%dim2,cmplx_1,A%zval,1,1,A%iaux1,C%zval,1,1,C%iaux1)
  end select

end subroutine m_ztransform

!================================================!
! recovery of coefficients                       !
!================================================!
subroutine m_back_transform(A,C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: A ! matrix A
  type(matrix), intent(inout) :: C ! matrix C

  !**********************************************!

  if ((C%is_real) .and. &
      (C%is_real)) then
    call m_dback_transform(A,C,label)
  else if ((.not. C%is_real) .and. &
           (.not. C%is_real)) then
    call m_zback_transform(A,C,label)
  else
    call die('m_back_transform: invalid implementation')
  end if

end subroutine m_back_transform

subroutine m_dback_transform(A,C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: A ! matrix A
  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot

  !**********************************************!

  if (.not. A%is_square) call die('m_dback_transform: matrix C is not square')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'lap') then
          ot=1
        else
          call die('m_dback_transform: invalid implementation')
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        else if (label .eq. 'psp') then
          ot=2
        else if (label .eq. 't1D') then
          ot=2
        else
          call die('m_dback_transform: invalid implementation')
        end if
      end if
  else
    call die('m_dback_transform: invalid implementation')
  end if

  select case (ot)
    case (1)
      call dtrsm('r','u','t','n',C%dim1,C%dim2,1.0_dp,A%dval,A%dim1,C%dval,C%dim1)
    case (2)
      call pdtrsm('r','u','t','n',C%dim1,C%dim2,1.0_dp,A%dval,1,1,A%iaux1,C%dval,1,1,C%iaux1)
  end select

end subroutine m_dback_transform

subroutine m_zback_transform(A,C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: A ! matrix A
  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot

  !**********************************************!

  if (.not. A%is_square) call die('m_zback_transform: matrix C is not square')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'lap') then
          ot=1
        else
          call die('m_zback_transform: invalid implementation')
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        else if (label .eq. 'psp') then
          ot=2
        else if (label .eq. 't1D') then
          ot=2
        else
          call die('m_zback_transform: invalid implementation')
        end if
      end if
  else
    call die('m_zback_transform: invalid implementation')
  end if

  select case (ot)
    case (1)
      call ztrsm('r','u','c','n',C%dim1,C%dim2,cmplx_1,A%zval,A%dim1,C%zval,C%dim1)
    case (2)
      call pztrsm('r','u','c','n',C%dim1,C%dim2,cmplx_1,A%zval,1,1,A%iaux1,C%zval,1,1,C%iaux1)
  end select

end subroutine m_zback_transform

!================================================!
! matrix inverse                                 !
! C := C^(-1)                                    !
!================================================!
subroutine m_inverse(C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**********************************************!

  if (C%is_real) then
    call m_dinverse(C,label)
  else
    call m_zinverse(C,label)
  end if

end subroutine m_inverse

subroutine m_dinverse(C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot, lwork, i, j, info
  integer, allocatable :: ipiv(:)
  integer :: liwork
  integer, allocatable :: iwork(:)

  real(dp), allocatable :: work(:)

  !**** EXTERNAL ********************************!

  integer, external :: ilaenv

  !**********************************************!

  if (.not. C%is_square) call die('m_dinverse: matrix C is not square')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'lap') then
          ot=1
        else if (label .eq. 'ref') then
          ot=1
        else
          call die('m_dinverse: invalid implementation')
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        else if (label .eq. 'psp') then
          ot=2
        else if (label .eq. 't1D') then
          ot=2
        else
          call die('m_dinverse: invalid implementation')
        end if
      end if
  else
    call die('m_dinverse: invalid implementation')
  end if

  select case (ot)
    case (1)
      allocate(ipiv(C%dim1))
      lwork=C%dim1*ilaenv(1,'dsytrf','u',C%dim1,-1,-1,-1)
      allocate(work(lwork))
      call dsytrf('u',C%dim1,C%dval,C%dim1,ipiv,work,lwork,info)
      if (info/=0) call die('m_dinverse: error in dsytrf')
      deallocate(work)
      allocate(work(C%dim1))
      call dsytri('u',C%dim1,C%dval,C%dim1,ipiv,work,info)
      if (info/=0) call die('m_dinverse: error in dsytri')
      deallocate(work)
      deallocate(ipiv)
      do i=1,C%dim1-1
        do j=i+1,C%dim1
          C%dval(j,i)=C%dval(i,j)
        end do
      end do
    case (2)
      allocate(ipiv(C%iaux1(3)+C%iaux1(5)))
      call pdgetrf(C%dim1,C%dim2,C%dval,1,1,C%iaux1,ipiv,info)
      if (info/=0) call die('m_dinverse: error in pdgetrf')
      allocate(work(1))
      allocate(iwork(1))
      call pdgetri(C%dim1,C%dval,1,1,C%iaux1,ipiv,work,-1,iwork,-1,info)
      if (info/=0) call die('m_dinverse: error in pdgetri')
      liwork=iwork(1)
      deallocate(iwork)
      lwork=int(work(1))
      deallocate(work)
      allocate(work(lwork))
      allocate(iwork(liwork))
      call pdgetri(C%dim1,C%dval,1,1,C%iaux1,ipiv,work,lwork,iwork,liwork,info)
      if (info/=0) call die('m_dinverse: error in pdgetri')
      deallocate(iwork)
      deallocate(work)
      deallocate(ipiv)
  end select

end subroutine m_dinverse

subroutine m_zinverse(C,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot, lwork, i, j, info
  integer, allocatable :: ipiv(:)
  integer :: liwork
  integer, allocatable :: iwork(:)

  complex(dp), allocatable :: work(:)

  !**** EXTERNAL ********************************!

  integer, external :: ilaenv

  !**********************************************!

  if (.not. C%is_square) call die('m_zinverse: matrix C is not square')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'lap') then
          ot=1
        else if (label .eq. 'ref') then
          ot=1
        else
          call die('m_zinverse: invalid implementation')
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        else if (label .eq. 'psp') then
          ot=2
        else if (label .eq. 't1D') then
          ot=2
        else
          call die('m_zinverse: invalid implementation')
        end if
      end if
  else
    call die('m_zinverse: invalid implementation')
  end if

  select case (ot)
    case (1)
      allocate(ipiv(C%dim1))
      lwork=C%dim1*ilaenv(1,'zhetrf','u',C%dim1,-1,-1,-1)
      allocate(work(lwork))
      call zhetrf('u',C%dim1,C%zval,C%dim1,ipiv,work,lwork,info)
      if (info/=0) call die('m_zinverse: error in zhetrf')
      deallocate(work)
      allocate(work(C%dim1))
      call zhetri('u',C%dim1,C%zval,C%dim1,ipiv,work,info)
      if (info/=0) call die('m_zinverse: error in zhetri')
      deallocate(work)
      deallocate(ipiv)
      do i=1,C%dim1-1
        do j=i+1,C%dim1
          C%zval(j,i)=conjg(C%zval(i,j))
        end do
      end do
    case (2)
      allocate(ipiv(C%iaux1(3)+C%iaux1(5)))
      call pzgetrf(C%dim1,C%dim2,C%zval,1,1,C%iaux1,ipiv,info)
      if (info/=0) call die('m_zinverse: error in pzgetrf')
      allocate(work(1))
      allocate(iwork(1))
      call pzgetri(C%dim1,C%zval,1,1,C%iaux1,ipiv,work,-1,iwork,-1,info)
      if (info/=0) call die('m_zinverse: error in pzgetri')
      liwork=iwork(1)
      deallocate(iwork)
      lwork=int(work(1))
      deallocate(work)
      allocate(work(lwork))
      allocate(iwork(liwork))
      call pzgetri(C%dim1,C%zval,1,1,C%iaux1,ipiv,work,lwork,iwork,liwork,info)
      if (info/=0) call die('m_zinverse: error in pzgetri')
      deallocate(iwork)
      deallocate(work)
      deallocate(ipiv)
  end select

end subroutine m_zinverse

!================================================!
! plane-wave preconditioner                      !
!================================================!
subroutine calc_PW_precon(T,scale_T,P)
  implicit none

  !**** INPUT ***********************************!

  real(dp), intent(in) :: scale_T ! kinetic energy scale for the preconditioning

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: T ! kinetic energy matrix
  type(matrix), intent(inout) :: P ! preconditioner matrix

  !**********************************************!

  if (T%is_real) then
    call dcalc_PW_precon(T,scale_T,P)
  else
    call zcalc_PW_precon(T,scale_T,P)
  end if

end subroutine calc_PW_precon

subroutine dcalc_PW_precon(T,scale_T,P)
  implicit none

  !**** INPUT ***********************************!

  real(dp), intent(in) :: scale_T ! kinetic energy scale for the preconditioning

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: T ! kinetic energy matrix
  type(matrix), intent(inout) :: P ! preconditioner matrix

  !**** INTERNAL ********************************!

  integer :: ot

  !**********************************************!

  ! operation table
  if (T%str_type .eq. 'csc') then
    ot=1
  else if (T%str_type .eq. 'csr') then
    ot=1
  else
    call die('dcalc_PW_precon: invalid implementation')
  end if

  select case (ot)
    case (1)
      call die('dcalc_PW_precon: compile with pspBLAS')
  end select

end subroutine dcalc_PW_precon

subroutine zcalc_PW_precon(T,scale_T,P)
  implicit none

  !**** INPUT ***********************************!

  real(dp), intent(in) :: scale_T ! kinetic energy scale for the preconditioning

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: T ! kinetic energy matrix
  type(matrix), intent(inout) :: P ! preconditioner matrix

  !**** INTERNAL ********************************!

  integer :: ot

  !**********************************************!

  ! operation table
  if (T%str_type .eq. 'csc') then
    ot=1
  else if (T%str_type .eq. 'csr') then
    ot=1
  else
    call die('zcalc_PW_precon: invalid implementation')
  end if

  select case (ot)
    case (1)
      call die('zcalc_PW_precon: compile with pspBLAS')
  end select

end subroutine zcalc_PW_precon

subroutine die(message)
  implicit none

  !**** INPUT ***********************************!

  character(*), intent(in), optional :: message

  !**** INTERNAL ********************************!

  integer :: err_unit=377

  !**********************************************!

  open(unit=err_unit,file='libOMM.err')
  write(err_unit,'(a)') 'FATAL ERROR in libOMM!'
  if (present(message)) write(err_unit,'(a)') message
  write(err_unit,'(a,1x,i5)') 'MPI rank:', ms_mpi_rank
  close(err_unit)
  stop

end subroutine die

end module omm_ops
