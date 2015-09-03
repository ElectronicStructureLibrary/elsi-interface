module MatrixSwitch
implicit none

private

!**** PARAMS ************************************!

integer, parameter :: dp=selected_real_kind(15,300)

complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

!**** VARIABLES *********************************!

#ifdef MPI
character(1), save :: ms_lap_order

integer, save :: ms_mpi_size
integer, save :: ms_lap_nprow
integer, save :: ms_lap_npcol
integer, save :: ms_lap_bs_def
integer, save :: ms_lap_bs_num
integer, save :: ms_lap_icontxt ! BLACS context handle used by MatrixSwitch
integer, allocatable, save :: ms_lap_bs_list(:,:)
#endif

!**** TYPES *************************************!

! This is the derived type that encapsulates all matrix storage possibilities and hides the details from the user.
type matrix
  character(3) :: str_type ! label identifying the storage format

  logical :: is_initialized=.false. ! has the matrix been initialized?
  logical :: is_serial ! is the matrix serial or parallel distributed?
  logical :: is_real ! is the matrix real or complex (both kind dp)?
  logical :: is_square ! is the matrix square?

  integer :: dim1 ! (global) row dimension size of the matrix
  integer :: dim2 ! (global) column dimension size of the matrix
  integer, pointer :: iaux1(:) => null() ! auxiliary information for certain storage formats
  integer, pointer :: iaux2(:) => null() ! auxiliary information for certain storage formats

  real(dp), pointer :: dval(:,:) => null() ! matrix elements for a real matrix

  complex(dp), pointer :: zval(:,:) => null() ! matrix elements for a complex matrix
end type matrix

!**** INTERFACES ********************************!

interface mm_multiply
  module procedure mm_dmultiply
  module procedure mm_zmultiply
end interface mm_multiply

interface m_add
  module procedure m_dadd
  module procedure m_zadd
end interface m_add

interface m_trace
  module procedure m_dtrace
  module procedure m_ztrace
end interface m_trace

interface mm_trace
  module procedure mm_dtrace
  module procedure mm_ztrace
end interface mm_trace

interface m_scale
  module procedure m_dscale
  module procedure m_zscale
end interface m_scale

interface m_set
  module procedure m_dset
  module procedure m_zset
end interface m_set

interface m_set_element
  module procedure m_dset_element
  module procedure m_zset_element
end interface m_set_element

interface m_get_element
  module procedure m_dget_element
  module procedure m_zget_element
end interface m_get_element

interface m_register_sden
  module procedure m_register_sdden
  module procedure m_register_szden
end interface m_register_sden

interface process_opM
  module procedure process_lopM
  module procedure process_iopM
end interface process_opM

#ifdef MPI
interface m_register_pdbc
  module procedure m_register_pddbc
  module procedure m_register_pzdbc
end interface m_register_pdbc
#endif

!************************************************!

public :: matrix
public :: mm_multiply
public :: m_add
public :: m_trace
public :: mm_trace
public :: m_scale
public :: m_set
public :: m_set_element
public :: m_get_element
public :: m_register_sden
public :: m_allocate
public :: m_deallocate
#ifdef MPI
public :: m_register_pdbc
public :: ms_scalapack_setup
public :: ms_lap_icontxt
#endif

contains

!================================================!
! allocate matrix                                !
!================================================!
subroutine m_allocate(m_name,i,j,label)
  implicit none

  !**** INPUT ***********************************!

  character(5), intent(in), optional :: label ! storage format to use (see documentation)

  integer, intent(in) :: i ! (global) row dimension size of the matrix
  integer, intent(in) :: j ! (global) column dimension size of the matrix

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: m_name ! matrix to be allocated

  !**** INTERNAL ********************************!

  character(1) :: c1, c2

  integer :: st

  !**********************************************!

  m_name%dim1=i
  m_name%dim2=j
  if (i==j) then
    m_name%is_square=.true.
  else
    m_name%is_square=.false.
  end if

  if (present(label)) then
    read(label,'(a1,a1,a3)') c1, c2, m_name%str_type
  else
    c1='s'
    c2='d'
    m_name%str_type='den'
  end if

  if (c1 .eq. 's') then
    m_name%is_serial=.true.
  else if (c1 .eq. 'p') then
    m_name%is_serial=.false.
#ifndef MPI
    call die('m_allocate: compile with MPI')
#endif
  else
    call die('m_allocate: invalid label')
  end if
  if (c2 .eq. 'd') then
    m_name%is_real=.true.
  else if (c2 .eq. 'z') then
    m_name%is_real=.false.
  else
    call die('m_allocate: invalid label')
  end if

  ! storage type
  if ((m_name%str_type .eq. 'den') .and. &
      (m_name%is_serial)) then
      st=1
  else if ((m_name%str_type .eq. 'dbc') .and. &
      (.not. m_name%is_serial)) then
      st=2
  else
    call die('m_allocate: invalid label')
  end if

  select case (st)
    case (1)
      if (m_name%is_real) then
        allocate(m_name%dval(m_name%dim1,m_name%dim2))
        m_name%dval=0.0_dp
      else
        allocate(m_name%zval(m_name%dim1,m_name%dim2))
        m_name%zval=cmplx_0
      end if
    case (2)
#if defined(MPI) && defined(SLAP)
      call ms_scalapack_allocate(m_name)
#else
      call die('m_allocate: compile with ScaLAPACK')
#endif
  end select

  m_name%is_initialized=.true.

end subroutine m_allocate

!================================================!
! deallocate matrix                              !
!================================================!
subroutine m_deallocate(m_name)
  implicit none

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: m_name ! matrix to be deallocated

  !**********************************************!

  if (associated(m_name%iaux1)) nullify(m_name%iaux1)
  if (associated(m_name%iaux2)) nullify(m_name%iaux2)
  if (associated(m_name%dval)) nullify(m_name%dval)
  if (associated(m_name%zval)) nullify(m_name%zval)

  m_name%is_initialized=.false.

end subroutine m_deallocate

!================================================!
! matrix-matrix multiplication                   !
! C := alpha*op(A)*op(B) + beta*C, where         !
! op(M) can be M^T (transpose) or                !
!              M^H (Hermitian transpose)         !
!================================================!
subroutine mm_dmultiply(A,opA,B,opB,C,alpha,beta,label)
  implicit none

  !**** INPUT ***********************************!

  character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
  character(1), intent(in) :: opB ! form of op(B)
  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  real(dp), intent(in) :: alpha ! scalar alpha
  real(dp), intent(in) :: beta ! scalar beta

  type(matrix), intent(in) :: A ! matrix A
  type(matrix), intent(in) :: B ! matrix B

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  logical :: trA, trB

  integer :: ot, i

#ifdef CONV
  complex(dp) :: cmplx_alpha, cmplx_beta
#endif

  !**********************************************!

#ifdef CONV
  if ((.not. A%is_real) .and. (.not. B%is_real) .and. (.not. C%is_real)) then
    cmplx_alpha=cmplx(alpha,0.0_dp,dp)
    cmplx_beta=cmplx(beta,0.0_dp,dp)
    call mm_zmultiply(A,opA,B,opB,C,cmplx_alpha,cmplx_beta,label)
    return
  end if
#endif
  if (.not. A%is_real) call die('mm_dmultiply: matrix A is complex')
  if (.not. B%is_real) call die('mm_dmultiply: matrix B is complex')
  if (.not. C%is_real) call die('mm_dmultiply: matrix C is complex')
  call process_opM(opA,trA)
  call process_opM(opB,trB)
  if ((.not. trA) .and. (.not. trB)) then
    if ((A%dim1/=C%dim1) .or. &
        (B%dim2/=C%dim2) .or. &
        (A%dim2/=B%dim1)) call die('mm_dmultiply: matrices A, B and C are not compatible')
  else if (trA .and. (.not. trB)) then
    if ((A%dim2/=C%dim1) .or. &
        (B%dim2/=C%dim2) .or. &
        (A%dim1/=B%dim1)) call die('mm_dmultiply: matrices A, B and C are not compatible')
  else if ((.not. trA) .and. trB) then
    if ((A%dim1/=C%dim1) .or. &
        (B%dim1/=C%dim2) .or. &
        (A%dim2/=B%dim2)) call die('mm_dmultiply: matrices A, B and C are not compatible')
  else if (trA .and. trB) then
    if ((A%dim2/=C%dim1) .or. &
        (B%dim1/=C%dim2) .or. &
        (A%dim1/=B%dim2)) call die('mm_dmultiply: matrices A, B and C are not compatible')
  end if

  ! operation table
  if ((A%str_type .eq. 'den') .and. &
      (A%is_serial) .and. &
      (B%str_type .eq. 'den') .and. &
      (B%is_serial) .and. &
      (C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else if ((A%str_type .eq. 'dbc') .and. &
           (.not. A%is_serial) .and. &
           (B%str_type .eq. 'dbc') .and. &
           (.not. B%is_serial) .and. &
           (C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=3
      else
        if (label .eq. 'lap') then
          ot=3
        end if
      end if
  else
    call die('mm_dmultiply: invalid implementation')
  end if

  select case (ot)
    case (1)
      call mm_multiply_sddenref(A,trA,B,trB,C,alpha,beta)
    case (2)
#ifdef LAP
      if (trA) then
        i=A%dim1
      else
        i=A%dim2
      end if
      call dgemm(opA,opB,C%dim1,C%dim2,i,alpha,A%dval,A%dim1,B%dval,B%dim1,beta,C%dval,C%dim1)
#else
      call die('mm_dmultiply: compile with LAPACK')
#endif
    case (3)
#ifdef SLAP
      if (trA) then
        i=A%dim1
      else
        i=A%dim2
      end if
      call pdgemm(opA,opB,C%dim1,C%dim2,i,alpha,A%dval,1,1,A%iaux1,B%dval,1,1,B%iaux1,beta,C%dval,1,1,C%iaux1)
#else
      call die('mm_dmultiply: compile with ScaLAPACK')
#endif
  end select

end subroutine mm_dmultiply

subroutine mm_zmultiply(A,opA,B,opB,C,alpha,beta,label)
  implicit none

  !**** INPUT ***********************************!

  character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T' for A^T, 'c/C' for A^H
  character(1), intent(in) :: opB ! form of op(B)
  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  complex(dp), intent(in) :: alpha ! scalar alpha
  complex(dp), intent(in) :: beta ! scalar beta

  type(matrix), intent(in) :: A ! matrix A
  type(matrix), intent(in) :: B ! matrix B

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: tcA, tcB, ot, i

#ifdef CONV
  real(dp) :: real_alpha, real_beta
#endif

  !**********************************************!

#ifdef CONV
  if (A%is_real .and. B%is_real .and. C%is_real) then
    real_alpha=real(alpha,dp)
    real_beta=real(beta,dp)
    call mm_dmultiply(A,opA,B,opB,C,real_alpha,real_beta,label)
    return
  end if
#endif
  if (A%is_real) call die('mm_zmultiply: matrix A is real')
  if (B%is_real) call die('mm_zmultiply: matrix B is real')
  if (C%is_real) call die('mm_zmultiply: matrix C is real')
  call process_opM(opA,tcA)
  call process_opM(opB,tcB)
  if ((tcA==0) .and. (tcB==0)) then
    if ((A%dim1/=C%dim1) .or. &
        (B%dim2/=C%dim2) .or. &
        (A%dim2/=B%dim1)) call die('mm_zmultiply: matrices A, B and C are not compatible')
  else if ((tcA>0) .and. (tcB==0)) then
    if ((A%dim2/=C%dim1) .or. &
        (B%dim2/=C%dim2) .or. &
        (A%dim1/=B%dim1)) call die('mm_zmultiply: matrices A, B and C are not compatible')
  else if ((tcA==0) .and. (tcB>0)) then
    if ((A%dim1/=C%dim1) .or. &
        (B%dim1/=C%dim2) .or. &
        (A%dim2/=B%dim2)) call die('mm_zmultiply: matrices A, B and C are not compatible')
  else if ((tcA>0) .and. (tcB>0)) then
    if ((A%dim2/=C%dim1) .or. &
        (B%dim1/=C%dim2) .or. &
        (A%dim1/=B%dim2)) call die('mm_zmultiply: matrices A, B and C are not compatible')
  end if

  ! operation table
  if ((A%str_type .eq. 'den') .and. &
      (A%is_serial) .and. &
      (B%str_type .eq. 'den') .and. &
      (B%is_serial) .and. &
      (C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else if ((A%str_type .eq. 'dbc') .and. &
           (.not. A%is_serial) .and. &
           (B%str_type .eq. 'dbc') .and. &
           (.not. B%is_serial) .and. &
           (C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=3
      else
        if (label .eq. 'lap') then
          ot=3
        end if
      end if
  else
    call die('mm_zmultiply: invalid implementation')
  end if

  select case (ot)
    case (1)
      call mm_multiply_szdenref(A,tcA,B,tcB,C,alpha,beta)
    case (2)
#ifdef LAP
      if (tcA>0) then
        i=A%dim1
      else
        i=A%dim2
      end if
      call zgemm(opA,opB,C%dim1,C%dim2,i,alpha,A%zval,A%dim1,B%zval,B%dim1,beta,C%zval,C%dim1)
#else
      call die('mm_zmultiply: compile with LAPACK')
#endif
    case (3)
#ifdef SLAP
      if (tcA>0) then
        i=A%dim1
      else
        i=A%dim2
      end if
      call pzgemm(opA,opB,C%dim1,C%dim2,i,alpha,A%zval,1,1,A%iaux1,B%zval,1,1,B%iaux1,beta,C%zval,1,1,C%iaux1)
#else
      call die('mm_zmultiply: compile with ScaLAPACK')
#endif
  end select

end subroutine mm_zmultiply

!================================================!
! matrix addition                                !
! C := alpha*op(A) + beta*C, where               !
! op(M) can be M^T or M^H                        !
!================================================!
subroutine m_dadd(A,opA,C,alpha,beta,label)
  implicit none

  !**** INPUT ***********************************!

  character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  real(dp), intent(in) :: alpha ! scalar alpha
  real(dp), intent(in) :: beta ! scalar beta

  type(matrix), intent(in) :: A ! matrix A

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  logical :: trA

  integer :: ot

#ifdef CONV
  complex(dp) :: cmplx_alpha, cmplx_beta
#endif

  !**********************************************!

#ifdef CONV
  if ((.not. A%is_real) .and. (.not. C%is_real)) then
    cmplx_alpha=cmplx(alpha,0.0_dp,dp)
    cmplx_beta=cmplx(beta,0.0_dp,dp)
    call m_zadd(A,opA,C,cmplx_alpha,cmplx_beta,label)
    return
  end if
#endif
  if (.not. A%is_real) call die('m_dadd: matrix A is complex')
  if (.not. C%is_real) call die('m_dadd: matrix C is complex')
  call process_opM(opA,trA)
  if (.not. (trA)) then
    if ((A%dim1/=C%dim1) .or. &
        (A%dim2/=C%dim2)) call die('m_dadd: matrices A and C are not compatible')
  else
    if ((A%dim1/=C%dim2) .or. &
        (A%dim2/=C%dim1)) call die('m_dadd: matrices A and C are not compatible')
  end if

  ! operation table
  if ((A%str_type .eq. 'den') .and. &
      (A%is_serial) .and. &
      (C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((A%str_type .eq. 'dbc') .and. &
           (.not. A%is_serial) .and. &
           (C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else
    call die('m_dadd: invalid implementation')
  end if

  select case (ot)
    case (1)
      call m_add_sddenref(A,trA,C,alpha,beta)
    case (2)
#ifdef SLAP
      call pdgeadd(opA,C%dim1,C%dim2,alpha,A%dval,1,1,A%iaux1,beta,C%dval,1,1,C%iaux1)
#else
      call die('m_dadd: compile with ScaLAPACK')
#endif
  end select

end subroutine m_dadd

subroutine m_zadd(A,opA,C,alpha,beta,label)
  implicit none

  !**** INPUT ***********************************!

  character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T' for A^T, 'c/C' for A^H
  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  complex(dp), intent(in) :: alpha ! scalar alpha
  complex(dp), intent(in) :: beta ! scalar beta

  type(matrix), intent(in) :: A ! matrix A

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: tcA, ot

#ifdef CONV
  real(dp) :: real_alpha, real_beta
#endif

  !**********************************************!

#ifdef CONV
  if (A%is_real .and. C%is_real) then
    real_alpha=real(alpha,dp)
    real_beta=real(beta,dp)
    call m_dadd(A,opA,C,real_alpha,real_beta,label)
    return
  end if
#endif
  if (A%is_real) call die('m_zadd: matrix A is real')
  if (C%is_real) call die('m_zadd: matrix C is real')
  call process_opM(opA,tcA)
  if (tcA==0) then
    if ((A%dim1/=C%dim1) .or. &
        (A%dim2/=C%dim2)) call die('m_dadd: matrices A and C are not compatible')
  else
    if ((A%dim1/=C%dim2) .or. &
        (A%dim2/=C%dim1)) call die('m_dadd: matrices A and C are not compatible')
  end if

  ! operation table
  if ((A%str_type .eq. 'den') .and. &
      (A%is_serial) .and. &
      (C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((A%str_type .eq. 'dbc') .and. &
           (.not. A%is_serial) .and. &
           (C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else
    call die('m_zadd: invalid implementation')
  end if

  select case (ot)
    case (1)
      call m_add_szdenref(A,tcA,C,alpha,beta)
    case (2)
#ifdef SLAP
      call pzgeadd(opA,C%dim1,C%dim2,alpha,A%zval,1,1,A%iaux1,beta,C%zval,1,1,C%iaux1)
#else
      call die('m_zadd: compile with ScaLAPACK')
#endif
  end select

end subroutine m_zadd

!================================================!
! matrix trace                                   !
! alpha := tr(A)                                 !
!================================================!
subroutine m_dtrace(A,alpha,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  type(matrix), intent(in) :: A ! matrix A

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: alpha ! scalar alpha

  !**** INTERNAL ********************************!

  integer :: ot, i

#ifdef CONV
  complex(dp) :: cmplx_alpha
#endif

  !**** EXTERNAL ********************************!

#if defined(MPI) && defined(SLAP)
  real(dp), external :: pdlatra
#endif

  !**********************************************!

#ifdef CONV
  if (.not. A%is_real) then
    call m_ztrace(A,cmplx_alpha,label)
    alpha=real(cmplx_alpha,dp)
    return
  end if
#else
  if (.not. A%is_real) call die('m_dtrace: matrix A is complex')
#endif
  if (.not. A%is_square) call die('m_dtrace: matrix A is not square')

  ! operation table
  if ((A%str_type .eq. 'den') .and. &
      (A%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((A%str_type .eq. 'dbc') .and. &
           (.not. A%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else
    call die('m_dtrace: invalid implementation')
  end if

  select case (ot)
    case (1)
      alpha=0.0_dp
      do i=1,A%dim1
        alpha=alpha+A%dval(i,i)
      end do
    case (2)
#ifdef SLAP
      alpha=pdlatra(A%dim1,A%dval,1,1,A%iaux1)
#else
      call die('m_dtrace: compile with ScaLAPACK')
#endif
  end select

end subroutine m_dtrace

subroutine m_ztrace(A,alpha,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  type(matrix), intent(in) :: A ! matrix A

  !**** OUTPUT **********************************!

  complex(dp), intent(out) :: alpha ! scalar alpha

  !**** INTERNAL ********************************!

  integer :: ot, i

#ifdef CONV
  real(dp) :: real_alpha
#endif

  !**** EXTERNAL ********************************!

#if defined(MPI) && defined(SLAP)
  real(dp), external :: pzlatra
#endif

  !**********************************************!

#ifdef CONV
  if (A%is_real) then
    call m_dtrace(A,real_alpha,label)
    alpha=cmplx(real_alpha,0.0_dp,dp)
    return
  end if
#else
  if (A%is_real) call die('m_ztrace: matrix A is real')
#endif
  if (.not. A%is_square) call die('m_ztrace: matrix A is not square')

  ! operation table
  if ((A%str_type .eq. 'den') .and. &
      (A%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((A%str_type .eq. 'dbc') .and. &
           (.not. A%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else
    call die('m_ztrace: invalid implementation')
  end if

  select case (ot)
    case (1)
      alpha=cmplx_0
      do i=1,A%dim1
        alpha=alpha+A%zval(i,i)
      end do
    case (2)
#ifdef SLAP
      alpha=pzlatra(A%dim1,A%zval,1,1,A%iaux1)
#else
      call die('m_ztrace: compile with ScaLAPACK')
#endif
  end select

end subroutine m_ztrace

!================================================!
! matrix product trace                           !
! alpha := tr(A^H*B) = tr(B*A^H)                 !
!================================================!
subroutine mm_dtrace(A,B,alpha,label)
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  type(matrix), intent(in) :: A ! matrix A
  type(matrix), intent(in) :: B ! matrix B

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: alpha ! scalar alpha

  !**** INTERNAL ********************************!

  integer :: ot, i, j

#ifdef CONV
  complex(dp) :: cmplx_alpha
#endif

  !**** EXTERNAL ********************************!

#ifdef MPI
  integer :: info

  real(dp) :: alpha_loc
#endif
#ifdef LAP
  real(dp), external :: ddot
#endif

  !**********************************************!

#ifdef CONV
  if ((.not. A%is_real) .and. (.not. B%is_real)) then
    call mm_ztrace(A,B,cmplx_alpha,label)
    alpha=real(cmplx_alpha,dp)
    return
  end if
#endif
  if (.not. A%is_real) call die('mm_dtrace: matrix A is complex')
  if (.not. B%is_real) call die('mm_dtrace: matrix B is complex')
  if ((A%dim1/=B%dim1) .or. &
      (A%dim2/=B%dim2)) call die('mm_dtrace: matrices A and B are not compatible')

  ! operation table
  if ((A%str_type .eq. 'den') .and. &
      (A%is_serial) .and. &
      (B%str_type .eq. 'den') .and. &
      (B%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else if ((A%str_type .eq. 'dbc') .and. &
           (.not. A%is_serial) .and. &
           (B%str_type .eq. 'dbc') .and. &
           (.not. B%is_serial)) then
      if (.not. present(label)) then
        ot=3
      else
        if (label .eq. 'lap') then
          ot=3
        end if
      end if
  else
    call die('mm_dtrace: invalid implementation')
  end if

  select case (ot)
    case (1)
      alpha=0.0_dp
      do i=1,A%dim1
        do j=1,A%dim2
          alpha=alpha+A%dval(i,j)*B%dval(i,j)
        end do
      end do
    case (2)
#ifdef LAP
      alpha=ddot(A%dim1*A%dim2,A%dval,1,B%dval,1)
#else
      call die('mm_dtrace: compile with LAPACK')
#endif
    case (3)
#if defined(MPI) && defined(LAP)
      if ((A%iaux2(1)/=B%iaux2(1)) .or. &
          (A%iaux2(2)/=B%iaux2(2))) call die('mm_dtrace: matrices A and B must have identical parallel distributions')
      alpha_loc=ddot(A%iaux2(1)*A%iaux2(2),A%dval,1,B%dval,1)
      call mpi_allreduce(alpha_loc,alpha,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
      if (info/=0) call die('mm_dtrace: error in mpi_allreduce')
#else
      call die('mm_dtrace: compile with MPI + LAPACK')
#endif
  end select

end subroutine mm_dtrace

subroutine mm_ztrace(A,B,alpha,label)
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  type(matrix), intent(in) :: A ! matrix A
  type(matrix), intent(in) :: B ! matrix B

  !**** OUTPUT **********************************!

  complex(dp), intent(out) :: alpha ! scalar alpha

  !**** INTERNAL ********************************!

  integer :: ot, i, j

#ifdef CONV
  real(dp) :: real_alpha
#endif

  !**** EXTERNAL ********************************!

#ifdef MPI
  integer :: info

  complex(dp) :: alpha_loc
#endif

  !**********************************************!

#ifdef CONV
  if (A%is_real .and. B%is_real) then
    call mm_dtrace(A,B,real_alpha,label)
    alpha=cmplx(real_alpha,0.0_dp,dp)
    return
  end if
#endif
  if (A%is_real) call die('mm_ztrace: matrix A is real')
  if (B%is_real) call die('mm_ztrace: matrix B is real')
  if ((A%dim1/=B%dim1) .or. &
      (A%dim2/=B%dim2)) call die('mm_ztrace: matrices A and B are not compatible')

  ! operation table
  if ((A%str_type .eq. 'den') .and. &
      (A%is_serial) .and. &
      (B%str_type .eq. 'den') .and. &
      (B%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((A%str_type .eq. 'dbc') .and. &
           (.not. A%is_serial) .and. &
           (B%str_type .eq. 'dbc') .and. &
           (.not. B%is_serial)) then
      if (.not. present(label)) then
        ot=3
      else
        if (label .eq. 'lap') then
          ot=3
        end if
      end if
  else
    call die('mm_ztrace: invalid implementation')
  end if

  select case (ot)
    case (1)
      alpha=cmplx_0
      do i=1,A%dim1
        do j=1,A%dim2
          alpha=alpha+conjg(A%zval(i,j))*B%zval(i,j)
        end do
      end do
    case (3)
#ifdef MPI
      if ((A%iaux2(1)/=B%iaux2(1)) .or. &
          (A%iaux2(2)/=B%iaux2(2))) call die('mm_ztrace: matrices A and B must have identical parallel distributions')
      alpha_loc=cmplx_0
      do i=1,A%iaux2(1)
        do j=1,A%iaux2(2)
          alpha_loc=alpha_loc+conjg(A%zval(i,j))*B%zval(i,j)
        end do
      end do
      call mpi_allreduce(alpha_loc,alpha,1,mpi_double_complex,mpi_sum,mpi_comm_world,info)
      if (info/=0) call die('mm_ztrace: error in mpi_allreduce')
#else
      call die('mm_ztrace: compile with MPI')
#endif
  end select

end subroutine mm_ztrace

!================================================!
! scale matrix                                   !
! C := beta*C                                    !
!================================================!
subroutine m_dscale(C,beta,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  real(dp), intent(in) :: beta ! scalar beta

  !**** INTOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot

#ifdef CONV
  complex(dp) :: cmplx_beta
#endif

  !**********************************************!

#ifdef CONV
  if (.not. C%is_real) then
    cmplx_beta=cmplx(beta,0.0_dp,dp)
    call m_zscale(C,cmplx_beta,label)
    return
  end if
#else
  if (.not. C%is_real) call die('m_dscale: matrix C is complex')
#endif

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else
    call die('m_dscale: invalid implementation')
  end if

  select case (ot)
    case (1)
      C%dval=beta*C%dval
  end select

end subroutine m_dscale

subroutine m_zscale(C,beta,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  complex(dp), intent(in) :: beta ! scalar beta

  !**** INTOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot

#ifdef CONV
  real(dp) :: real_beta
#endif

  !**********************************************!

#ifdef CONV
  if (C%is_real) then
    real_beta=real(beta,dp)
    call m_dscale(C,real_beta,label)
    return
  end if
#else
  if (C%is_real) call die('m_zscale: matrix C is complex')
#endif

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else
    call die('m_zscale: invalid implementation')
  end if

  select case (ot)
    case (1)
      C%zval=beta*C%zval
  end select

end subroutine m_zscale

!================================================!
! set matrix                                     !
! C_ij := alpha (off-diagonal elements) and      !
!         beta (diagonal elements)               !
!================================================!
subroutine m_dset(C,seC,alpha,beta,label)
  implicit none

  !**** INPUT ***********************************!

  character(1), intent(in) :: seC ! part of matrix to set: 'l/L' for lower, 'u/U' for upper, other for all
  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  real(dp), intent(in) :: alpha ! scalar alpha
  real(dp), intent(in) :: beta ! scalar beta

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot, i, j

#ifdef CONV
  complex(dp) :: cmplx_alpha, cmplx_beta
#endif

  !**********************************************!

#ifdef CONV
  if (.not. C%is_real) then
    cmplx_alpha=cmplx(alpha,0.0_dp,dp)
    cmplx_beta=cmplx(beta,0.0_dp,dp)
    call m_zset(C,seC,cmplx_alpha,cmplx_beta,label)
    return
  end if
#else
  if (.not. C%is_real) call die('m_dset: matrix C is complex')
#endif

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else
    call die('m_dset: invalid implementation')
  end if

  select case (ot)
    case (1)
      call m_set_sddenref(C,seC,alpha,beta)
    case (2)
#ifdef SLAP
      call pdlaset(seC,C%dim1,C%dim2,alpha,beta,C%dval,1,1,C%iaux1)
#else
      call die('m_dset: compile with ScaLAPACK')
#endif
  end select

end subroutine m_dset

subroutine m_zset(C,seC,alpha,beta,label)
  implicit none

  !**** INPUT ***********************************!

  character(1), intent(in) :: seC ! part of matrix to set: 'l/L' for lower, 'u/U' for upper, other for all
  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  complex(dp), intent(in) :: alpha ! scalar alpha
  complex(dp), intent(in) :: beta ! scalar beta

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot, i, j

#ifdef CONV
  real(dp) :: real_alpha, real_beta
#endif

  !**********************************************!

#ifdef CONV
  if (C%is_real) then
    real_alpha=real(alpha,dp)
    real_beta=real(beta,dp)
    call m_dset(C,seC,real_alpha,real_beta,label)
    return
  end if
#else
  if (C%is_real) call die('m_zset: matrix C is real')
#endif

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else
    call die('m_zset: invalid implementation')
  end if

  select case (ot)
    case (1)
      call m_set_szdenref(C,seC,alpha,beta)
    case (2)
#ifdef SLAP
      call pzlaset(seC,C%dim1,C%dim2,alpha,beta,C%zval,1,1,C%iaux1)
#else
      call die('m_zset: compile with ScaLAPACK')
#endif
  end select

end subroutine m_zset

!================================================!
! set matrix element                             !
! C_ij := alpha                                  !
!================================================!
subroutine m_dset_element(C,i,j,alpha,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  integer, intent(in) :: i ! row index of element
  integer, intent(in) :: j ! column index of element

  real(dp), intent(in) :: alpha ! scalar alpha

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot

#ifdef CONV
  complex(dp) :: cmplx_alpha
#endif

  !**********************************************!

#ifdef CONV
  if (.not. C%is_real) then
    cmplx_alpha=cmplx(alpha,0.0_dp,dp)
    call m_zset_element(C,i,j,cmplx_alpha,label)
    return
  end if
#else
  if (.not. C%is_real) call die('m_dset_element: matrix C is complex')
#endif
  if ((i<1) .or. &
      (i>C%dim1) .or. &
      (j<1) .or. &
      (j>C%dim2)) call die('m_dset_element: element out of range')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else
    call die('m_dset_element: invalid implementation')
  end if

  select case (ot)
    case (1)
      C%dval(i,j)=alpha
    case (2)
#ifdef SLAP
      call pdelset(C%dval,i,j,C%iaux1,alpha)
#else
      call die('m_dset_element: compile with ScaLAPACK')
#endif
  end select

end subroutine m_dset_element

subroutine m_zset_element(C,i,j,alpha,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  integer, intent(in) :: i ! row index of element
  integer, intent(in) :: j ! column index of element

  complex(dp), intent(in) :: alpha ! scalar alpha

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: C ! matrix C

  !**** INTERNAL ********************************!

  integer :: ot

#ifdef CONV
  real(dp) :: real_alpha
#endif

  !**********************************************!

#ifdef CONV
  if (C%is_real) then
    real_alpha=real(alpha,dp)
    call m_dset_element(C,i,j,real_alpha,label)
    return
  end if
#else
  if (C%is_real) call die('m_zset_element: matrix C is real')
#endif
  if ((i<1) .or. &
      (i>C%dim1) .or. &
      (j<1) .or. &
      (j>C%dim2)) call die('m_zset_element: element out of range')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else
    call die('m_zset_element: invalid implementation')
  end if

  select case (ot)
    case (1)
      C%zval(i,j)=alpha
    case (2)
#ifdef SLAP
      call pzelset(C%zval,i,j,C%iaux1,alpha)
#else
      call die('m_zset_element: compile with ScaLAPACK')
#endif
  end select

end subroutine m_zset_element

!================================================!
! get matrix element                             !
! alpha := C_ij                                  !
!================================================!
subroutine m_dget_element(C,i,j,alpha,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  integer, intent(in) :: i ! row index of element
  integer, intent(in) :: j ! column index of element

  type(matrix), intent(in) :: C ! matrix C

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: alpha ! scalar alpha

  !**** INTERNAL ********************************!

  integer :: ot

#ifdef CONV
  complex(dp) :: cmplx_alpha
#endif

  !**********************************************!

#ifdef CONV
  if (.not. C%is_real) then
    call m_zget_element(C,i,j,cmplx_alpha,label)
    alpha=real(cmplx_alpha,dp)
    return
  end if
#else
  if (.not. C%is_real) call die('m_dget_element: matrix C is complex')
#endif
  if ((i<1) .or. &
      (i>C%dim1) .or. &
      (j<1) .or. &
      (j>C%dim2)) call die('m_dget_element: element out of range')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else
    call die('m_dget_element: invalid implementation')
  end if

  select case (ot)
    case (1)
      alpha=C%dval(i,j)
    case (2)
#ifdef SLAP
      call pdelget('a',' ',alpha,C%dval,i,j,C%iaux1)
#else
      call die('m_dget_element: compile with ScaLAPACK')
#endif
  end select

end subroutine m_dget_element

subroutine m_zget_element(C,i,j,alpha,label)
  implicit none

  !**** INPUT ***********************************!

  character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

  integer, intent(in) :: i ! row index of element
  integer, intent(in) :: j ! column index of element

  type(matrix), intent(in) :: C ! matrix C

  !**** OUTPUT **********************************!

  complex(dp), intent(out) :: alpha ! scalar alpha

  !**** INTERNAL ********************************!

  integer :: ot

#ifdef CONV
  real(dp) :: real_alpha
#endif

  !**********************************************!

#ifdef CONV
  if (C%is_real) then
    call m_dget_element(C,i,j,real_alpha,label)
    alpha=cmplx(real_alpha,0.0_dp,dp)
    return
  end if
#else
  if (C%is_real) call die('m_zget_element: matrix C is real')
#endif
  if ((i<1) .or. &
      (i>C%dim1) .or. &
      (j<1) .or. &
      (j>C%dim2)) call die('m_zget_element: element out of range')

  ! operation table
  if ((C%str_type .eq. 'den') .and. &
      (C%is_serial)) then
      if (.not. present(label)) then
        ot=1
      else
        if (label .eq. 'ref') then
          ot=1
        else if (label .eq. 'lap') then
          ot=1
        end if
      end if
  else if ((C%str_type .eq. 'dbc') .and. &
           (.not. C%is_serial)) then
      if (.not. present(label)) then
        ot=2
      else
        if (label .eq. 'lap') then
          ot=2
        end if
      end if
  else
    call die('m_zget_element: invalid implementation')
  end if

  select case (ot)
    case (1)
      alpha=C%zval(i,j)
    case (2)
#ifdef SLAP
      call pzelget('a',' ',alpha,C%zval,i,j,C%iaux1)
#else
      call die('m_zget_element: compile with ScaLAPACK')
#endif
  end select

end subroutine m_zget_element

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

  m_name%zval => A

  m_name%is_initialized=.true.

end subroutine m_register_szden

!================================================!
! register matrix: dense block cyclic parallel   !
!================================================!
#ifdef MPI
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

  m_name%dval => A

  m_name%is_initialized=.true.

end subroutine m_register_pddbc
#endif

#ifdef MPI
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
  m_name%iaux1=desc
  m_name%dim1=desc(3)
  m_name%dim2=desc(4)
  allocate(m_name%iaux2(2))
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

  m_name%zval => A

  m_name%is_initialized=.true.

end subroutine m_register_pzdbc
#endif

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

subroutine process_lopM(opM,trM)
  implicit none

  !**** INPUT ***********************************!

  character(1), intent(in) :: opM

  !**** INOUT ***********************************!

  logical, intent(inout) :: trM

  !**********************************************!

  if ((opM .eq. 'T') .or. &
      (opM .eq. 't') .or. &
      (opM .eq. 'C') .or. &
      (opM .eq. 'c')) then
    trM=.true.
  else if ((opM .eq. 'N') .or. &
           (opM .eq. 'n')) then
    trM=.false.
  else
    call die('process_lopM: invalid opM')
  end if

end subroutine process_lopM

subroutine process_iopM(opM,tcM)
  implicit none

  !**** INPUT ***********************************!

  character(1), intent(in) :: opM

  !**** INOUT ***********************************!

  integer, intent(inout) :: tcM

  !**********************************************!

  if ((opM .eq. 'T') .or. &
      (opM .eq. 't')) then
    tcM=2
  else if ((opM .eq. 'C') .or. &
           (opM .eq. 'c')) then
    tcM=1
  else if ((opM .eq. 'N') .or. &
           (opM .eq. 'n')) then
    tcM=0
  else
    call die('process_iopM: invalid opM')
  end if

end subroutine process_iopM

subroutine process_seM(seM,luM)
  implicit none

  !**** INPUT ***********************************!

  character(1), intent(in) :: seM

  !**** INOUT ***********************************!

  integer, intent(inout) :: luM

  !**********************************************!

  if ((seM .eq. 'L') .or. &
      (seM .eq. 'l')) then
    luM=2
  else if ((seM .eq. 'U') .or. &
           (seM .eq. 'u')) then
    luM=1
  else
    luM=0
  end if

end subroutine process_seM

!================================================!
! implementation: ScaLAPACK                      !
!================================================!
#ifdef MPI
subroutine ms_scalapack_setup(mpi_size,nprow,order,bs_def,bs_list,icontxt)
  implicit none

  !**** INPUT ***********************************!

  character(1), intent(in) :: order ! ordering of processor grid: 'r/R' or other for row-major, 'c/C' for column-major

  integer, intent(in) :: mpi_size ! total number of MPI processes for the processor grid
  integer, intent(in) :: nprow ! number of rows in the processor grid
  integer, intent(in) :: bs_def ! default block size
  ! This is a list of exceptions to the default block size for specific matrix dimension sizes. The list has to be formatted as:
  !   * matrix dimension size 1
  !   * block size to use for matrix dimension size 1
  !   * matrix dimension size 2
  !   * block size to use for matrix dimension size 2
  !   * etc.
  integer, intent(in), optional :: bs_list(:)
  integer, intent(in), optional :: icontxt ! existing BLACS context handle in case ScaLAPACK is already initialized

  !**** INTERNAL ********************************!

  integer :: i

  !**********************************************!

  ms_mpi_size=mpi_size
  ms_lap_nprow=nprow
  ms_lap_npcol=mpi_size/nprow
  ms_lap_order=order
  ms_lap_bs_def=bs_def
  if (present(bs_list)) then
    if ((size(bs_list)/2)*2/=size(bs_list)) call die('ms_scalapack_setup: invalid bs_list')
    ms_lap_bs_num=size(bs_list)/2
    allocate(ms_lap_bs_list(2,ms_lap_bs_num))
    do i=1,ms_lap_bs_num
      ms_lap_bs_list(1,i)=bs_list((i-1)*2+1)
      ms_lap_bs_list(2,i)=bs_list((i-1)*2+2)
    end do
  end if

  if (present(icontxt)) then
    ms_lap_icontxt=icontxt
  else
    call blacs_get(-1,0,ms_lap_icontxt)
    call blacs_gridinit(ms_lap_icontxt,ms_lap_order,ms_lap_nprow,ms_lap_npcol)
  end if

end subroutine ms_scalapack_setup
#endif

#if defined(MPI) && defined(SLAP)
subroutine ms_scalapack_allocate(A)
  implicit none

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: A

  !**** INTERNAL ********************************!

  integer :: i, j, k, l, bs1, bs2, info

  !**** EXTERNAL ********************************!

  integer, external :: numroc

  !**********************************************!

  allocate(A%iaux1(9))
  allocate(A%iaux2(2))
  call blacs_gridinfo(ms_lap_icontxt,i,j,k,l)
  bs1=ms_lap_bs_def
  bs2=ms_lap_bs_def
  do i=1,ms_lap_bs_num
    if (ms_lap_bs_list(1,i)==A%dim1) then
      bs1=ms_lap_bs_list(2,i)
      exit
    end if
  end do
  do i=1,ms_lap_bs_num
    if (ms_lap_bs_list(1,i)==A%dim2) then
      bs2=ms_lap_bs_list(2,i)
      exit
    end if
  end do
  A%iaux2(1)=numroc(A%dim1,bs1,k,0,ms_lap_nprow)
  A%iaux2(2)=numroc(A%dim2,bs2,l,0,ms_lap_npcol)
  call descinit(A%iaux1,A%dim1,A%dim2,bs1,bs2,0,0,ms_lap_icontxt,A%iaux2(1),info)
  if (info/=0) call die('ms_scalapack_allocate: error in descinit')
  if (A%is_real) then
    allocate(A%dval(A%iaux2(1),A%iaux2(2)))
    A%dval=0.0_dp
  else
    allocate(A%zval(A%iaux2(1),A%iaux2(2)))
    A%zval=cmplx_0
  end if

end subroutine ms_scalapack_allocate
#endif

subroutine die(message)
  implicit none

  !**** INPUT ***********************************!

  character(*), intent(in), optional :: message

  !**********************************************!

  print('(a)'), 'FATAL ERROR in matrix_switch!'
  print('(a)'), message
  stop

end subroutine die

end module MatrixSwitch
