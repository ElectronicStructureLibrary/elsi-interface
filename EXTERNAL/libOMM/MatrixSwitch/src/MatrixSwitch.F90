module MatrixSwitch
  use MatrixSwitch_ops
  use MatrixSwitch_mm_multiply
  use MatrixSwitch_m_add
  use MatrixSwitch_m_set
  use MatrixSwitch_m_copy
  use MatrixSwitch_m_register
#ifdef PSP
  use pspBLAS
#endif

  implicit none

  private

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

  !************************************************!

  public :: matrix
  public :: m_allocate
  public :: m_deallocate
  public :: m_copy
  public :: m_convert
  public :: mm_multiply
  public :: m_add
  public :: m_trace
  public :: mm_trace
  public :: m_scale
  public :: m_set
  public :: m_set_element
  public :: m_get_element
  public :: m_register_sden
#ifdef MPI
  public :: m_register_pdbc
  public :: ms_scalapack_setup
  public :: ms_lap_icontxt
#endif
#ifdef PSP
  public :: m_register_psp_thre
  public :: m_register_psp_st
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
    if (m_name%is_serial) then
       if (m_name%str_type .eq. 'den') then
          m_name%is_sparse=.false.
          st=1
       else if (m_name%str_type .eq. 'coo') then
          m_name%is_sparse=.true.
          st=3
       else if (m_name%str_type .eq. 'csc') then
          m_name%is_sparse=.true.
          st=3
       else if (m_name%str_type .eq. 'csr') then
          m_name%is_sparse=.true.
          st=3
       else
          call die('m_allocate: invalid label')
       end if
    else
       if (m_name%str_type .eq. 'dbc') then
          m_name%is_sparse=.false.
          st=2
       else
          call die('m_allocate: invalid label')
       end if
    end if

    select case (st)
    case (1)
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          m_name%dval_is_allocated=.true.
          m_name%dval=0.0_dp
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          m_name%zval_is_allocated=.true.
          m_name%zval=cmplx_0
       end if
    case (2)
#if defined(MPI) && defined(SLAP)
       call ms_scalapack_allocate(m_name)
#else
       call die('m_allocate: compile with ScaLAPACK')
#endif
    case (3)
       allocate(m_name%iaux2(1))
       m_name%iaux2_is_allocated=.true.
       m_name%iaux2(1)=0
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

    if (associated(m_name%iaux1)) then
       if (m_name%iaux1_is_allocated) then
          deallocate(m_name%iaux1)
          m_name%iaux1_is_allocated=.false.
       else
          nullify(m_name%iaux1)
       end if
    end if
    if (associated(m_name%iaux2)) then
       if (m_name%iaux2_is_allocated) then
          deallocate(m_name%iaux2)
          m_name%iaux2_is_allocated=.false.
       else
          nullify(m_name%iaux2)
       end if
    end if
    if (associated(m_name%iaux3)) then
       if (m_name%iaux3_is_allocated) then
          deallocate(m_name%iaux3)
          m_name%iaux3_is_allocated=.false.
       else
          nullify(m_name%iaux3)
       end if
    end if
    if (associated(m_name%iaux4)) then
       if (m_name%iaux4_is_allocated) then
          deallocate(m_name%iaux4)
          m_name%iaux4_is_allocated=.false.
       else
          nullify(m_name%iaux4)
       end if
    end if
    if (associated(m_name%dval)) then
       if (m_name%dval_is_allocated) then
          deallocate(m_name%dval)
          m_name%dval_is_allocated=.false.
       else
          nullify(m_name%dval)
       end if
    end if
    if (associated(m_name%zval)) then
       if (m_name%zval_is_allocated) then
          deallocate(m_name%zval)
          m_name%zval_is_allocated=.false.
       else
          nullify(m_name%zval)
       end if
    end if

#ifdef PSP
    if (((m_name%str_type .eq. 'coo') .or. &
         (m_name%str_type .eq. 'csc')) .and. &
        (.not. m_name%is_serial)) call psp_deallocate_spm(m_name%spm)
#endif

    m_name%is_initialized=.false.

  end subroutine m_deallocate

  !================================================!
  ! copy matrix                                    !
  !================================================!
  subroutine m_copy(m_name,A,label,threshold,threshold_is_soft)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label ! storage format to use for m_name (see documentation)

    logical, intent(in), optional :: threshold_is_soft ! soft or hard thresholding

    real(dp), intent(in), optional :: threshold ! threshold for zeroing elements

    type(matrix), intent(inout) :: A ! matrix to copy from

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name ! matrix to copy onto

    !**** INTERNAL ********************************!

    character(1) :: c1, c2

    integer :: st, i, j, k

    real(dp) :: abs_threshold, soft_threshold

    !**********************************************!

    m_name%dim1=A%dim1
    m_name%dim2=A%dim2
    m_name%is_square=A%is_square

    if (present(label)) then
       read(label,'(a1,a1,a3)') c1, c2, m_name%str_type
       if (c1 .eq. 's') then
          m_name%is_serial=.true.
       else if (c1 .eq. 'p') then
          m_name%is_serial=.false.
#ifndef MPI
          call die('m_copy: compile with MPI')
#endif
       else
          call die('m_copy: invalid label')
       end if
       if (c2 .eq. 'd') then
          m_name%is_real=.true.
       else if (c2 .eq. 'z') then
          m_name%is_real=.false.
       else
          call die('m_copy: invalid label')
       end if
    else
       m_name%is_serial=A%is_serial
       m_name%is_real=A%is_real
       m_name%str_type=A%str_type
    end if

    if (m_name%is_real .neqv. A%is_real) call die('m_copy: invalid label')

    ! storage type
    if ((m_name%str_type .eq. 'den') .and. &
         (m_name%is_serial)) then
       m_name%is_sparse=.false.
       if ((A%str_type .eq. 'coo') .and. &
           (A%is_serial)) then
          st=9
       else if ((A%str_type .eq. 'csc') .and. &
                (A%is_serial)) then
          st=10
       else if ((A%str_type .eq. 'csr') .and. &
                (A%is_serial)) then
          st=11
       else if ((A%str_type .eq. 'den') .and. &
                (A%is_serial)) then
          st=1
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'dbc') .and. &
             (.not. m_name%is_serial) .and. &
             (A%str_type .eq. 'dbc') .and. &
             (.not. A%is_serial)) then
       m_name%is_sparse=.false.
       st=2
    else if ((m_name%str_type .eq. 'coo') .and. &
             (m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'coo') .and. &
           (A%is_serial)) then
          st=3
       else if ((A%str_type .eq. 'csc') .and. &
                (A%is_serial)) then
          st=12
       else if ((A%str_type .eq. 'csr') .and. &
                (A%is_serial)) then
          st=13
       else if ((A%str_type .eq. 'den') .and. &
                (A%is_serial)) then
          st=4
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'coo') .and. &
             (.not. m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'dbc') .and. &
           (.not. A%is_serial)) then
          st=16
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'csc') .and. &
             (m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'coo') .and. &
           (A%is_serial)) then
          st=14
       else if ((A%str_type .eq. 'csc') .and. &
           (A%is_serial)) then
          st=5
       else if ((A%str_type .eq. 'den') .and. &
                (A%is_serial)) then
          st=6
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'csc') .and. &
             (.not. m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'dbc') .and. &
           (.not. A%is_serial)) then
          st=17
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'csr') .and. &
             (m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'coo') .and. &
           (A%is_serial)) then
          st=15
       else if ((A%str_type .eq. 'csr') .and. &
           (A%is_serial)) then
          st=7
       else if ((A%str_type .eq. 'den') .and. &
                (A%is_serial)) then
          st=8
       else
          call die('m_copy: invalid label')
       end if
    else
       call die('m_copy: invalid label')
    end if

    if (present(threshold)) then
       abs_threshold=abs(threshold)
       if (present(threshold_is_soft)) then
          if (threshold_is_soft) then
             soft_threshold=abs_threshold
          else
             soft_threshold=0.0_dp
          end if
       else
          soft_threshold=0.0_dp
       end if
    else
       abs_threshold=0.0_dp
       soft_threshold=0.0_dp
    end if

    select case (st)
    case (1)
       if (present(threshold)) then
          call m_copy_sdensdenref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_sdensdenref(m_name,A)
       end if
    case (2)
       if (present(threshold)) then
          call m_copy_pdbcpdbcref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_pdbcpdbcref(m_name,A)
       end if
    case (3)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scooscooref(m_name,A)
       end if
    case (4)
       if (present(threshold)) then
          call m_copy_sdenscooref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_sdenscooref(m_name,A)
       end if
    case (5)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scscscscref(m_name,A)
       end if
    case (6)
       if (present(threshold)) then
          call m_copy_sdenscscref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_sdenscscref(m_name,A)
       end if
    case (7)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scsrscsrref(m_name,A)
       end if
    case (8)
       if (present(threshold)) then
          call m_copy_sdenscsrref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_sdenscsrref(m_name,A)
       end if
    case (9)
       if (present(threshold)) then
          call m_copy_scoosdenref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_scoosdenref(m_name,A)
       end if
    case (10)
       if (present(threshold)) then
          call m_copy_scscsdenref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_scscsdenref(m_name,A)
       end if
    case (11)
       if (present(threshold)) then
          call m_copy_scsrsdenref_thre(m_name,A,abs_threshold,soft_threshold)
       else
          call m_copy_scsrsdenref(m_name,A)
       end if
    case (12)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scscscooref(m_name,A)
       end if
    case (13)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scsrscooref(m_name,A)
       end if
    case (14)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scooscscref(m_name,A)
       end if
    case (15)
       if (present(threshold)) then
          call die('m_copy: thresholding not yet implemented')
       else
          call m_copy_scooscsrref(m_name,A)
       end if
    case (16)
       if (.not. present(threshold)) then
          call die('m_copy: threshold must be specified')
       else
          if (present(threshold_is_soft) .and. (threshold_is_soft)) then
             call die('m_copy: soft thresholding not yet implemented')
          else
#ifdef PSP
             if (A%is_real) then
                call m_register_pdsp_thre(m_name,A%dval,A%iaux1,'coo',threshold)
             else
                call m_register_pzsp_thre(m_name,A%zval,A%iaux1,'coo',threshold)
             end if
#else
             call die('mm_dmultiply: compile with pspBLAS')
#endif
          end if
       end if
    case (17)
       if (.not. present(threshold)) then
          call die('m_copy: threshold must be specified')
       else
          if (present(threshold_is_soft) .and. (threshold_is_soft)) then
             call die('m_copy: soft thresholding not yet implemented')
          else
#ifdef PSP
             if (A%is_real) then
                call m_register_pdsp_thre(m_name,A%dval,A%iaux1,'csc',threshold)
             else
                call m_register_pzsp_thre(m_name,A%zval,A%iaux1,'csc',threshold)
             end if
#else
             call die('mm_dmultiply: compile with pspBLAS')
#endif
          end if
       end if
    end select

    m_name%is_initialized=.true.

  end subroutine m_copy

  !================================================!
  ! wrapper for in-place matrix type conversion    !
  !================================================!
  subroutine m_convert(m_name,label,threshold,threshold_is_soft)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label ! storage format to use for m_name (see documentation)

    logical, intent(in), optional :: threshold_is_soft ! soft or hard thresholding

    real(dp), intent(in), optional :: threshold ! threshold for zeroing elements

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name ! matrix to convert

    !**** INTERNAL ********************************!

    type(matrix) :: temp_matrix

    !**********************************************!

    if (present(label)) then
       if (present(threshold)) then
          if (present(threshold_is_soft)) then
             call m_copy(temp_matrix,m_name,label,threshold,threshold_is_soft)
          else
             call m_copy(temp_matrix,m_name,label,threshold)
          end if
       else
          if (present(threshold_is_soft)) then
             call m_copy(temp_matrix,m_name,label,threshold_is_soft=threshold_is_soft)
          else
             call m_copy(temp_matrix,m_name,label)
          end if
       end if
    else
       if (present(threshold)) then
          if (present(threshold_is_soft)) then
             call m_copy(temp_matrix,m_name,threshold=threshold,threshold_is_soft=threshold_is_soft)
          else
             call m_copy(temp_matrix,m_name,threshold=threshold)
          end if
       else
          if (present(threshold_is_soft)) then
             call m_copy(temp_matrix,m_name,threshold_is_soft=threshold_is_soft)
          else
             call m_copy(temp_matrix,m_name)
          end if
       end if
    end if
    call m_deallocate(m_name)
    call m_copy(m_name,temp_matrix)
    call m_deallocate(temp_matrix)

  end subroutine m_convert

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
          else
             call die('mm_dmultiply: invalid implementation')
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
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'dbc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=4
       else
          if (label .eq. 'psp') then
             ot=4
          else if (label .eq. 't1D') then
             ot=6
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'csc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=5
       else
          if (label .eq. 'psp') then
             ot=5
          else if (label .eq. 't1D') then
             ot=7
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'dbc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'csc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=14
       else
          if (label .eq. 't1D') then
             ot=14
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=8
       else
          if (label .eq. 'ref') then
             ot=8
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'csc') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=9
       else
          if (label .eq. 'ref') then
             ot=9
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'csc') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=10
       else
          if (label .eq. 'ref') then
             ot=10
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csr') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=11
       else
          if (label .eq. 'ref') then
             ot=11
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'csr') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=12
       else
          if (label .eq. 'ref') then
             ot=12
          else
             call die('mm_dmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'csr') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=13
       else
          if (label .eq. 'ref') then
             ot=13
          else
             call die('mm_dmultiply: invalid implementation')
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
#if defined(MPI) && defined(SLAP)
       if (trA) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call pdgemm(opA,opB,C%dim1,C%dim2,i,alpha,A%dval,1,1,A%iaux1,B%dval,1,1,B%iaux1,beta,C%dval,1,1,C%iaux1)
#else
       call die('mm_dmultiply: compile with ScaLAPACK')
#endif
    case (4)
#ifdef PSP
       if (trA) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call psp_gespmm(C%dim1,C%dim2,i,A%spm,opA,B%dval,opB,C%dval,alpha,beta)
#else
       call die('mm_dmultiply: compile with pspBLAS')
#endif
    case (5)
#ifdef PSP
       if (trA) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call psp_gemspm(C%dim1,C%dim2,i,A%dval,opA,B%spm,opB,C%dval,alpha,beta)
#else
       call die('mm_dmultiply: compile with pspBLAS')
#endif
    case (6)
#ifdef PSP
       if (trB) then
          call die('mm_dmultiply: not implemented for transposed B')
       else
          call mm_multiply_pdcscpddbcpddbct1D(A,trA,B,C,alpha,beta)
       end if
#else
       call die('mm_dmultiply: compile with pspBLAS')
#endif
    case (7)
#ifdef PSP
       if (trA) then
          call die('mm_dmultiply: not implemented for transposed A')
       else
          call mm_multiply_pddbcpdcscpddbct1D(A,B,trB,C,alpha,beta)
       end if
#else
       call die('mm_dmultiply: compile with pspBLAS')
#endif
    case (8)
       call mm_multiply_sdcscsddensddenref(A,trA,B,trB,C,alpha,beta)
    case (9)
       call mm_multiply_sddensdcscsddenref(A,trA,B,trB,C,alpha,beta)
    case (10)
       call mm_multiply_sddensddensdcscref(A,trA,B,trB,C,alpha,beta)
    case (11)
       call mm_multiply_sdcsrsddensddenref(A,trA,B,trB,C,alpha,beta)
    case (12)
       call mm_multiply_sddensdcsrsddenref(A,trA,B,trB,C,alpha,beta)
    case (13)
       call mm_multiply_sddensddensdcsrref(A,trA,B,trB,C,alpha,beta)
    case (14)
#ifdef PSP
       if (trA .eqv. trB) then
          call die('mm_dmultiply: not implemented for combination of op(A) and op(B)')
       else
          call mm_multiply_pddbcpddbcpdcsct1D(A,trA,B,trB,C,alpha,beta)
       end if
#else
       call die('mm_dmultiply: compile with pspBLAS')
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
          else
             call die('mm_zmultiply: invalid implementation')
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
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'dbc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=4
       else
          if (label .eq. 'psp') then
             ot=4
          else if (label .eq. 't1D') then
             ot=6
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'csc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=5
       else
          if (label .eq. 'psp') then
             ot=5
          else if (label .eq. 't1D') then
             ot=7
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial) .and. &
         (B%str_type .eq. 'dbc') .and. &
         (.not. B%is_serial) .and. &
         (C%str_type .eq. 'csc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=14
       else
          if (label .eq. 't1D') then
             ot=14
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=8
       else
          if (label .eq. 'ref') then
             ot=8
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'csc') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=9
       else
          if (label .eq. 'ref') then
             ot=9
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'csc') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=10
       else
          if (label .eq. 'ref') then
             ot=10
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csr') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=11
       else
          if (label .eq. 'ref') then
             ot=11
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'csr') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=12
       else
          if (label .eq. 'ref') then
             ot=12
          else
             call die('mm_zmultiply: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'den') .and. &
         (A%is_serial) .and. &
         (B%str_type .eq. 'den') .and. &
         (B%is_serial) .and. &
         (C%str_type .eq. 'csr') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=13
       else
          if (label .eq. 'ref') then
             ot=13
          else
             call die('mm_zmultiply: invalid implementation')
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
#if defined(MPI) && defined(SLAP)
       if (tcA>0) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call pzgemm(opA,opB,C%dim1,C%dim2,i,alpha,A%zval,1,1,A%iaux1,B%zval,1,1,B%iaux1,beta,C%zval,1,1,C%iaux1)
#else
       call die('mm_zmultiply: compile with ScaLAPACK')
#endif
    case (4)
#ifdef PSP
       if (tcA>0) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call psp_gespmm(C%dim1,C%dim2,i,A%spm,opA,B%zval,opB,C%zval,alpha,beta)
#else
       call die('mm_zmultiply: compile with pspBLAS')
#endif
    case (5)
#ifdef PSP
       if (tcA>0) then
          i=A%dim1
       else
          i=A%dim2
       end if
       call psp_gemspm(C%dim1,C%dim2,i,A%zval,opA,B%spm,opB,C%zval,alpha,beta)
#else
       call die('mm_zmultiply: compile with pspBLAS')
#endif
    case (6)
#ifdef PSP
       if (tcB>0) then
          call die('mm_zmultiply: not implemented for transposed B')
       else
          call mm_multiply_pzcscpzdbcpzdbct1D(A,tcA,B,C,alpha,beta)
       end if
#else
       call die('mm_zmultiply: compile with pspBLAS')
#endif
    case (7)
#ifdef PSP
       if (tcA>0) then
          call die('mm_zmultiply: not implemented for transposed A')
       else
          call mm_multiply_pzdbcpzcscpzdbct1D(A,B,tcB,C,alpha,beta)
       end if
#else
       call die('mm_zmultiply: compile with pspBLAS')
#endif
    case (8)
       call mm_multiply_szcscszdenszdenref(A,tcA,B,tcB,C,alpha,beta)
    case (9)
       call mm_multiply_szdenszcscszdenref(A,tcA,B,tcB,C,alpha,beta)
    case (10)
       call mm_multiply_szdenszdenszcscref(A,tcA,B,tcB,C,alpha,beta)
    case (11)
       call mm_multiply_szcsrszdenszdenref(A,tcA,B,tcB,C,alpha,beta)
    case (12)
       call mm_multiply_szdenszcsrszdenref(A,tcA,B,tcB,C,alpha,beta)
    case (13)
       call mm_multiply_szdenszdenszcsrref(A,tcA,B,tcB,C,alpha,beta)
    case (14)
#ifdef PSP
       if ((tcA==tcB) .or. (tcA+tcB>2)) then
          call die('mm_zmultiply: not implemented for combination of op(A) and op(B)')
       else
          call mm_multiply_pzdbcpzdbcpzcsct1D(A,tcA,B,tcB,C,alpha,beta)
       end if
#else
       call die('mm_zmultiply: compile with pspBLAS')
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
          else
             call die('m_dadd: invalid implementation')
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
          else if (label .eq. 'psp') then
             ot=2
          else if (label .eq. 't1D') then
             ot=2
          else
             call die('m_dadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (.not. A%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('m_dadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (A%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=4
       else
          if (label .eq. 'ref') then
             ot=4
          else
             call die('m_dadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csr') .and. &
         (A%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=5
       else
          if (label .eq. 'ref') then
             ot=5
          else
             call die('m_dadd: invalid implementation')
          end if
       end if
    else
       call die('m_dadd: invalid implementation')
    end if

    select case (ot)
    case (1)
       call m_add_sddenref(A,trA,C,alpha,beta)
    case (2)
#if defined(MPI) && defined(SLAP)
       call pdgeadd(opA,C%dim1,C%dim2,alpha,A%dval,1,1,A%iaux1,beta,C%dval,1,1,C%iaux1)
#else
       call die('m_dadd: compile with ScaLAPACK')
#endif
    case (3)
#ifdef PSP
       if (trA) call die('m_dadd: implementation only valid for opA=''n''')
       if ((A%spm%loc_dim1/=C%iaux2(1)) .or. &
           (A%spm%loc_dim2/=C%iaux2(2))) call die('m_dadd: matrices A and C must have identical parallel distributions')
       call m_add_pdcscpddbcref(A,C,alpha,beta)
#else
       call die('m_dadd: compile with pspBLAS')
#endif
    case (4)
       call m_add_sdcscsddenref(A,trA,C,alpha,beta)
    case (5)
       call m_add_sdcsrsddenref(A,trA,C,alpha,beta)
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
          else
             call die('m_zadd: invalid implementation')
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
          else if (label .eq. 'psp') then
             ot=2
          else if (label .eq. 't1D') then
             ot=2
          else
             call die('m_zadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (.not. A%is_serial) .and. &
         (C%str_type .eq. 'dbc') .and. &
         (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('m_dadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csc') .and. &
         (A%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=4
       else
          if (label .eq. 'ref') then
             ot=4
          else
             call die('m_dadd: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'csr') .and. &
         (A%is_serial) .and. &
         (C%str_type .eq. 'den') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=5
       else
          if (label .eq. 'ref') then
             ot=5
          else
             call die('m_dadd: invalid implementation')
          end if
       end if
    else
       call die('m_zadd: invalid implementation')
    end if

    select case (ot)
    case (1)
       call m_add_szdenref(A,tcA,C,alpha,beta)
    case (2)
#if defined(MPI) && defined(SLAP)
       call pzgeadd(opA,C%dim1,C%dim2,alpha,A%zval,1,1,A%iaux1,beta,C%zval,1,1,C%iaux1)
#else
       call die('m_zadd: compile with ScaLAPACK')
#endif
    case (3)
#ifdef PSP
       if (tcA>0) call die('m_zadd: implementation only valid for opA=''n''')
       if ((A%spm%loc_dim1/=C%iaux2(1)) .or. &
           (A%spm%loc_dim2/=C%iaux2(2))) call die('m_zadd: matrices A and C must have identical parallel distributions')
       call m_add_pzcscpzdbcref(A,C,alpha,beta)
#else
       call die('m_dadd: compile with pspBLAS')
#endif
    case (4)
       call m_add_szcscszdenref(A,tcA,C,alpha,beta)
    case (5)
       call m_add_szcsrszdenref(A,tcA,C,alpha,beta)
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
          else
             call die('m_dtrace: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial)) then
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
             call die('m_dtrace: invalid implementation')
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
#if defined(MPI) && defined(SLAP)
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
    complex(dp), external :: pzlatra
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
          else
             call die('m_ztrace: invalid implementation')
          end if
       end if
    else if ((A%str_type .eq. 'dbc') .and. &
         (.not. A%is_serial)) then
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
             call die('m_ztrace: invalid implementation')
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
#if defined(MPI) && defined(SLAP)
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
          else
             call die('mm_dtrace: invalid implementation')
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
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('mm_dtrace: invalid implementation')
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
       call mpi_allreduce(alpha_loc,alpha,1,mpi_double_precision,mpi_sum,ms_mpi_comm,info)
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
          else
             call die('mm_ztrace: invalid implementation')
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
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('mm_ztrace: invalid implementation')
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
       call mpi_allreduce(alpha_loc,alpha,1,mpi_double_complex,mpi_sum,ms_mpi_comm,info)
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
          else
             call die('m_dscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'dbc') .and. &
             (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_dscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'coo') .and. &
             (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_dscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'csc') .and. &
             (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_dscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'csr') .and. &
             (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_dscale: invalid implementation')
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
          else
             call die('m_zscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'dbc') .and. &
             (.not. C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_zscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'coo') .and. &
             (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_dscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'csc') .and. &
             (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_dscale: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'csr') .and. &
             (C%is_serial)) then
       if (.not. present(label)) then
          ot=1
       else
          if (label .eq. 'lap') then
             ot=1
          else if (label .eq. 'psp') then
             ot=1
          else if (label .eq. 't1D') then
             ot=1
          else
             call die('m_dscale: invalid implementation')
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
          else
             call die('m_dset: invalid implementation')
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
             call die('m_dset: invalid implementation')
          end if
       end if
    else
       call die('m_dset: invalid implementation')
    end if

    select case (ot)
    case (1)
       call m_set_sddenref(C,seC,alpha,beta)
    case (2)
#if defined(MPI) && defined(SLAP)
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
          else
             call die('m_zset: invalid implementation')
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
             call die('m_zset: invalid implementation')
          end if
       end if
    else
       call die('m_zset: invalid implementation')
    end if

    select case (ot)
    case (1)
       call m_set_szdenref(C,seC,alpha,beta)
    case (2)
#if defined(MPI) && defined(SLAP)
       call pzlaset(seC,C%dim1,C%dim2,alpha,beta,C%zval,1,1,C%iaux1)
#else
       call die('m_zset: compile with ScaLAPACK')
#endif
    end select

  end subroutine m_zset

  !================================================!
  ! set matrix element                             !
  ! C_ij := alpha + beta*C_ij                      !
  !================================================!
  subroutine m_dset_element(C,i,j,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

    integer, intent(in) :: i ! row index of element
    integer, intent(in) :: j ! column index of element

    real(dp), intent(in) :: alpha ! scalar alpha
    real(dp), intent(in) :: beta ! scalar beta

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C ! matrix C

    !**** INTERNAL ********************************!

    logical :: el_present

    integer :: ot, k, buffer
    integer, allocatable :: iaux3_temp(:), iaux4_temp(:)

    real(dp) :: el
    real(dp), allocatable :: dval_temp(:,:)

#ifdef CONV
    complex(dp) :: cmplx_alpha, cmplx_beta
#endif

    !**********************************************!

#ifdef CONV
    if (.not. C%is_real) then
       cmplx_alpha=cmplx(alpha,0.0_dp,dp)
       cmplx_beta=cmplx(beta,0.0_dp,dp)
       call m_zset_element(C,i,j,cmplx_alpha,cmplx_beta,label)
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
          else
            call die('m_dset_element: invalid implementation')
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
            call die('m_dset_element: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'coo') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'ref') then
             ot=3
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
            call die('m_dset_element: invalid implementation')
          end if
       end if
    else
       call die('m_dset_element: invalid implementation')
    end if

    select case (ot)
    case (1)
       C%dval(i,j)=alpha+beta*C%dval(i,j)
    case (2)
#if defined(MPI) && defined(SLAP)
       if (beta==0.0_dp) then
          call pdelset(C%dval,i,j,C%iaux1,alpha)
       else
          call pdelget('a',' ',el,C%dval,i,j,C%iaux1)
          call pdelset(C%dval,i,j,C%iaux1,alpha+beta*el)
       end if
#else
       call die('m_dset_element: compile with ScaLAPACK')
#endif
    case (3)
       if (C%iaux2(1)==0) then
          C%iaux2(1)=1
          buffer=min(C%dim1,C%dim2)
          allocate(C%iaux3(buffer))
          C%iaux3_is_allocated=.true.
          allocate(C%iaux4(buffer))
          C%iaux4_is_allocated=.true.
          allocate(C%dval(buffer,1))
          C%dval_is_allocated=.true.
          C%iaux3(1)=i
          C%iaux4(1)=j
          C%dval(1,1)=alpha
       else
          el_present=.false.
          do k=1,C%iaux2(1)
             if ((C%iaux3(k)==i) .and. &
                 (C%iaux4(k)==j)) then
                C%dval(k,1)=alpha+beta*C%dval(k,1)
                el_present=.true.
                exit
             end if
          end do
          if (.not. el_present) then
             if (C%iaux2(1)==size(C%iaux3)) then
                allocate(iaux3_temp(C%iaux2(1)))
                allocate(iaux4_temp(C%iaux2(1)))
                allocate(dval_temp(C%iaux2(1),1))
                iaux3_temp=C%iaux3
                iaux4_temp=C%iaux4
                dval_temp=C%dval
                deallocate(C%dval)
                C%dval_is_allocated=.false.
                deallocate(C%iaux4)
                C%iaux4_is_allocated=.false.
                deallocate(C%iaux3)
                C%iaux3_is_allocated=.false.
                buffer=C%iaux2(1)+min(C%dim1,C%dim2)
                allocate(C%iaux3(buffer))
                C%iaux3_is_allocated=.true.
                allocate(C%iaux4(buffer))
                C%iaux4_is_allocated=.true.
                allocate(C%dval(buffer,1))
                C%dval_is_allocated=.true.
                C%iaux3(1:C%iaux2(1))=iaux3_temp(1:C%iaux2(1))
                C%iaux4(1:C%iaux2(1))=iaux4_temp(1:C%iaux2(1))
                C%dval(1:C%iaux2(1),1)=dval_temp(1:C%iaux2(1),1)
                deallocate(dval_temp)
                deallocate(iaux4_temp)
                deallocate(iaux3_temp)
             end if
             C%iaux2(1)=C%iaux2(1)+1
             C%iaux3(C%iaux2(1))=i
             C%iaux4(C%iaux2(1))=j
             C%dval(C%iaux2(1),1)=alpha
          end if
       end if
    end select

  end subroutine m_dset_element

  subroutine m_zset_element(C,i,j,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

    integer, intent(in) :: i ! row index of element
    integer, intent(in) :: j ! column index of element

    complex(dp), intent(in) :: alpha ! scalar alpha
    complex(dp), intent(in) :: beta ! scalar beta

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: C ! matrix C

    !**** INTERNAL ********************************!

    logical :: el_present

    integer :: ot, k, buffer
    integer, allocatable :: iaux3_temp(:), iaux4_temp(:)

#ifdef CONV
    real(dp) :: real_alpha, real_beta
#endif

    complex(dp) :: el
    complex(dp), allocatable :: zval_temp(:,:)

    !**********************************************!

#ifdef CONV
    if (C%is_real) then
       real_alpha=real(alpha,dp)
       real_beta=real(beta,dp)
       call m_dset_element(C,i,j,real_alpha,real_beta,label)
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
          else
             call die('m_zset_element: invalid implementation')
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
             call die('m_zset_element: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'coo') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'ref') then
             ot=3
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('m_zset_element: invalid implementation')
          end if
       end if
    else
       call die('m_zset_element: invalid implementation')
    end if

    select case (ot)
    case (1)
       C%zval(i,j)=alpha+beta*C%zval(i,j)
    case (2)
#if defined(MPI) && defined(SLAP)
       if (beta==cmplx_0) then
          call pzelset(C%zval,i,j,C%iaux1,alpha)
       else
          call pzelget('a',' ',el,C%zval,i,j,C%iaux1)
          call pzelset(C%zval,i,j,C%iaux1,alpha+beta*el)
       end if
#else
       call die('m_zset_element: compile with ScaLAPACK')
#endif
    case (3)
       if (C%iaux2(1)==0) then
          C%iaux2(1)=1
          buffer=min(C%dim1,C%dim2)
          allocate(C%iaux3(buffer))
          C%iaux3_is_allocated=.true.
          allocate(C%iaux4(buffer))
          C%iaux4_is_allocated=.true.
          allocate(C%zval(buffer,1))
          C%zval_is_allocated=.true.
          C%iaux3(1)=i
          C%iaux4(1)=j
          C%zval(1,1)=alpha
       else
          el_present=.false.
          do k=1,C%iaux2(1)
             if ((C%iaux3(k)==i) .and. &
                 (C%iaux4(k)==j)) then
                C%zval(k,1)=alpha+beta*C%zval(k,1)
                el_present=.true.
                exit
             end if
          end do
          if (.not. el_present) then
             if (C%iaux2(1)==size(C%iaux3)) then
                allocate(iaux3_temp(C%iaux2(1)))
                allocate(iaux4_temp(C%iaux2(1)))
                allocate(zval_temp(C%iaux2(1),1))
                iaux3_temp=C%iaux3
                iaux4_temp=C%iaux4
                zval_temp=C%zval
                deallocate(C%zval)
                C%zval_is_allocated=.false.
                deallocate(C%iaux4)
                C%iaux4_is_allocated=.false.
                deallocate(C%iaux3)
                C%iaux3_is_allocated=.false.
                buffer=C%iaux2(1)+min(C%dim1,C%dim2)
                allocate(C%iaux3(buffer))
                C%iaux3_is_allocated=.true.
                allocate(C%iaux4(buffer))
                C%iaux4_is_allocated=.true.
                allocate(C%zval(buffer,1))
                C%zval_is_allocated=.true.
                C%iaux3(1:C%iaux2(1))=iaux3_temp(1:C%iaux2(1))
                C%iaux4(1:C%iaux2(1))=iaux4_temp(1:C%iaux2(1))
                C%zval(1:C%iaux2(1),1)=zval_temp(1:C%iaux2(1),1)
                deallocate(zval_temp)
                deallocate(iaux4_temp)
                deallocate(iaux3_temp)
             end if
             C%iaux2(1)=C%iaux2(1)+1
             C%iaux3(C%iaux2(1))=i
             C%iaux4(C%iaux2(1))=j
             C%zval(C%iaux2(1),1)=alpha
          end if
       end if
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

    integer :: ot, k

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
          else
             call die('m_dget_element: invalid implementation')
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
             call die('m_dget_element: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'coo') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'ref') then
             ot=3
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('m_dget_element: invalid implementation')
          end if
       end if
    else
       call die('m_dget_element: invalid implementation')
    end if

    select case (ot)
    case (1)
       alpha=C%dval(i,j)
    case (2)
#if defined(MPI) && defined(SLAP)
       call pdelget('a',' ',alpha,C%dval,i,j,C%iaux1)
#else
       call die('m_dget_element: compile with ScaLAPACK')
#endif
    case (3)
       alpha=0.0_dp
       do k=1,C%iaux2(1)
          if ((C%iaux3(k)==i) .and. &
              (C%iaux4(k)==j)) then
             alpha=C%dval(k,1)
             exit
          end if
       end do
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

    integer :: ot, k

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
          else
             call die('m_zget_element: invalid implementation')
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
             call die('m_zget_element: invalid implementation')
          end if
       end if
    else if ((C%str_type .eq. 'coo') .and. &
         (C%is_serial)) then
       if (.not. present(label)) then
          ot=3
       else
          if (label .eq. 'ref') then
             ot=3
          else if (label .eq. 'psp') then
             ot=3
          else if (label .eq. 't1D') then
             ot=3
          else
             call die('m_zget_element: invalid implementation')
          end if
       end if
    else
       call die('m_zget_element: invalid implementation')
    end if

    select case (ot)
    case (1)
       alpha=C%zval(i,j)
    case (2)
#if defined(MPI) && defined(SLAP)
       call pzelget('a',' ',alpha,C%zval,i,j,C%iaux1)
#else
       call die('m_zget_element: compile with ScaLAPACK')
#endif
    case (3)
       alpha=cmplx_0
       do k=1,C%iaux2(1)
          if ((C%iaux3(k)==i) .and. &
              (C%iaux4(k)==j)) then
             alpha=C%zval(k,1)
             exit
          end if
       end do
    end select

  end subroutine m_zget_element

  !================================================!
  ! implementation: ScaLAPACK                      !
  !================================================!
#if defined(MPI) && defined(SLAP)
  subroutine ms_scalapack_setup(mpi_comm,nprow,order,bs_def,bs_list,icontxt)
    implicit none
    include 'mpif.h'

    !**** INPUT ***********************************!

    character(1), intent(in) :: order ! ordering of processor grid: 'r/R' or other for row-major, 'c/C' for column-major

    integer, intent(in) :: mpi_comm ! MPI communicator
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

    integer :: i, mpi_err

    !**********************************************!

    ms_mpi_comm=mpi_comm
    call mpi_comm_size(ms_mpi_comm,ms_mpi_size,mpi_err)
    call mpi_comm_rank(ms_mpi_comm,ms_mpi_rank,mpi_err)
    ms_lap_nprow=nprow
    ms_lap_npcol=ms_mpi_size/nprow
    if (ms_lap_nprow*ms_lap_npcol/=ms_mpi_size) call die('ms_scalapack_setup: invalid nprow')
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
       ms_lap_icontxt=ms_mpi_comm
       call blacs_gridinit(ms_lap_icontxt,ms_lap_order,ms_lap_nprow,ms_lap_npcol)
    end if

#ifdef PSP
    ! initialized grid information in pspBLAS
    call psp_gridinit_2D(ms_mpi_comm,ms_mpi_size,ms_lap_nprow,ms_lap_order,ms_lap_bs_def,ms_lap_bs_def,ms_lap_icontxt)
#endif

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
    A%iaux1_is_allocated=.true.
    allocate(A%iaux2(2))
    A%iaux2_is_allocated=.true.
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
       A%dval_is_allocated=.true.
       A%dval=0.0_dp
    else
       allocate(A%zval(A%iaux2(1),A%iaux2(2)))
       A%zval_is_allocated=.true.
       A%zval=cmplx_0
    end if

  end subroutine ms_scalapack_allocate
#endif

end module MatrixSwitch
