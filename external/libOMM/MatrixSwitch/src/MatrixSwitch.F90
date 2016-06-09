module MatrixSwitch
#ifdef PSP
  use pspBLAS
#endif

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
  integer, save :: ms_mpi_rank
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
     logical :: is_sparse ! is the matrix sparse?

     integer :: dim1 ! (global) row dimension size of the matrix
     integer :: dim2 ! (global) column dimension size of the matrix
     integer, pointer :: iaux1(:) => null() ! auxiliary information for certain storage formats
     integer, pointer :: iaux2(:) => null() ! auxiliary information for certain storage formats
     integer, pointer :: iaux3(:) => null() ! auxiliary information for certain storage formats
     integer, pointer :: iaux4(:) => null() ! auxiliary information for certain storage formats

     real(dp), pointer :: dval(:,:) => null() ! matrix elements for a real matrix

     complex(dp), pointer :: zval(:,:) => null() ! matrix elements for a complex matrix

#ifdef PSP
     type(psp_matrix_spm) :: spm ! a sparse matrix in pspBLAS
#endif
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

#ifdef PSP
  interface m_register_psp_thre
     module procedure m_register_pdsp_thre
     module procedure m_register_pzsp_thre
  end interface m_register_psp_thre

  interface m_register_psp_st
     module procedure m_register_pdsp_st
     module procedure m_register_pzsp_st
  end interface m_register_psp_st
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
  public :: m_copy
  public :: m_convert
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
    case (3)
       allocate(m_name%iaux2(1))
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

    if (associated(m_name%iaux1)) nullify(m_name%iaux1)
    if (associated(m_name%iaux2)) nullify(m_name%iaux2)
    if (associated(m_name%iaux3)) nullify(m_name%iaux3)
    if (associated(m_name%iaux4)) nullify(m_name%iaux4)
    if (associated(m_name%dval)) nullify(m_name%dval)
    if (associated(m_name%zval)) nullify(m_name%zval)

#ifdef PSP
    if ((m_name%str_type .eq. 'coo') .or. &
        (m_name%str_type .eq. 'csc')) call psp_deallocate_spm(m_name%spm)
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
         (m_name%is_serial) .and. &
         (A%str_type .eq. 'den') .and. &
         (A%is_serial)) then
       m_name%is_sparse=.false.
       st=1
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
       else if ((A%str_type .eq. 'den') .and. &
                (A%is_serial)) then
          st=4
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'csc') .and. &
             (m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'csc') .and. &
           (A%is_serial)) then
          st=5
       else if ((A%str_type .eq. 'den') .and. &
                (A%is_serial)) then
          st=6
       else
          call die('m_copy: invalid label')
       end if
    else if ((m_name%str_type .eq. 'csr') .and. &
             (m_name%is_serial)) then
       m_name%is_sparse=.true.
       if ((A%str_type .eq. 'csr') .and. &
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
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%dim1,m_name%dim2))
          if (present(threshold)) then
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
             m_name%dval=A%dval
          end if
       else
          allocate(m_name%zval(m_name%dim1,m_name%dim2))
          if (present(threshold)) then
             do i=1,m_name%dim1
                do j=1,m_name%dim2
                   if (abs(A%zval(i,j))>abs_threshold) then
                      m_name%zval(i,j)=A%zval(i,j)-soft_threshold*A%zval(i,j)/abs(A%zval(i,j))
                   else
                      m_name%zval(i,j)=cmplx_0
                   end if
                end do
             end do
          else
             m_name%zval=A%zval
          end if
       end if
    case (2)
       allocate(m_name%iaux1(9))
       allocate(m_name%iaux2(2))
       m_name%iaux1=A%iaux1
       m_name%iaux2=A%iaux2
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),m_name%iaux2(2)))
          if (present(threshold)) then
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
             m_name%dval=A%dval
          end if
       else
          allocate(m_name%zval(m_name%iaux2(1),m_name%iaux2(2)))
          if (present(threshold)) then
             do i=1,m_name%iaux2(1)
                do j=1,m_name%iaux2(2)
                   if (abs(A%zval(i,j))>abs_threshold) then
                      m_name%zval(i,j)=A%zval(i,j)-soft_threshold*A%zval(i,j)/abs(A%zval(i,j))
                   else
                      m_name%zval(i,j)=cmplx_0
                   end if
                end do
             end do
          else
             m_name%zval=A%zval
          end if
       end if
    case (3)
       allocate(m_name%iaux2(1))
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%iaux2(1)))
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux3=A%iaux3
       m_name%iaux4=A%iaux4
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval=A%zval
       end if
    case (4)
       allocate(m_name%iaux2(1))
       m_name%iaux2(1)=0
       if (m_name%is_real) then
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%dval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%iaux2(1)))
             allocate(m_name%iaux4(m_name%iaux2(1)))
             allocate(m_name%dval(m_name%iaux2(1),1))
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
             allocate(m_name%iaux4(m_name%iaux2(1)))
             allocate(m_name%zval(m_name%iaux2(1),1))
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
    case (5)
       allocate(m_name%iaux2(1))
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%iaux2(1)))
       allocate(m_name%iaux4(m_name%dim2+1))
       m_name%iaux3=A%iaux3
       m_name%iaux4=A%iaux4
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval=A%zval
       end if
    case (6)
       allocate(m_name%iaux2(1))
       m_name%iaux2(1)=0
       if (m_name%is_real) then
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%dval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%iaux2(1)))
             allocate(m_name%iaux4(m_name%dim2+1))
             allocate(m_name%dval(m_name%iaux2(1),1))
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
             allocate(m_name%iaux4(m_name%dim2+1))
             allocate(m_name%zval(m_name%iaux2(1),1))
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
    case (7)
       allocate(m_name%iaux2(1))
       m_name%iaux2(1)=A%iaux2(1)
       allocate(m_name%iaux3(m_name%dim1+1))
       allocate(m_name%iaux4(m_name%iaux2(1)))
       m_name%iaux3=A%iaux3
       m_name%iaux4=A%iaux4
       if (m_name%is_real) then
          allocate(m_name%dval(m_name%iaux2(1),1))
          m_name%dval=A%dval
       else
          allocate(m_name%zval(m_name%iaux2(1),1))
          m_name%zval=A%zval
       end if
    case (8)
       allocate(m_name%iaux2(1))
       m_name%iaux2(1)=0
       if (m_name%is_real) then
          do i=1,m_name%dim2
             do j=1,m_name%dim1
                if (abs(A%dval(j,i))>abs_threshold) m_name%iaux2(1)=m_name%iaux2(1)+1
             end do
          end do
          if (m_name%iaux2(1)>0) then
             allocate(m_name%iaux3(m_name%dim1+1))
             allocate(m_name%iaux4(m_name%iaux2(1)))
             allocate(m_name%dval(m_name%iaux2(1),1))
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
             allocate(m_name%iaux4(m_name%iaux2(1)))
             allocate(m_name%zval(m_name%iaux2(1),1))
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
       call mm_multiply_pdcscpddbcref(A,trA,B,trB,C,alpha,beta)
#else
       call die('mm_dmultiply: compile with pspBLAS')
#endif
    case (7)
#ifdef PSP
       call mm_multiply_pddbcpdcscref(A,trA,B,trB,C,alpha,beta)
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
       call mm_multiply_pzcscpzdbcref(A,tcA,B,tcB,C,alpha,beta)
#else
       call die('mm_zmultiply: compile with pspBLAS')
#endif
    case (7)
#ifdef PSP
       call mm_multiply_pzdbcpzcscref(A,tcA,B,tcB,C,alpha,beta)
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
          else if (label .eq. 'psp') then
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
#if defined(MPI) && defined(SLAP)
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
          else if (label .eq. 'psp') then
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
          else if (label .eq. 'psp') then
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
          else if (label .eq. 'psp') then
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
          else if (label .eq. 'psp') then
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
          else if (label .eq. 'psp') then
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

    logical :: el_present

    integer :: ot, k
    integer, allocatable :: iaux3_temp(:), iaux4_temp(:)

    real(dp), allocatable :: dval_temp(:,:)

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
          else if (label .eq. 'psp') then
             ot=2
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
          end if
       end if
    else
       call die('m_dset_element: invalid implementation')
    end if

    select case (ot)
    case (1)
       C%dval(i,j)=alpha
    case (2)
#if defined(MPI) && defined(SLAP)
       call pdelset(C%dval,i,j,C%iaux1,alpha)
#else
       call die('m_dset_element: compile with ScaLAPACK')
#endif
    case (3)
       el_present=.false.
       do k=1,C%iaux2(1)
          if ((C%iaux3(k)==i) .and. &
               (C%iaux4(k)==j)) then
             C%dval(k,1)=alpha
             el_present=.true.
          end if
       end do
       if (.not. el_present) then
          allocate(iaux3_temp(C%iaux2(1)))
          allocate(iaux4_temp(C%iaux2(1)))
          allocate(dval_temp(C%iaux2(1),1))
          iaux3_temp=C%iaux3
          iaux4_temp=C%iaux4
          dval_temp=C%dval
          deallocate(C%dval)
          deallocate(C%iaux4)
          deallocate(C%iaux3)
          C%iaux2(1)=C%iaux2(1)+1
          allocate(C%iaux3(C%iaux2(1)))
          allocate(C%iaux4(C%iaux2(1)))
          allocate(C%dval(C%iaux2(1),1))
          C%iaux3(1:C%iaux2(1)-1)=iaux3_temp(1:C%iaux2(1)-1)
          C%iaux4(1:C%iaux2(1)-1)=iaux4_temp(1:C%iaux2(1)-1)
          C%dval(1:C%iaux2(1)-1,1)=dval_temp(1:C%iaux2(1)-1,1)
          deallocate(dval_temp)
          deallocate(iaux4_temp)
          deallocate(iaux3_temp)
          C%iaux3(C%iaux2(1))=i
          C%iaux4(C%iaux2(1))=j
          C%dval(C%iaux2(1),1)=alpha
       end if
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

    logical :: el_present

    integer :: ot, k
    integer, allocatable :: iaux3_temp(:), iaux4_temp(:)

#ifdef CONV
    real(dp) :: real_alpha
#endif

    complex(dp), allocatable :: zval_temp(:,:)

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
          else if (label .eq. 'psp') then
             ot=2
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
          end if
       end if
    else
       call die('m_zset_element: invalid implementation')
    end if

    select case (ot)
    case (1)
       C%zval(i,j)=alpha
    case (2)
#if defined(MPI) && defined(SLAP)
       call pzelset(C%zval,i,j,C%iaux1,alpha)
#else
       call die('m_zset_element: compile with ScaLAPACK')
#endif
    case (3)
       el_present=.false.
       do k=1,C%iaux2(1)
          if ((C%iaux3(k)==i) .and. &
               (C%iaux4(k)==j)) then
             C%zval(k,1)=alpha
             el_present=.true.
          end if
       end do
       if (.not. el_present) then
          allocate(iaux3_temp(C%iaux2(1)))
          allocate(iaux4_temp(C%iaux2(1)))
          allocate(zval_temp(C%iaux2(1),1))
          iaux3_temp=C%iaux3
          iaux4_temp=C%iaux4
          zval_temp=C%zval
          deallocate(C%zval)
          deallocate(C%iaux4)
          deallocate(C%iaux3)
          C%iaux2(1)=C%iaux2(1)+1
          allocate(C%iaux3(C%iaux2(1)))
          allocate(C%iaux4(C%iaux2(1)))
          allocate(C%zval(C%iaux2(1),1))
          C%iaux3(1:C%iaux2(1)-1)=iaux3_temp(1:C%iaux2(1)-1)
          C%iaux4(1:C%iaux2(1)-1)=iaux4_temp(1:C%iaux2(1)-1)
          C%zval(1:C%iaux2(1)-1,1)=zval_temp(1:C%iaux2(1)-1,1)
          deallocate(zval_temp)
          deallocate(iaux4_temp)
          deallocate(iaux3_temp)
          C%iaux3(C%iaux2(1))=i
          C%iaux4(C%iaux2(1))=j
          C%zval(C%iaux2(1),1)=alpha
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
          else if (label .eq. 'psp') then
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
#if defined(MPI) && defined(SLAP)
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
          else if (label .eq. 'psp') then
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
#if defined(MPI) && defined(SLAP)
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
    m_name%is_sparse=.false.

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
    m_name%is_sparse=.false.

    m_name%zval => A

    m_name%is_initialized=.true.

  end subroutine m_register_pzdbc
#endif

  !======================================================!
  ! register matrix by thresholding                      !
  ! parallel distributed 2D block cyclic sparse matrix   !
  !======================================================!
#ifdef PSP
  subroutine m_register_pdsp_thre(m_name,A,desc,spm_storage,thre)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9) ! BLACS array descriptor
    character(3), intent(in), target :: spm_storage ! storage format of sparse matrices, 'coo' or 'csc'
    real(dp), intent(in), target :: A(:,:) ! two-dimensional array containing the local matrix elements

    !**** OPTIONAL INPUT ***********************************!
    real(dp), optional :: thre ! non-negative threshold
    ! If thre=0, generate a sparse matrix with nonzero entries.
    ! If thre>0, generate a sparse matrix with entries with an absolute value >= thre.

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name ! matrix to be allocated

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!
    if (m_name%is_initialized .EQV. .true.) then
       call m_deallocate(m_name)
    end if
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
    m_name%str_type=spm_storage
    m_name%is_serial=.false.
    m_name%is_real=.true.
    m_name%is_sparse=.true.
    m_name%is_initialized=.true.

    call psp_den2sp_m(A,desc,m_name%spm,spm_storage,thre)

  end subroutine m_register_pdsp_thre
#endif

#ifdef PSP
  subroutine m_register_pzsp_thre(m_name,A,desc,spm_storage,thre)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9) ! BLACS array descriptor
    character(3), intent(in), target :: spm_storage ! storage format of sparse matrices, 'coo' or 'csc'
    complex(dp), intent(in), target :: A(:,:) ! two-dimensional array containing the local matrix elements

    !**** OPTIONAL INPUT ***********************************!
    real(dp), optional :: thre ! non-negative threshold
    ! If thre=0, generate a sparse matrix with nonzero entries.
    ! If thre>0, generate a sparse matrix with entries with an absolute value >= thre.

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name ! matrix to be allocated

    !**** INTERNAL ********************************!

    integer :: dim(2)

    !**********************************************!
    if (m_name%is_initialized .EQV. .true.) then
       call m_deallocate(m_name)
    end if
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
    m_name%str_type=spm_storage
    m_name%is_serial=.false.
    m_name%is_real=.false.
    m_name%is_sparse=.true.

    m_name%is_initialized=.true.

    call psp_den2sp_m(A,desc,m_name%spm,spm_storage,thre)

  end subroutine m_register_pzsp_thre
#endif

  !======================================================!
  ! register matrix using the Sparse Triplet format      !
  ! parallel distributed 2D block cyclic sparse matrix   !
  !======================================================!
#ifdef PSP
  subroutine m_register_pdsp_st(m_name,idx1,idx2,val,desc,spm_storage,nprow,npcol)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9) ! BLACS array descriptor
    character(3), intent(in), target :: spm_storage ! storage format of sparse matrices, 'coo' or 'csc'
    integer, intent(in), target :: idx1(:) ! one-dimensional array for row indices, local
    integer, intent(in), target :: idx2(:) ! one-dimensional array for column indices in 'coo' format or column pointers in 'csc' format, local
    real(dp), intent(in), target :: val(:) ! one-dimensional array containing the nonzero local matrix elements
    integer :: nprow, npcol

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name ! matrix to be allocated

    !**** INTERNAL ********************************!

    integer :: dim(2), iprow, ipcol

    !**********************************************!

    integer, external :: numroc

    !***** COMMON BLOCK ***************************!
    integer :: psp_bs_def_row, psp_bs_def_col, psp_icontxt

    common /coeff/ psp_bs_def_row, psp_bs_def_col, psp_icontxt

    !**********************************************!
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
    if (m_name%is_initialized .EQV. .true.) then
       call m_deallocate(m_name)
    end if
    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    dim(1)=numroc(m_name%dim1,psp_bs_def_row,iprow,0,nprow)
    dim(2)=numroc(m_name%dim2,psp_bs_def_col,ipcol,0,npcol)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type=spm_storage
    m_name%is_serial=.false.
    m_name%is_real=.true.
    m_name%is_sparse=.true.
    m_name%is_initialized=.true.

    call psp_register_spm(m_name%spm,idx1,idx2,val,desc,spm_storage,dim,nprow,npcol)

  end subroutine m_register_pdsp_st
#endif

#ifdef PSP
  subroutine m_register_pzsp_st(m_name,idx1,idx2,val,desc,spm_storage,nprow,npcol)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in), target :: desc(9) ! BLACS array descriptor
    character(3), intent(in), target :: spm_storage ! storage format of sparse matrices, 'coo' or 'csc'
    integer, intent(in), target :: idx1(:) ! one-dimensional array for row indices, local
    integer, intent(in), target :: idx2(:) ! one-dimensional array for column indices in 'coo' format or column pointers in 'csc' format, local
    complex(dp), intent(in), target :: val(:) ! one-dimensional array containing the nonzero local matrix elements
    integer :: nprow, npcol

    !**** INOUT ***********************************!

    type(matrix), intent(inout) :: m_name ! matrix to be allocated

    !**** INTERNAL ********************************!

    integer :: dim(2), iprow, ipcol

    !**********************************************!

    integer, external :: numroc

    !***** COMMON BLOCK ***************************!
    integer :: psp_bs_def_row, psp_bs_def_col, psp_icontxt

    common /coeff/ psp_bs_def_row, psp_bs_def_col, psp_icontxt

    !**********************************************!
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
    if (m_name%is_initialized .EQV. .true.) then
       call m_deallocate(m_name)
    end if
    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    dim(1)=numroc(m_name%dim1,psp_bs_def_row,iprow,0,nprow)
    dim(2)=numroc(m_name%dim2,psp_bs_def_col,ipcol,0,npcol)
    m_name%iaux2(1)=dim(1)
    m_name%iaux2(2)=dim(2)
    if (m_name%dim1==m_name%dim2) then
       m_name%is_square=.true.
    else
       m_name%is_square=.false.
    end if
    m_name%str_type=spm_storage
    m_name%is_serial=.false.
    m_name%is_real=.true.
    m_name%is_sparse=.true.
    m_name%is_initialized=.true.

    call psp_register_spm(m_name%spm,idx1,idx2,val,desc,spm_storage,dim,nprow,npcol)

  end subroutine m_register_pzsp_st
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

#ifdef PSP
  subroutine mm_multiply_pddbcpdcscref(A,trA,B,trB,C,alpha,beta)
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

    C%dval=beta*C%dval

    do n_comm=0,ms_mpi_size-1
       if (n_comm==ms_mpi_rank) then
          loc_dim_recv=B%spm%loc_dim2
          nnz_recv=B%spm%nnz
       end if
       call mpi_bcast(loc_dim_recv,1,mpi_integer,n_comm,mpi_comm_world,info)
       call mpi_bcast(nnz_recv,    1,mpi_integer,n_comm,mpi_comm_world,info)
       allocate(col_ptr_recv(loc_dim_recv+1))
       allocate(row_ind_recv(nnz_recv))
       allocate(dval_recv(nnz_recv))
       if (n_comm==ms_mpi_rank) then
          col_ptr_recv(1:loc_dim_recv+1)=B%spm%col_ptr(1:loc_dim_recv+1)
          row_ind_recv(1:nnz_recv)=B%spm%row_ind(1:nnz_recv)
          dval_recv(1:nnz_recv)=B%spm%dval(1:nnz_recv)
       end if
       call mpi_bcast(col_ptr_recv(1),loc_dim_recv+1,mpi_integer,         n_comm,mpi_comm_world,info)
       call mpi_bcast(row_ind_recv(1),nnz_recv,      mpi_integer,         n_comm,mpi_comm_world,info)
       call mpi_bcast(dval_recv(1),   nnz_recv,      mpi_double_precision,n_comm,mpi_comm_world,info)
       do i=1,loc_dim_recv
          do j=0,col_ptr_recv(i+1)-col_ptr_recv(i)-1
             l=col_ptr_recv(i)+j
             m=indxl2g(row_ind_recv(l),B%spm%desc(6),n_comm,B%spm%desc(8),ms_mpi_size)
             C%dval(:,m)=C%dval(:,m)+alpha*A%dval(:,i)*dval_recv(l)
          end do
       end do
       deallocate(dval_recv)
       deallocate(row_ind_recv)
       deallocate(col_ptr_recv)
    end do

  end subroutine mm_multiply_pddbcpdcscref

  subroutine mm_multiply_pdcscpddbcref(A,trA,B,trB,C,alpha,beta)
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

    C%dval=beta*C%dval

    do n_comm=0,ms_mpi_size-1
       if (n_comm==ms_mpi_rank) then
          loc_dim_recv=A%spm%loc_dim2
          nnz_recv=A%spm%nnz
       end if
       call mpi_bcast(loc_dim_recv,1,mpi_integer,n_comm,mpi_comm_world,info)
       call mpi_bcast(nnz_recv,    1,mpi_integer,n_comm,mpi_comm_world,info)
       allocate(col_ptr_recv(loc_dim_recv+1))
       allocate(row_ind_recv(nnz_recv))
       allocate(dval_recv(nnz_recv))
       if (n_comm==ms_mpi_rank) then
          col_ptr_recv(1:loc_dim_recv+1)=A%spm%col_ptr(1:loc_dim_recv+1)
          row_ind_recv(1:nnz_recv)=A%spm%row_ind(1:nnz_recv)
          dval_recv(1:nnz_recv)=A%spm%dval(1:nnz_recv)
       end if
       call mpi_bcast(col_ptr_recv(1),loc_dim_recv+1,mpi_integer,         n_comm,mpi_comm_world,info)
       call mpi_bcast(row_ind_recv(1),nnz_recv,      mpi_integer,         n_comm,mpi_comm_world,info)
       call mpi_bcast(dval_recv(1),   nnz_recv,      mpi_double_precision,n_comm,mpi_comm_world,info)
       do i=1,loc_dim_recv
          do j=0,col_ptr_recv(i+1)-col_ptr_recv(i)-1
             l=col_ptr_recv(i)+j
             m=indxl2g(i,A%spm%desc(6),n_comm,A%spm%desc(8),ms_mpi_size)
             C%dval(m,:)=C%dval(m,:)+alpha*dval_recv(l)*B%dval(row_ind_recv(l),:)
          end do
       end do
       deallocate(dval_recv)
       deallocate(row_ind_recv)
       deallocate(col_ptr_recv)
    end do

  end subroutine mm_multiply_pdcscpddbcref

  subroutine mm_multiply_pzdbcpzcscref(A,tcA,B,tcB,C,alpha,beta)
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

    C%zval=beta*C%zval

    do n_comm=0,ms_mpi_size-1
       if (n_comm==ms_mpi_rank) then
          loc_dim_recv=B%spm%loc_dim2
          nnz_recv=B%spm%nnz
       end if
       call mpi_bcast(loc_dim_recv,1,mpi_integer,n_comm,mpi_comm_world,info)
       call mpi_bcast(nnz_recv,    1,mpi_integer,n_comm,mpi_comm_world,info)
       allocate(col_ptr_recv(loc_dim_recv+1))
       allocate(row_ind_recv(nnz_recv))
       allocate(zval_recv(nnz_recv))
       if (n_comm==ms_mpi_rank) then
          col_ptr_recv(1:loc_dim_recv+1)=B%spm%col_ptr(1:loc_dim_recv+1)
          row_ind_recv(1:nnz_recv)=B%spm%row_ind(1:nnz_recv)
          zval_recv(1:nnz_recv)=B%spm%zval(1:nnz_recv)
       end if
       call mpi_bcast(col_ptr_recv(1),loc_dim_recv+1,mpi_integer,       n_comm,mpi_comm_world,info)
       call mpi_bcast(row_ind_recv(1),nnz_recv,      mpi_integer,       n_comm,mpi_comm_world,info)
       call mpi_bcast(zval_recv(1),   nnz_recv,      mpi_double_complex,n_comm,mpi_comm_world,info)
       do i=1,loc_dim_recv
          do j=0,col_ptr_recv(i+1)-col_ptr_recv(i)-1
             l=col_ptr_recv(i)+j
             m=indxl2g(row_ind_recv(l),B%spm%desc(6),n_comm,B%spm%desc(8),ms_mpi_size)
             C%zval(:,m)=C%zval(:,m)+alpha*A%zval(:,i)*conjg(zval_recv(l))
          end do
       end do
       deallocate(zval_recv)
       deallocate(row_ind_recv)
       deallocate(col_ptr_recv)
    end do

  end subroutine mm_multiply_pzdbcpzcscref

  subroutine mm_multiply_pzcscpzdbcref(A,tcA,B,tcB,C,alpha,beta)
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

    C%zval=beta*C%zval

    do n_comm=0,ms_mpi_size-1
       if (n_comm==ms_mpi_rank) then
          loc_dim_recv=A%spm%loc_dim2
          nnz_recv=A%spm%nnz
       end if
       call mpi_bcast(loc_dim_recv,1,mpi_integer,n_comm,mpi_comm_world,info)
       call mpi_bcast(nnz_recv,    1,mpi_integer,n_comm,mpi_comm_world,info)
       allocate(col_ptr_recv(loc_dim_recv+1))
       allocate(row_ind_recv(nnz_recv))
       allocate(zval_recv(nnz_recv))
       if (n_comm==ms_mpi_rank) then
          col_ptr_recv(1:loc_dim_recv+1)=A%spm%col_ptr(1:loc_dim_recv+1)
          row_ind_recv(1:nnz_recv)=A%spm%row_ind(1:nnz_recv)
          zval_recv(1:nnz_recv)=A%spm%zval(1:nnz_recv)
       end if
       call mpi_bcast(col_ptr_recv(1),loc_dim_recv+1,mpi_integer,       n_comm,mpi_comm_world,info)
       call mpi_bcast(row_ind_recv(1),nnz_recv,      mpi_integer,       n_comm,mpi_comm_world,info)
       call mpi_bcast(zval_recv(1),   nnz_recv,      mpi_double_complex,n_comm,mpi_comm_world,info)
       do i=1,loc_dim_recv
          do j=0,col_ptr_recv(i+1)-col_ptr_recv(i)-1
             l=col_ptr_recv(i)+j
             m=indxl2g(i,A%spm%desc(6),n_comm,A%spm%desc(8),ms_mpi_size)
             C%zval(m,:)=C%zval(m,:)+alpha*conjg(zval_recv(l))*B%zval(row_ind_recv(l),:)
          end do
       end do
       deallocate(zval_recv)
       deallocate(row_ind_recv)
       deallocate(col_ptr_recv)
    end do

  end subroutine mm_multiply_pzcscpzdbcref
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
  subroutine ms_scalapack_setup(mpi_rank,mpi_size,nprow,order,bs_def,bs_list,icontxt)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: order ! ordering of processor grid: 'r/R' or other for row-major, 'c/C' for column-major

    integer, intent(in) :: mpi_rank ! rank of local MPI process
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
    ms_mpi_rank=mpi_rank
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

#ifdef PSP
    ! initialized grid information in pspBLAS
    call psp_gridinit(mpi_size,nprow,order,bs_def,bs_def,icontxt)
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

    !**** INTERNAL ********************************!

    logical, save :: log_start=.false.

    integer :: log_unit

    !**********************************************!

    if (log_start) then
       open(newunit=log_unit,file='MatrixSwitch.log',position='append')
    else
       open(newunit=log_unit,file='MatrixSwitch.log',status='replace')
       log_start=.true.
    end if
    write(log_unit,'(a)'), 'FATAL ERROR in matrix_switch!'
    write(log_unit,'(a)'), message
    close(log_unit)
    stop

  end subroutine die

end module MatrixSwitch
