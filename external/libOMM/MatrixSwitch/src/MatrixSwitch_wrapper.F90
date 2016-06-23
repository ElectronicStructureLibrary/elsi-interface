module MatrixSwitch_wrapper
  use MatrixSwitch_wrapper_params
  use MatrixSwitch, only: &
    mm_multiply_orig => mm_multiply, &
    m_add_orig => m_add, &
    m_trace_orig => m_trace, &
    mm_trace_orig => mm_trace, &
    m_scale_orig => m_scale, &
    m_set_orig => m_set, &
    m_set_element_orig => m_set_element, &
    m_get_element_orig => m_get_element, &
    m_register_sden_orig => m_register_sden, &
    m_allocate_orig => m_allocate, &
    m_deallocate_orig => m_deallocate, &
    m_copy_orig => m_copy, &
    m_convert_orig => m_convert
#ifdef MPI
  use MatrixSwitch, only: &
    m_register_pdbc_orig => m_register_pdbc, &
    ms_scalapack_setup, &
    ms_lap_icontxt
#endif
#ifdef PSP
  use MatrixSwitch, only: &
    m_register_psp_thre_orig => m_register_psp_thre, &
    m_register_psp_st_orig => m_register_psp_st
#endif

  implicit none

  private

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300)

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

  public :: ms_wrapper_open
  public :: ms_wrapper_close
  public :: ms_is_initialized
  public :: ms_is_serial
  public :: ms_is_real
  public :: ms_is_square
  public :: ms_is_sparse
  public :: ms_dim1
  public :: ms_dim2
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
  ! setup the MatrixSwitch wrapper                 !
  !================================================!
  subroutine ms_wrapper_open(num_matrices,keys)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: keys(num_matrices) ! mapping keys

    integer, intent(in) :: num_matrices ! number of matrices to use

    !**********************************************!

    allocate(ms_matrices(num_matrices))
    allocate(ms_keys(num_matrices))

    ms_num_matrices=num_matrices
    ms_keys=keys

  end subroutine ms_wrapper_open

  !================================================!
  ! close the MatrixSwitch wrapper                 !
  !================================================!
  subroutine ms_wrapper_close()
    implicit none

    !**** INTERNAL ********************************!

    integer :: i

    !**********************************************!

    if (allocated(ms_keys)) deallocate(ms_keys)
    if (allocated(ms_matrices)) then
      do i=1,ms_num_matrices
        call m_deallocate_orig(ms_matrices(i))
      end do
      deallocate(ms_matrices)
    end if

  end subroutine ms_wrapper_close

  !================================================!
  ! functions to unpack the matrix data type       !
  !================================================!
  logical function ms_is_initialized(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name ! matrix to unpack

    !**********************************************!

    ms_is_initialized=ms_matrices(ms_lookup(m_name))%is_initialized

  end function ms_is_initialized

  logical function ms_is_serial(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name ! matrix to unpack

    !**********************************************!

    ms_is_serial=ms_matrices(ms_lookup(m_name))%is_serial

  end function ms_is_serial

  logical function ms_is_real(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name ! matrix to unpack

    !**********************************************!

    ms_is_real=ms_matrices(ms_lookup(m_name))%is_real

  end function ms_is_real

  logical function ms_is_square(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name ! matrix to unpack

    !**********************************************!

    ms_is_square=ms_matrices(ms_lookup(m_name))%is_square

  end function ms_is_square

  logical function ms_is_sparse(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name ! matrix to unpack

    !**********************************************!

    ms_is_sparse=ms_matrices(ms_lookup(m_name))%is_sparse

  end function ms_is_sparse

  integer function ms_dim1(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name ! matrix to unpack

    !**********************************************!

    ms_dim1=ms_matrices(ms_lookup(m_name))%dim1

  end function ms_dim1

  integer function ms_dim2(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name ! matrix to unpack

    !**********************************************!

    ms_dim2=ms_matrices(ms_lookup(m_name))%dim2

  end function ms_dim2

  !================================================!
  ! allocate matrix                                !
  !================================================!
  subroutine m_allocate(m_name,i,j,label)
    implicit none

    !**** INPUT ***********************************!

    character(5), intent(in), optional :: label ! storage format to use (see documentation)
    character(*), intent(in) :: m_name ! matrix to be allocated

    integer, intent(in) :: i ! (global) row dimension size of the matrix
    integer, intent(in) :: j ! (global) column dimension size of the matrix

    !**********************************************!

    call m_allocate_orig(ms_matrices(ms_lookup(m_name)),i,j,label)

  end subroutine m_allocate

  !================================================!
  ! deallocate matrix                              !
  !================================================!
  subroutine m_deallocate(m_name)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: m_name ! matrix to be deallocated

    !**********************************************!

    call m_deallocate_orig(ms_matrices(ms_lookup(m_name)))

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

    character(*), intent(in) :: A ! matrix to copy from
    character(*), intent(in) :: m_name ! matrix to copy onto

    !**********************************************!

    call m_copy_orig(ms_matrices(ms_lookup(m_name)),ms_matrices(ms_lookup(A)),label,threshold,threshold_is_soft)

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

    character(*), intent(in) :: m_name ! matrix to convert

    !**********************************************!

    call m_convert_orig(ms_matrices(ms_lookup(m_name)),label,threshold,threshold_is_soft)

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

    character(*), intent(in) :: A ! matrix A
    character(*), intent(in) :: B ! matrix B
    character(*), intent(in) :: C ! matrix C

    !**********************************************!

    call mm_multiply_orig(ms_matrices(ms_lookup(A)),opA,ms_matrices(ms_lookup(B)),opB,ms_matrices(ms_lookup(C)),alpha,beta,label)

  end subroutine mm_dmultiply

  subroutine mm_zmultiply(A,opA,B,opB,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T' for A^T, 'c/C' for A^H
    character(1), intent(in) :: opB ! form of op(B)
    character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

    complex(dp), intent(in) :: alpha ! scalar alpha
    complex(dp), intent(in) :: beta ! scalar beta

    character(*), intent(in) :: A ! matrix A
    character(*), intent(in) :: B ! matrix B
    character(*), intent(in) :: C ! matrix C

    !**********************************************!

    call mm_multiply_orig(ms_matrices(ms_lookup(A)),opA,ms_matrices(ms_lookup(B)),opB,ms_matrices(ms_lookup(C)),alpha,beta,label)

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

    character(*), intent(in) :: A ! matrix A
    character(*), intent(in) :: C ! matrix C

    !**********************************************!

    call m_add_orig(ms_matrices(ms_lookup(A)),opA,ms_matrices(ms_lookup(C)),alpha,beta,label)

  end subroutine m_dadd

  subroutine m_zadd(A,opA,C,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T' for A^T, 'c/C' for A^H
    character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

    complex(dp), intent(in) :: alpha ! scalar alpha
    complex(dp), intent(in) :: beta ! scalar beta

    character(*), intent(in) :: A ! matrix A
    character(*), intent(in) :: C ! matrix C

    !**********************************************!

    call m_add_orig(ms_matrices(ms_lookup(A)),opA,ms_matrices(ms_lookup(C)),alpha,beta,label)

  end subroutine m_zadd

  !================================================!
  ! matrix trace                                   !
  ! alpha := tr(A)                                 !
  !================================================!
  subroutine m_dtrace(A,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

    character(*), intent(in) :: A ! matrix A

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha ! scalar alpha

    !**********************************************!

    call m_trace_orig(ms_matrices(ms_lookup(A)),alpha,label)

  end subroutine m_dtrace

  subroutine m_ztrace(A,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

    character(*), intent(in) :: A ! matrix A

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha ! scalar alpha

    !**********************************************!

    call m_trace_orig(ms_matrices(ms_lookup(A)),alpha,label)

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

    character(*), intent(in) :: A ! matrix A
    character(*), intent(in) :: B ! matrix B

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha ! scalar alpha

    !**********************************************!

    call mm_trace_orig(ms_matrices(ms_lookup(A)),ms_matrices(ms_lookup(B)),alpha,label)

  end subroutine mm_dtrace

  subroutine mm_ztrace(A,B,alpha,label)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

    character(*), intent(in) :: A ! matrix A
    character(*), intent(in) :: B ! matrix B

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha ! scalar alpha

    !**********************************************!

    call mm_trace_orig(ms_matrices(ms_lookup(A)),ms_matrices(ms_lookup(B)),alpha,label)

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

    character(*), intent(in) :: C ! matrix C

    !**********************************************!

    call m_scale_orig(ms_matrices(ms_lookup(C)),beta,label)

  end subroutine m_dscale

  subroutine m_zscale(C,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

    complex(dp), intent(in) :: beta ! scalar beta

    character(*), intent(in) :: C ! matrix C

    !**********************************************!

    call m_scale_orig(ms_matrices(ms_lookup(C)),beta,label)

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

    character(*), intent(in) :: C ! matrix C

    !**********************************************!

    call m_set_orig(ms_matrices(ms_lookup(C)),seC,alpha,beta,label)

  end subroutine m_dset

  subroutine m_zset(C,seC,alpha,beta,label)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: seC ! part of matrix to set: 'l/L' for lower, 'u/U' for upper, other for all
    character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

    complex(dp), intent(in) :: alpha ! scalar alpha
    complex(dp), intent(in) :: beta ! scalar beta

    character(*), intent(in) :: C ! matrix C

    !**********************************************!

    call m_set_orig(ms_matrices(ms_lookup(C)),seC,alpha,beta,label)

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

    character(*), intent(in) :: C ! matrix C

    !**********************************************!

    call m_set_element_orig(ms_matrices(ms_lookup(C)),i,j,alpha,label)

  end subroutine m_dset_element

  subroutine m_zset_element(C,i,j,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

    integer, intent(in) :: i ! row index of element
    integer, intent(in) :: j ! column index of element

    complex(dp), intent(in) :: alpha ! scalar alpha

    character(*), intent(in) :: C ! matrix C

    !**********************************************!

    call m_set_element_orig(ms_matrices(ms_lookup(C)),i,j,alpha,label)

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

    character(*), intent(in) :: C ! matrix C

    !**** OUTPUT **********************************!

    real(dp), intent(out) :: alpha ! scalar alpha

    !**********************************************!

    call m_get_element_orig(ms_matrices(ms_lookup(C)),i,j,alpha,label)

  end subroutine m_dget_element

  subroutine m_zget_element(C,i,j,alpha,label)
    implicit none

    !**** INPUT ***********************************!

    character(3), intent(in), optional :: label ! implementation of the operation to use (see documentation)

    integer, intent(in) :: i ! row index of element
    integer, intent(in) :: j ! column index of element

    character(*), intent(in) :: C ! matrix C

    !**** OUTPUT **********************************!

    complex(dp), intent(out) :: alpha ! scalar alpha

    !**********************************************!

    call m_get_element_orig(ms_matrices(ms_lookup(C)),i,j,alpha,label)

  end subroutine m_zget_element

  !================================================!
  ! register matrix: simple dense serial           !
  !================================================!
  subroutine m_register_sdden(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    real(dp), intent(in), target :: A(:,:) ! two-dimensional array containing the matrix elements

    character(*), intent(in) :: m_name ! matrix to be allocated

    !**********************************************!

    call m_register_sden_orig(ms_matrices(ms_lookup(m_name)),A)

  end subroutine m_register_sdden

  subroutine m_register_szden(m_name,A)
    implicit none

    !**** INPUT ***********************************!

    complex(dp), intent(in), target :: A(:,:) ! two-dimensional array containing the matrix elements

    character(*), intent(in) :: m_name ! matrix to be allocated

    !**********************************************!

    call m_register_sden_orig(ms_matrices(ms_lookup(m_name)),A)

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

    character(*), intent(in) :: m_name ! matrix to be allocated

    !**********************************************!

    call m_register_pdbc_orig(ms_matrices(ms_lookup(m_name)),A,desc)

  end subroutine m_register_pddbc
#endif

#ifdef MPI
  subroutine m_register_pzdbc(m_name,A,desc)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: desc(9) ! BLACS array descriptor

    complex(dp), intent(in), target :: A(:,:) ! two-dimensional array containing the local matrix elements

    character(*), intent(in) :: m_name ! matrix to be allocated

    !**********************************************!

    call m_register_pdbc_orig(ms_matrices(ms_lookup(m_name)),A,desc)

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

    character(*), intent(in) :: m_name ! matrix to be allocated

    !**********************************************!

    call m_register_psp_thre_orig(ms_matrices(ms_lookup(m_name)),A,desc,spm_storage,thre)

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

    character(*), intent(in) :: m_name ! matrix to be allocated

    !**********************************************!

    call m_register_psp_thre_orig(ms_matrices(ms_lookup(m_name)),A,desc,spm_storage,thre)

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

    character(*), intent(in) :: m_name ! matrix to be allocated

    !**********************************************!

    call m_register_psp_st_orig(ms_matrices(ms_lookup(m_name)),idx1,idx2,val,desc,spm_storage,nprow,npcol)

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

    character(*), intent(in) :: m_name ! matrix to be allocated

    !**********************************************!

    call m_register_psp_st_orig(ms_matrices(ms_lookup(m_name)),idx1,idx2,val,desc,spm_storage,nprow,npcol)

  end subroutine m_register_pzsp_st
#endif

end module MatrixSwitch_wrapper