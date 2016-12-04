module MatrixSwitch_m_register
  use MatrixSwitch_ops

  implicit none

  !**** INTERFACES ********************************!

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

contains

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
    m_name%iaux2_is_allocated=.true.
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
    m_name%iaux1_is_allocated=.true.
    m_name%iaux1=desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
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
    !if (m_name%is_initialized .EQV. .true.) then
    !   call m_deallocate(m_name)
    !end if
    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
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
    !if (m_name%is_initialized .EQV. .true.) then
    !   call m_deallocate(m_name)
    !end if
    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
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
    !if (m_name%is_initialized .EQV. .true.) then
    !   call m_deallocate(m_name)
    !end if
    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
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
    !if (m_name%is_initialized .EQV. .true.) then
    !   call m_deallocate(m_name)
    !end if
    m_name%iaux1 => desc
    m_name%dim1=desc(3)
    m_name%dim2=desc(4)
    allocate(m_name%iaux2(2))
    m_name%iaux2_is_allocated=.true.
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

end module MatrixSwitch_m_register
