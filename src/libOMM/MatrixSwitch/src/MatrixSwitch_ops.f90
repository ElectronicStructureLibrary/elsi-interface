module MatrixSwitch_ops

  implicit none

  public

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** VARIABLES *********************************!

  character(1), save :: ms_lap_order

  integer, save :: ms_mpi_comm
  integer, save :: ms_mpi_size
  integer, save :: ms_mpi_rank
  integer, save :: ms_lap_nprow
  integer, save :: ms_lap_npcol
  integer, save :: ms_lap_bs_def
  integer, save :: ms_lap_bs_num
  integer, save :: ms_lap_icontxt ! BLACS context handle used by MatrixSwitch
  integer, allocatable, save :: ms_lap_bs_list(:,:)

  !**** TYPES *************************************!

  ! This is the derived type that encapsulates all matrix storage possibilities and hides the details from the user.
  type matrix
     character(3) :: str_type ! label identifying the storage format

     logical :: is_initialized=.false. ! has the matrix been initialized?
     logical :: is_serial ! is the matrix serial or parallel distributed?
     logical :: is_real ! is the matrix real or complex (both kind dp)?
     logical :: is_square ! is the matrix square?
     logical :: is_sparse ! is the matrix sparse?
     logical :: iaux1_is_allocated=.false. ! is iaux1 directly allocated?
     logical :: iaux2_is_allocated=.false. ! is iaux2 directly allocated?
     logical :: iaux3_is_allocated=.false. ! is iaux3 directly allocated?
     logical :: iaux4_is_allocated=.false. ! is iaux4 directly allocated?
     logical :: dval_is_allocated=.false. ! is dval directly allocated?
     logical :: zval_is_allocated=.false. ! is zval directly allocated?

     integer :: dim1 ! (global) row dimension size of the matrix
     integer :: dim2 ! (global) column dimension size of the matrix
     integer, pointer :: iaux1(:) => null() ! auxiliary information for certain storage formats
     integer, pointer :: iaux2(:) => null() ! auxiliary information for certain storage formats
     integer, pointer :: iaux3(:) => null() ! auxiliary information for certain storage formats
     integer, pointer :: iaux4(:) => null() ! auxiliary information for certain storage formats

     real(dp), pointer :: dval(:,:) => null() ! matrix elements for a real matrix

     complex(dp), pointer :: zval(:,:) => null() ! matrix elements for a complex matrix

  end type matrix

  !**** INTERFACES ********************************!

  interface process_opM
     module procedure process_lopM
     module procedure process_iopM
  end interface process_opM

  !************************************************!

contains

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

  subroutine die(message)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in), optional :: message

    !**** INTERNAL ********************************!

    integer :: err_unit=377

    !**********************************************!

    open(unit=err_unit,file='MatrixSwitch.err')
    write(err_unit,'(a)') 'FATAL ERROR in matrix_switch!'
    if (present(message)) write(err_unit,'(a)') message
    write(err_unit,'(a,1x,i5)') 'MPI rank:', ms_mpi_rank
    close(err_unit)
    stop

  end subroutine die

end module MatrixSwitch_ops
