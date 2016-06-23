module MatrixSwitch_wrapper_params
  use MatrixSwitch, only: matrix

  implicit none

  private

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300)
  integer, parameter :: max_key_length=10

  !**** VARIABLES *********************************!

  character(max_key_length), allocatable :: ms_keys(:)

  integer :: ms_num_matrices

  type(matrix), allocatable :: ms_matrices(:)

  !************************************************!

  public :: matrix
  public :: ms_keys
  public :: ms_num_matrices
  public :: ms_matrices
  public :: ms_lookup

contains

  !================================================!
  ! map key to index of ms_matrices array          !
  !================================================!
  integer function ms_lookup(key)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in) :: key ! matrix key

    !**** INTERNAL ********************************!

    integer :: i

    !**********************************************!

    do i=1,ms_num_matrices
      if (key .eq. ms_keys(i)) then
        ms_lookup=i
        exit
      end if
    end do

  end function ms_lookup

end module MatrixSwitch_wrapper_params