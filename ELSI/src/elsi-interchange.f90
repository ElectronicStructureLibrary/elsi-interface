!> ELSI Interchange
!! 
!! Copyright of the original code rests with the authors inside the ELSI
!! consortium. The copyright of any additional modifications shall rest
!! with their original authors, but shall adhere to the licensing terms
!! distributed along with the original code in the file "COPYING".

module ELSI

  use iso_c_binding

  implicit none

  PRIVATE ! By default, all routines contained are private

  ! The following routines are public:

  public :: set_matrix_element

contains

function set_matrix_element(i_row, i_col, element)

!>
!!  This routine sets a single element of the matrix
!!
!!  Parameters
!!
!!  i_row      row of the matrix
!!
!!  i_col      column of the matrix
!!
!!  element    value to be set
!!

   implicit none

   !! empty so far

   return 0
end function

end module ELSI
