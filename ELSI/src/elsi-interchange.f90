!> ELSI Interchange
!! 
!! Copyright of the original code rests with the authors inside the ELSI
!! consortium. The copyright of any additional modifications shall rest
!! with their original authors, but shall adhere to the licensing terms
!! distributed along with the original code in the file "COPYING".

module ELSI

  use iso_c_binding
  use HDF5

  implicit none

  PRIVATE ! By default, all routines contained are private

  integer method = -1 !<Method for EV Solver (ELPA=1,OMM=2,PEXSI=3)
  integer mode = -1 !<Mode for EV Solver (REAL_VALUES=1,COMPLEX_VALUES=2)

  real*8, allocatable H_real, S_real 
  complex*8, allocatable H_complex, S_complex
  
  ! The following routines are public:


  public :: set_matrix_element

  enum, bind( C )
    enumerator :: ELPA, OMM, PEXSI
  end enum

  enum, bind( C )
    enumerator :: REAL_VALUES, COMPLEX_VALUES
  end enum   


contains

subroutine set_method(i_method)

   !>
   !!  This routine sets the method of choice for solving the eigenvalue problem
   !!

   implicit none

   integer, intent(in) :: i_method !<one of (ELPA,OMM,PEXSI)

   method = i_method

end subroutine set_method

function allocate_matrices(n_rows, n_cols)

!>
!!  This routine prepares the matrices based on the dimension and method
!!  specified

   implicit none

   integer, intent(in) :: n_rows !< Number of rows 
   integer, intent(in) :: n_cols !< Number of cols
   
   select case (method)
      case ELPA
         select case (mode)
            case COMPLEX_VALUES
               allocate(H_complex (n_rows, n_cols))      
               allocate(S_complex (n_rows, n_cols))
            case REAL_VALUES
               allocate(H_real (n_rows, n_cols))      
               allocate(S_real (n_rows, n_cols))
         end select
      case OMM
         write(*,'a') "OMM not implemented yet!"
         stop
      case PEXSI
         write(*,'a') "PEXSI not implemented yet!"
         stop
      case DEFAULT
         write(*,'2a') "No method has been chosen. " &
            "Please choose method ELPA, OMM, or PEXSI"
         stop
      end select

   return 0
end function

function deallocate_matrices()

!>
!!  This routine prepares the matrices based on the dimension and method
!!  specified

   implicit none

   select case (method)
      case ELPA
         select case (mode)
            case COMPLEX_VALUES
               deallocate(H_complex)      
               deallocate(S_complex)
            case REAL_VALUES
               deallocate(H_real)      
               deallocate(S_real)
         end select
      case OMM
         write(*,'a') "OMM not implemented yet!"
         stop
      case PEXSI
         write(*,'a') "PEXSI not implemented yet!"
         stop
      case DEFAULT
         write(*,'2a') "No method has been chosen. " &
            "Please choose method ELPA, OMM, or PEXSI"
         stop
      end select

   return 0
end function


end module ELSI
