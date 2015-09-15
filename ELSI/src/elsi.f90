!> ELSI Interchange
!! 
!! Copyright of the original code rests with the authors inside the ELSI
!! consortium. The copyright of any additional modifications shall rest
!! with their original authors, but shall adhere to the licensing terms
!! distributed along with the original code in the file "COPYING".

module ELSI

  use iso_c_binding
  use HDF5
  use ELPA1
  use ELPA2
  use MPI_TOOLS

  implicit none

  integer :: method = -1 !<Method for EV Solver (ELPA=1,OMM=2,PEXSI=3)
  integer :: mode = -1 !<Mode for EV Solver (REAL_VALUES=1,COMPLEX_VALUES=2)

  real*8, allocatable    :: H_real(:,:), S_real(:,:) 
  complex*8, allocatable :: H_complex(:,:), S_complex(:,:)
  
  enum, bind( C )
    enumerator :: ELPA, OMM, PEXSI
  end enum

  enum, bind( C )
    enumerator :: REAL_VALUES, COMPLEX_VALUES
  end enum   

  ! The following variables are private
  private :: H_real, S_real, H_complex, S_complex

  ! The following routines are public:

  public :: elsi_set_method          !< Set routine for method
  public :: elsi_set_mode            !< Set routine for mode
  public :: elsi_allocate_matrices   !< Initialize matrices
  public :: elsi_deallocate_matrices !< Cleanup matrix memory 

contains

subroutine elsi_set_method(i_method)

   !>
   !!  This routine sets the method of choice for solving the eigenvalue problem
   !!

   implicit none

   integer, intent(in) :: i_method !<one of (ELPA,OMM,PEXSI)

   method = i_method

end subroutine

subroutine elsi_set_mode(i_mode)

   !>
   !!  This routine sets the mode of the eigenvalue solver (Real or Complex)
   !!

   implicit none

   integer, intent(in) :: i_mode !<one of (REAL_VALUES,COMPLEX_VALUES)

   mode = i_mode


end subroutine 


subroutine elsi_allocate_matrices(n_rows, n_cols)

!>
!!  This routine prepares the matrices based on the dimension and method
!!  specified

   implicit none

   integer, intent(in) :: n_rows !< Number of rows 
   integer, intent(in) :: n_cols !< Number of cols
   
   select case (method)
      case (ELPA)
         select case (mode)
            case (COMPLEX_VALUES)
               allocate(H_complex (n_rows, n_cols))      
               allocate(S_complex (n_rows, n_cols))
            case (REAL_VALUES)
               allocate(H_real (n_rows, n_cols))      
               allocate(S_real (n_rows, n_cols))
            case DEFAULT
               write(*,'(2a)') "No mode has been chosen. ", &
               "Please choose method REAL_VALUES or COMPLEX_VALUES"
               stop

         end select
      case (OMM)
         write(*,'(a)') "OMM not implemented yet!"
         stop
      case (PEXSI)
         write(*,'(a)') "PEXSI not implemented yet!"
         stop
      case DEFAULT
         write(*,'(2a)') "No method has been chosen. ", &
            "Please choose method ELPA, OMM, or PEXSI"
         stop
   end select
end subroutine

subroutine elsi_set_real_hamiltonian(h,n_rows,n_cols)

!>
!!  This routine sets the real hamiltonian matrix 

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   real*8, intent(in) :: h(n_rows,n_cols) !< hamiltonian 
   
   select case (method)
      case (ELPA)
         if (mode == REAL_VALUES) then
            H_real(:,:) = h(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Complex valued hamiltonian to be written in real storage"
            stop
         end if
      case (OMM)
         write(*,'(a)') "OMM not implemented yet!"
         stop
      case (PEXSI)
         write(*,'(a)') "PEXSI not implemented yet!"
         stop
      case DEFAULT
         write(*,'(2a)') "No method has been chosen. ", &
            "Please choose method ELPA, OMM, or PEXSI"
         stop
   end select


end subroutine

subroutine elsi_set_complex_hamiltonian(h,n_rows,n_cols)

!>
!!  This routine sets the complex hamiltonian matrix 

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   real*8, intent(in) :: h(n_rows,n_cols) !< hamiltonian 
   
   select case (method)
      case (ELPA)
         if (mode == COMPLEX_VALUES) then
            H_complex(:,:) = h(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Real valued hamiltonian to be written in complex storage"
            stop
         end if
      case (OMM)
         write(*,'(a)') "OMM not implemented yet!"
         stop
      case (PEXSI)
         write(*,'(a)') "PEXSI not implemented yet!"
         stop
      case DEFAULT
         write(*,'(2a)') "No method has been chosen. ", &
            "Please choose method ELPA, OMM, or PEXSI"
         stop
   end select


end subroutine

subroutine elsi_set_real_overlap(s,n_rows,n_cols)

!>
!!  This routine sets the real overlap matrix 

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   real*8, intent(in) :: s(n_rows,n_cols) !< overlap 
   
   select case (method)
      case (ELPA)
         if (mode == REAL_VALUES) then
            S_real(:,:) = s(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Complex valued overlap to be written in real storage"
            stop
         end if
      case (OMM)
         write(*,'(a)') "OMM not implemented yet!"
         stop
      case (PEXSI)
         write(*,'(a)') "PEXSI not implemented yet!"
         stop
      case DEFAULT
         write(*,'(2a)') "No method has been chosen. ", &
            "Please choose method ELPA, OMM, or PEXSI"
         stop
   end select


end subroutine

subroutine elsi_set_complex_overlap(s,n_rows,n_cols)

!>
!!  This routine sets the hamiltonian matrix 

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   complex*8, intent(in) :: s(n_rows,n_cols) !< overlap 
   
   select case (method)
      case (ELPA)
         if (mode == REAL_VALUES) then
            S_complex(:,:) = s(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:", &
               " Real valued overlap to be written in complex storage"
            stop
         end if
      case (OMM)
         write(*,'(a)') "OMM not implemented yet!"
         stop
      case (PEXSI)
         write(*,'(a)') "PEXSI not implemented yet!"
         stop
      case DEFAULT
         write(*,'(2a)') "No method has been chosen. ", &
            "Please choose method ELPA, OMM, or PEXSI"
         stop
   end select


end subroutine

subroutine elsi_get_real_hamiltonian(h,n_rows,n_cols)

!>
!!  This routine gets the real hamiltonian matrix 

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   real*8, intent(out) :: h(n_rows,n_cols) !< hamiltonian 
   
   select case (method)
      case (ELPA)
         if (mode == COMPLEX_VALUES) then
            h(:,:) = H_real(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Real valued hamiltonian to be written in complex storage"
            stop
         end if
      case (OMM)
         write(*,'(a)') "OMM not implemented yet!"
         stop
      case (PEXSI)
         write(*,'(a)') "PEXSI not implemented yet!"
         stop
      case DEFAULT
         write(*,'(2a)') "No method has been chosen. ", &
            "Please choose method ELPA, OMM, or PEXSI"
         stop
   end select


end subroutine


subroutine elsi_get_complex_hamiltonian(h,n_rows,n_cols)

!>
!!  This routine gets the complex hamiltonian matrix 

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   complex*8, intent(out) :: h(n_rows,n_cols) !< hamiltonian 
   
   select case (method)
      case (ELPA)
         if (mode == COMPLEX_VALUES) then
            h(:,:) = H_complex(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Complex valued hamiltonian to be written in real storage"
            stop
         end if
      case (OMM)
         write(*,'(a)') "OMM not implemented yet!"
         stop
      case (PEXSI)
         write(*,'(a)') "PEXSI not implemented yet!"
         stop
      case DEFAULT
         write(*,'(2a)') "No method has been chosen. ", &
            "Please choose method ELPA, OMM, or PEXSI"
         stop
   end select


end subroutine

subroutine elsi_get_real_overlap(s,n_rows,n_cols)

!>
!!  This routine gets the real overlap matrix 

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   real*8, intent(out) :: s(n_rows,n_cols) !< overlap 
   
   select case (method)
      case (ELPA)
         if (mode == REAL_VALUES) then
            s(:,:) = S_real(:,:)    
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Real valued overlap to be written in complex storage"
            stop
         end if
      case (OMM)
         write(*,'(a)') "OMM not implemented yet!"
         stop
      case (PEXSI)
         write(*,'(a)') "PEXSI not implemented yet!"
         stop
      case DEFAULT
         write(*,'(2a)') "No method has been chosen. ", &
            "Please choose method ELPA, OMM, or PEXSI"
         stop
   end select


end subroutine

subroutine elsi_get_complex_overlap(s,n_rows,n_cols)

!>
!!  This routine gets the complex overlap matrix 

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   complex*8, intent(out) :: s(n_rows,n_cols) !< overlap 
   
   select case (method)
      case (ELPA)
         if (mode == COMPLEX_VALUES) then
            s(:,:) = S_complex(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Complex valued overlap to be written in real storage"
            stop
         end if
      case (OMM)
         write(*,'(a)') "OMM not implemented yet!"
         stop
      case (PEXSI)
         write(*,'(a)') "PEXSI not implemented yet!"
         stop
      case DEFAULT
         write(*,'(2a)') "No method has been chosen. ", &
            "Please choose method ELPA, OMM, or PEXSI"
         stop
   end select


end subroutine

subroutine elsi_write_ev_problem(filename)

!>
!!  This routine writes the complete eigenvalue problem into a file

   implicit none
   character(len=*), intent(in) :: filename

   integer :: h5err
   integer(hid_t) :: plist_id
   integer(hid_t) :: file_id

   call h5open_f(h5err)
   if (h5err) then
      write(*,'(a)') "HDF5: Failed to initialize Fortran interface."
      stop
   end if

   call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
   if (h5err) then
      write(*,'(a)') "HDF5: Failed to create property list for parallel access."
      stop
   end if

   call h5pset_fapl_mpio_f(plist_id, mpi_comm_world, mpi_info_null, h5err)
   if (h5err) then
      write(*,'(a)') "HDF5: Failed to set property list for parallel access."
      stop
   end if

   call h5fcreate_f( trim(filename),H5F_ACC_TRUNC_F, file_id, h5err, & 
                     access_prp = plist_id)
   if (h5err) then
      write(*,'(a)') "HDF5: Failed to create file for parallel access."
      stop
   end if

   call h5fclose_f(file_id, h5err)
   if (h5err) then
      write(*,'(a)') "HDF5: Failed to close file for parallel access."
      stop
   end if

   call h5close_f(h5err)
   if (h5err) then
      write(*,'(a)') "HDF5: Failed to finalize Fortran interface."
      stop
   end if

end subroutine

subroutine elsi_read_ev_problem(filename)

!>
!!  This routine writes the complete eigenvalue problem into a file

   implicit none
   character(len=*), intent(in) :: filename

   integer :: h5err
   integer(hid_t) :: plist_id
   integer(hid_t) :: file_id

   call h5open_f(h5err)
   if (h5err) then
      write(*,'(a)') "HDF5: Failed to initialize Fortran interface."
      stop
   end if

   call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
   if (h5err) then
      write(*,'(a)') "HDF5: Failed to create property list for parallel access."
      stop
   end if

   call h5pset_fapl_mpio_f(plist_id, mpi_comm_world, mpi_info_null, h5err)
   if (h5err) then
      write(*,'(a)') "HDF5: Failed to set property list for parallel access."
      stop
   end if

   call h5fcreate_f( filename,H5F_ACC_TRUNC_F, file_id, h5err, & 
                     access_prp = plist_id)
   if (h5err) then
      write(*,'(a)') "HDF5: Failed to create file for parallel access."
      stop
   end if

   call h5fclose_f(file_id, h5err)
   if (h5err) then
      write(*,'(a)') "HDF5: Failed to close file for parallel access."
      stop
   end if

   call h5close_f(h5err)
   if (h5err) then
      write(*,'(a)') "HDF5: Failed to finalize Fortran interface."
      stop
   end if

end subroutine


subroutine elsi_deallocate_matrices()

!>
!!  This routine prepares the matrices based on the dimension and method
!!  specified

   implicit none

   select case (method)
      case (ELPA)
         select case (mode)
            case (COMPLEX_VALUES)
               deallocate(H_complex)      
               deallocate(S_complex)
            case (REAL_VALUES)
               deallocate(H_real)      
               deallocate(S_real)
         end select
      case (OMM)
         write(*,'(a)') "OMM not implemented yet!"
         stop
      case (PEXSI)
         write(*,'(a)') "PEXSI not implemented yet!"
         stop
      case DEFAULT
         write(*,'(2a)') "No method has been chosen. ", &
            "Please choose method ELPA, OMM, or PEXSI"
         stop
   end select

end subroutine


end module ELSI
