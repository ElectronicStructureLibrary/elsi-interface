!> ELSI Interchange
!! 
!! Copyright of the original code rests with the authors inside the ELSI
!! consortium. The copyright of any additional modifications shall rest
!! with their original authors, but shall adhere to the licensing terms
!! distributed along with the original code in the file "COPYING".

module ELSI

  !> ELSI Interface
  !! 

  use iso_c_binding
  use DIMENSIONS
  use MPI_TOOLS
  use HDF5_TOOLS
  use ELPA1
  use ELPA2

  implicit none
  private

  integer :: method = -1 !<Method for EV Solver (ELPA=1,OMM=2,PEXSI=3)
  integer :: mode = -1 !<Mode for EV Solver (REAL_VALUES=1,COMPLEX_VALUES=2)

  real*8, allocatable       :: H_real(:,:) !< Real Hamiltonian Matrix 
  real*8, allocatable       :: S_real(:,:) !< Real Overlap Matrix
  complex*16, allocatable   :: H_complex(:,:) !< Complex Hamiltonian Matrix
  complex*16, allocatable   :: S_complex(:,:) !< Complex Overlap Matrix
  real*8, allocatable       :: eigenvalues(:) !< Eigenvalues
  real*8, allocatable       :: vectors_real(:,:) !< Real Eigenvectors
  complex*16, allocatable   :: vectors_complex(:,:) !< Complex Eigenvectors
  
  enum, bind( C )
    enumerator :: ELPA, OMM, PEXSI
  end enum

  enum, bind( C )
    enumerator :: REAL_VALUES, COMPLEX_VALUES
  end enum   

  ! The following variables are public
  public :: ELPA, OMM, PEXSI
  public :: REAL_VALUES, COMPLEX_VALUES

  ! The following routines are public:
  public :: elsi_initialize_mpi      !< ELSI MPI initializer if not done 
                                     !! elsewhere  
  public :: elsi_initialize_problem  !< Set dimensions in code
  public :: elsi_initialize_problem_from_file !< Get dimensions from HDF5 file
  public :: elsi_set_method          !< Set routine for method
  public :: elsi_set_mode            !< Set routine for mode
  public :: elsi_allocate_matrices   !< Initialize matrices
  public :: elsi_deallocate_matrices !< Cleanup matrix memory
  public :: elsi_set_hamiltonian     !< Set Hamilitonian by Reference
  public :: elsi_set_hamiltonian_element !< Set Hamilitonian element
  public :: elsi_symmetrize_hamiltonian !< Symmetrize Hamiltonian matrix
  public :: elsi_get_hamiltonian     !< Get Hamiltonian by Reference
  public :: elsi_set_overlap         !< Set Overlap by Reference
  public :: elsi_set_overlap_element !< Set Overlap by Reference
  public :: elsi_symmetrize_overlap   !< Symmetrize overlap matrix
  public :: elsi_get_overlap         !< Get Overlap by Reference
  public :: elsi_write_ev_problem    !< Write eigenvalue problem to HDF5
  public :: elsi_read_ev_problem     !< Read eigenvalue problem from HDF5
  public :: elsi_solve_ev_problem    !< Solve eigenvalue problem 
  public :: elsi_get_eigenvalues     !< Get the eigenvalues 
  public :: elsi_finalize            !< Finalize and cleanup ELSI

  interface elsi_set_hamiltonian
     module procedure elsi_set_real_hamiltonian, &
                      elsi_set_complex_hamiltonian
  end interface 

  interface elsi_set_hamiltonian_element
     module procedure elsi_set_real_hamiltonian_element, &
                      elsi_set_complex_hamiltonian_element
  end interface 

  
  interface elsi_get_hamiltonian
     module procedure elsi_get_real_hamiltonian, &
                      elsi_get_complex_hamiltonian
  end interface 
  
  interface elsi_set_overlap
     module procedure elsi_set_real_overlap, &
                      elsi_set_complex_overlap
  end interface 

   interface elsi_set_overlap_element
     module procedure elsi_set_real_overlap_element, &
                      elsi_set_complex_overlap_element
  end interface 

  
  interface elsi_get_overlap
     module procedure elsi_set_real_overlap, &
                      elsi_set_complex_overlap
  end interface 

contains

subroutine elsi_initialize_problem(matrixsize, block_rows, block_cols)

   !>
   !!  This routine sets the method of choice for solving the eigenvalue problem
   !!

   implicit none

   integer, intent(in) :: matrixsize  !< global dimension of matrix
   integer, intent(in) :: block_rows  !< block rows of matrix
   integer, intent(in) :: block_cols  !< block cols of matrix

   n_g_rank = matrixsize
   n_b_rows = block_rows
   n_b_cols = block_cols

end subroutine


subroutine elsi_set_method(i_method)

   !>
   !!  This routine sets the method of choice for solving the eigenvalue problem
   !!

   implicit none

   integer, intent(in) :: i_method !<one of (ELPA,OMM,PEXSI)

   method = i_method
   select case (method)
      case (ELPA)
         call elsi_initialize_blacs()
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

subroutine elsi_set_mode(i_mode)

   !>
   !!  This routine sets the mode of the eigenvalue solver (Real or Complex)
   !!

   implicit none

   integer, intent(in) :: i_mode !<one of (REAL_VALUES,COMPLEX_VALUES)

   mode = i_mode


end subroutine 

subroutine elsi_allocate_matrices()

!>
!!  This routine prepares the matrices based on the dimension and method
!!  specified

   implicit none

   select case (method)
      case (ELPA)
         allocate(eigenvalues(n_g_rank))
         select case (mode)
            case (COMPLEX_VALUES)
               allocate(H_complex (n_l_rows, n_l_cols))      
               allocate(S_complex (n_l_rows, n_l_cols))
               allocate(vectors_complex(n_l_rows, n_l_cols))
            case (REAL_VALUES)
               allocate(H_real (n_l_rows, n_l_cols))      
               allocate(S_real (n_l_rows, n_l_cols))
               allocate(vectors_real(n_l_rows, n_l_cols))
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

subroutine elsi_set_real_hamiltonian_element(element,i_row,i_col)

!>
!!  This routine sets the real hamiltonian matrix 

   implicit none

   integer, intent(in) :: i_row   !<  row position
   integer, intent(in) :: i_col   !<  col position
   real*8, intent(in)  :: element !< value 
   
   select case (method)
      case (ELPA)
         if (mode == REAL_VALUES) then
            H_real(i_row,i_col) = element      
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

subroutine elsi_set_complex_hamiltonian_element(element,i_row,i_col)

!>
!!  This routine sets the real hamiltonian matrix 

   implicit none

   integer, intent(in)     :: i_row   !<  row position
   integer, intent(in)     :: i_col   !<  col position
   complex*16, intent(in)  :: element !< value 
   
   select case (method)
      case (ELPA)
         if (mode == COMPLEX_VALUES) then
            H_complex(i_row,i_col) = element      
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
   complex*16, intent(in) :: h(n_rows,n_cols) !< hamiltonian 
   
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

subroutine elsi_set_real_overlap_element(element,i_row,i_col)

!>
!!  This routine sets the real overlap element 

   implicit none

   integer, intent(in) :: i_row   !<  row position
   integer, intent(in) :: i_col   !<  col position
   real*8, intent(in)  :: element !< value 
   
   select case (method)
      case (ELPA)
         if (mode == REAL_VALUES) then
            S_real(i_row,i_col) = element      
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

subroutine elsi_symmetrize_hamiltonian()

  !>
  !!  This routine symmetrizes an upper or lower triangle hamiltonian

  implicit none

  real*8,     allocatable :: buffer_real   (:,:)
  complex*16, allocatable :: buffer_complex(:,:)
  complex*16, parameter :: CONE = (1d0,0d0)

  select case (method)
      case (ELPA)
         if (mode == REAL_VALUES) then
            allocate(buffer_real (n_l_rows, n_l_cols))
            buffer_real = H_real
            call pdtran(n_g_rank, n_g_rank, &
                  1.d0, buffer_real, 1, 1, sc_desc, &
                  1.d0, H_real, 1, 1, sc_desc)
            deallocate(buffer_real)
         else  
            allocate(buffer_complex (n_l_rows, n_l_cols))
            buffer_complex = H_complex
            call pztranc(n_g_rank, n_g_rank, &
                  CONE, buffer_complex, 1, 1, sc_desc, &
                  CONE, H_complex, 1, 1, sc_desc)
            deallocate(buffer_complex)
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

subroutine elsi_symmetrize_overlap()

   !>
  !!  This routine symmetrizes an upper or lower triangle overlap

  implicit none

  real*8,     allocatable :: buffer_real   (:,:)
  complex*16, allocatable :: buffer_complex(:,:)
  complex*16, parameter :: CONE = (1d0,0d0)

  select case (method)
      case (ELPA)
         if (mode == REAL_VALUES) then
            allocate(buffer_real (n_l_rows, n_l_cols))
            buffer_real = S_real
            call pdtran(n_g_rank, n_g_rank, &
                  1.d0, buffer_real, 1, 1, sc_desc, &
                  1.d0, S_real, 1, 1, sc_desc)
            deallocate(buffer_real)
         else  
            allocate(buffer_complex (n_l_rows, n_l_cols))
            buffer_complex = S_complex
            call pztranc(n_g_rank, n_g_rank, &
                  CONE, buffer_complex, 1, 1, sc_desc, &
                  CONE, S_complex, 1, 1, sc_desc)
            deallocate(buffer_complex)
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

subroutine elsi_set_complex_overlap_element(element,i_row,i_col)

!>
!!  This routine sets the real overlap element 

   implicit none

   integer, intent(in)     :: i_row   !<  row position
   integer, intent(in)     :: i_col   !<  col position
   complex*16, intent(in)  :: element !< value 
   
   select case (method)
      case (ELPA)
         if (mode == COMPLEX_VALUES) then
            S_complex(i_row,i_col) = element      
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


subroutine elsi_set_complex_overlap(s,n_rows,n_cols)

!>
!!  This routine sets the hamiltonian matrix 

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   complex*16, intent(in) :: s(n_rows,n_cols) !< overlap 
   
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
         if (mode == REAL_VALUES) then
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
   complex*16, intent(out) :: h(n_rows,n_cols) !< hamiltonian 
   
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
   complex*16, intent(out) :: s(n_rows,n_cols) !< overlap 
   
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

subroutine elsi_write_ev_problem(file_name)

!>
!!  This routine writes the complete eigenvalue problem into a file

   implicit none
   include "mpif.h"

   character(len=*), intent(in) :: file_name

   integer :: file_id
   integer :: group_id

   logical :: pattern

   call hdf5_initialize ()

   call hdf5_create_file (file_name, mpi_comm_global, mpi_info_null, file_id)

   ! The Hamiltonian
   call hdf5_create_group (file_id, "hamiltonian", group_id)
   
   ! Matrix dimension
   call hdf5_write_attribute (group_id, "n_matrix_rows", n_g_rank)
   call hdf5_write_attribute (group_id, "n_matrix_cols", n_g_rank)
   call hdf5_write_attribute (group_id, "n_block_rows", n_b_rows)
   call hdf5_write_attribute (group_id, "n_block_cols", n_b_cols)

   !TODO Hamiltonian Write
   call hdf5_get_scalapack_pattern()
   
   call hdf5_write_matrix_parallel (group_id, "matrix", H_real)
   
   call hdf5_close_group (group_id)

   ! The Overlap Matrix
   call hdf5_create_group (file_id, "overlap", group_id)
   
   ! Matrix dimension
   call hdf5_write_attribute (group_id, "n_matrix_rows", n_g_rank)
   call hdf5_write_attribute (group_id, "n_matrix_cols", n_g_rank)
   call hdf5_write_attribute (group_id, "n_block_rows", n_b_rows)
   call hdf5_write_attribute (group_id, "n_block_cols", n_b_cols)

   !TODO Overlap Write
   call hdf5_write_matrix_parallel (group_id, "matrix", S_real)
   call hdf5_close_group (group_id)

   call hdf5_close_file (file_id)

   call hdf5_finalize()

end subroutine

subroutine elsi_read_ev_problem(file_name)

!>
!!  This routine reads the eigenvalue problem from a file

   implicit none
   include "mpif.h"
   
   character(len=*), intent(in) :: file_name

   integer :: file_id
   integer :: group_id

   logical :: pattern

   call hdf5_initialize ()

   call hdf5_open_file (file_name, mpi_comm_global, mpi_info_null, file_id)

   ! The Hamiltonian
   call hdf5_open_group (file_id, "hamiltonian", group_id)
   
   ! Matrix dimension
   call hdf5_read_attribute (group_id, "n_matrix_rows", n_g_rank)
   call hdf5_read_attribute (group_id, "n_matrix_cols", n_g_rank)
   call hdf5_read_attribute (group_id, "n_block_rows", n_b_rows)
   call hdf5_read_attribute (group_id, "n_block_cols", n_b_cols)

   ! Hamiltonian Read
   call hdf5_get_scalapack_pattern()
   call hdf5_read_matrix_parallel (group_id, "matrix", H_real)
   
   call hdf5_close_group (group_id)

   ! The Overlap Matrix
   call hdf5_open_group (file_id, "overlap", group_id)
   
   ! Matrix dimension
   call hdf5_read_attribute (group_id, "n_matrix_rows", n_g_rank)
   call hdf5_read_attribute (group_id, "n_matrix_cols", n_g_rank)
   call hdf5_read_attribute (group_id, "n_block_rows", n_b_rows)
   call hdf5_read_attribute (group_id, "n_block_cols", n_b_cols)

   ! Overlap Read
   call hdf5_read_matrix_parallel (group_id, "matrix", S_real)
   
   call hdf5_close_group (group_id)

   call hdf5_close_file (file_id)

   call hdf5_finalize()

end subroutine

subroutine elsi_initialize_problem_from_file(file_name)

   !>
   !!  This routine sets the method of choice for solving the eigenvalue problem
   !!

   implicit none
   include "mpif.h"

   character(len=*), intent(in) :: file_name

   integer :: file_id
   integer :: group_id

   call hdf5_initialize ()

   call hdf5_open_file (file_name, mpi_comm_global, mpi_info_null, file_id)

   ! The Hamiltonian
   call hdf5_open_group (file_id, "hamiltonian", group_id)

   ! Matrix dimension
   call hdf5_read_attribute (group_id, "n_matrix_rows", n_g_rank)
   call hdf5_read_attribute (group_id, "n_matrix_cols", n_g_rank)
   call hdf5_read_attribute (group_id, "n_block_rows",  n_b_rows)
   call hdf5_read_attribute (group_id, "n_block_cols",  n_b_cols)

   call hdf5_close_group (group_id)

   call hdf5_close_file (file_id)

   call hdf5_finalize()


end subroutine

subroutine elsi_solve_ev_problem(n_vectors)

   !>
   !!  This routine interfaces to the possible eigenvalue solvers

   implicit none
   include "mpif.h"

   integer, intent(in) :: n_vectors !< Number of eigenvectors to be calculated

   logical :: success

   select case (method)
      case (ELPA)
         ! TODO if S is not 1 we need the cholesky decomposition first
         !      also we need to choose if 1stage or 2 stage elpa
         select case (mode)
            case (COMPLEX_VALUES)
               success = solve_evp_complex_2stage( &
                     n_g_rank, n_vectors, H_complex, &
                     n_l_rows, eigenvalues, vectors_complex, &
                     n_l_rows, n_b_rows, &
                     mpi_comm_row, mpi_comm_col, mpi_comm_global)
            case (REAL_VALUES)
               success = solve_evp_real_2stage( &
                     n_g_rank, n_vectors, H_real, &
                     n_l_rows, eigenvalues, vectors_real, &
                     n_l_rows, n_b_rows,&
                     mpi_comm_row, mpi_comm_col, mpi_comm_global)
         end select
       
         if (.not.success) then
            write(*,'(a)') "ELPA failed."
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

   call MPI_BARRIER(mpi_comm_global, mpierr)

end subroutine

subroutine elsi_get_eigenvalues(eigenvalues_out,n_eigenvalues)

!>
!!  This routine gets the eigenvalues 

   implicit none

   integer, intent(in) :: n_eigenvalues !< Number of eigenvalues
   real*8, intent(out) :: eigenvalues_out(n_eigenvalues) !< eigenvalues 
   
   select case (method)
      case (ELPA)
            eigenvalues_out(:) = eigenvalues(:)      
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



subroutine elsi_deallocate_matrices()

!>
!!  This routine prepares the matrices based on the dimension and method
!!  specified

   implicit none

   select case (method)
      case (ELPA)
         deallocate(eigenvalues) 
         select case (mode)
            case (COMPLEX_VALUES)
               deallocate(H_complex)      
               deallocate(S_complex)
               deallocate(vectors_complex)
            case (REAL_VALUES)
               deallocate(H_real)      
               deallocate(S_real)
               deallocate(vectors_real)
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

subroutine elsi_finalize()

   !>
   !!  This routine prepares the matrices based on the dimension and method
   !!  specified

   implicit none
   include "mpif.h"

   call MPI_BARRIER(mpi_comm_global,mpierr)

   call elsi_deallocate_matrices()

   call elsi_finalize_blacs()

   call elsi_finalize_mpi()

end subroutine

end module ELSI
