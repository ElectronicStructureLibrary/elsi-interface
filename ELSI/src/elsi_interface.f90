!Copyright (c) 2015, ELSI consortium
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
! * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer in the
!   documentation and/or other materials provided with the distribution.
! * Neither the name of the ELSI project nor the
!   names of its contributors may be used to endorse or promote products
!   derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!> ELSI Interchange
!! This module is the actual ELSI Interface, providing functions for setting up
!! and solving or circumventing an eigenvalue problem using ELPA, OMM, or PEXSI
!!

module ELSI

  use iso_c_binding
  use ELSI_DIMENSIONS
  use ELSI_TIMERS
  use ELSI_MPI_TOOLS
  use ELSI_HDF5_TOOLS
  use ELSI_MATRIX_CONVERSION
  use ELPA1
  use ELPA2
  use matrixswitch
  use f_ppexsi_interface

  implicit none
  private

  !< Internal Storage
  !< Real Hamiltonian Matrix 
  real*8, target, allocatable     :: H_real_target(:,:) 
  !< Complex Hamiltonian Matrix 
  complex*16, target, allocatable :: H_complex_target(:,:) 
  !< Real Overlap Matrix
  real*8, target, allocatable     :: S_real_target(:,:)   
  !< Complex Overlap Matrix
  complex*16, target, allocatable :: S_complex_target(:,:)   
  !< Eigenvalues   
  real*8, allocatable             :: eigenvalues(:)
  !< Real Eigenvectors   
  real*8, allocatable             :: vectors_real(:,:)
  !< Complex Eigenvectors 
  complex*16, allocatable         :: vectors_complex(:,:) 

  ! Pointers
  !< Real Hamiltonian Matrix 
  real*8, pointer                 :: H_real(:,:)
  !< External Complex Hamiltonian Matrix
  complex*16, pointer             :: H_complex(:,:) 
  !< External Real Overlap Matrix
  real*8, pointer                 :: S_real(:,:)   
  !< External Complex Overlap Matrix
  complex*16, pointer             :: S_complex(:,:)
  
  !OMM
  !< OMM Hamiltonian matrix
  type(Matrix) :: OMM_H_matrix
  !< OMM Overlap matrix
  type(Matrix) :: OMM_S_matrix
  !< OMM coefficient matrix
  type(Matrix) :: OMM_C_matrix
  !< OMM density matrix
  type(Matrix) :: OMM_D_matrix
  !< OMM kinetic energy density matrix
  type(Matrix) :: OMM_T_matrix

  !PESXI
  !< Sparse Hamiltonian
  real*8,  allocatable  :: H_real_sparse(:)
  !< Sparse Overlap 
  real*8,  allocatable  :: S_real_sparse(:)
  !< Sparse Density matrix 
  real*8,  allocatable  :: D_real_sparse(:)
  !< Sparse Energydensity matrix 
  real*8,  allocatable  :: ED_real_sparse(:)
  !< Sparse Freeenergydensity matrix 
  real*8,  allocatable  :: FD_real_sparse(:)
  !< Sparse index array
  integer, allocatable  :: sparse_index(:)
  !< Sparse pointer array
  integer, allocatable  :: sparse_pointer(:)
  

  !< The following variables from ELSI Dimensions are public
  public :: ELPA, OMM_DENSE, PEXSI
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
  public :: elsi_get_total_energy    !< Get the sum of occupied eigenvalues 
  public :: elsi_finalize            !< Finalize and cleanup ELSI

  ! from other modules
  public :: elsi_initialize_blacs
  public :: elsi_get_local_dimensions
  public :: elsi_get_myid
  public :: elsi_get_global_row
  public :: elsi_get_global_col
  public :: elsi_set_mpi
  public :: elsi_set_blacs
  public :: elsi_get_global_dimensions

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

!>
!!  This routine sets the matrix dimensions for the eigenvalue problem
!!
subroutine elsi_initialize_problem(matrixsize, block_rows, block_cols)


   implicit none

   integer, intent(in) :: matrixsize  !< global dimension of matrix
   integer, intent(in) :: block_rows  !< block rows of matrix
   integer, intent(in) :: block_cols  !< block cols of matrix

   n_g_rank = matrixsize
   n_b_rows = block_rows
   n_b_cols = block_cols

end subroutine


!>
!!  This routine sets the method of choice for solving the eigenvalue problem
!!
subroutine elsi_set_method(i_method)


   implicit none

   integer, intent(in) :: i_method !<one of (ELPA,OMM,PEXSI)
   
   call elsi_initialize_timers()
   call elsi_start_total_time()
   call hdf5_initialize ()

   method = i_method
   select case (method)
      case (ELPA)
      case (OMM_DENSE)
      case (PEXSI)
      case DEFAULT
         call elsi_stop("No method has been chosen."//&
               " Please choose method ELPA, OMM_DENSE, or PEXSI", "elsi_set_method")
   end select

end subroutine

!>
!!  This routine sets the mode of the eigenvalue solver (Real or Complex)
!!
subroutine elsi_set_mode(i_mode)


   implicit none

   integer, intent(in) :: i_mode !<one of (REAL_VALUES,COMPLEX_VALUES)

   mode = i_mode

   if (i_mode == COMPLEX_VALUES) then
     call elsi_stop("COMPLEX VALUES not yet supported", "elsi_set_mode")
   end if

end subroutine 

!>
!!  This routine prepares the matrices based on the dimension and method
!!  specified
subroutine elsi_allocate_matrices()


   implicit none

   character*200 :: message

   select case (method)
      case (ELPA)
         allocate(eigenvalues(n_g_rank))
         select case (mode)
            case (COMPLEX_VALUES)
               allocate(H_complex_target (n_l_rows, n_l_cols))      
               allocate(S_complex_target (n_l_rows, n_l_cols))
               allocate(vectors_complex(n_l_rows, n_l_cols))
               H_complex => H_complex_target
               S_complex => S_complex_target
               H_complex = CMPLX(0d0,0d0)
               S_complex = CMPLX(0d0,0d0)
            case (REAL_VALUES)
               allocate(H_real_target (n_l_rows, n_l_cols))      
               allocate(S_real_target (n_l_rows, n_l_cols))
               allocate(vectors_real(n_l_rows, n_l_cols))
               H_real => H_real_target
               S_real => S_real_target
               H_real = 0d0
               S_real = 0d0
            case DEFAULT
               call elsi_stop("No mode has been chosen. "// &
                  "Please choose method REAL_VALUES or COMPLEX_VALUES", &
                  "elsi_allocate_matrices")
         end select
      case (OMM_DENSE)
         select case (mode)
            case (COMPLEX_VALUES)
               allocate(H_complex_target (n_l_rows, n_l_cols))      
               allocate(S_complex_target (n_l_rows, n_l_cols))
               allocate(vectors_complex(n_l_rows, n_l_cols))
               H_complex => H_complex_target
               S_complex => S_complex_target
               H_complex = CMPLX(0d0,0d0)
               S_complex = CMPLX(0d0,0d0)
               call m_register_pdbc (OMM_H_matrix, H_complex, sc_desc)
               call m_register_pdbc (OMM_S_matrix, S_complex, sc_desc)
               call m_allocate (OMM_D_matrix, n_g_rank, n_g_rank, "pzdbc")
               call m_allocate (OMM_T_matrix, n_g_rank, n_g_rank, "pzdbc")
               OMM_T_matrix%zval = CMPLX(0d0,0d0)
               OMM_D_matrix%zval = CMPLX(0d0,0d0)

            case (REAL_VALUES)
               allocate(H_real_target (n_l_rows, n_l_cols))      
               allocate(S_real_target (n_l_rows, n_l_cols))
               allocate(vectors_real(n_l_rows, n_l_cols))
               H_real => H_real_target
               S_real => S_real_target
               H_real = 0d0
               S_real = 0d0

               call m_register_pdbc (OMM_H_matrix, H_real, sc_desc)
               call m_register_pdbc (OMM_S_matrix, S_real, sc_desc)
               call m_allocate (OMM_D_matrix, n_g_rank, n_g_rank, "pddbc")
               call m_allocate (OMM_T_matrix, n_g_rank, n_g_rank, "pddbc")
               OMM_T_matrix%dval = 0d0
               OMM_D_matrix%dval = 0d0

            case DEFAULT
               call elsi_stop("No mode has been chosen. "// &
                  "Please choose method REAL_VALUES or COMPLEX_VALUES", &
                  "elsi_allocate_matrices")
         end select

      case (PEXSI)
         select case (mode)
            case (REAL_VALUES)
               ! We do not allocate anything here, yet.
            case DEFAULT
               call elsi_stop("No mode has been chosen. "// &
                  "Please choose method REAL_VALUES or COMPLEX_VALUES", &
                  "elsi_allocate_matrices")
         end select

      case DEFAULT
         call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_allocate_matrices")
   end select
end subroutine

!>
!!  This routine sets the element (i_row, i_col) real hamiltonian matrix 
!!
subroutine elsi_set_real_hamiltonian_element(element,i_row,i_col)


   implicit none

   integer, intent(in) :: i_row   !<  row position
   integer, intent(in) :: i_col   !<  col position
   real*8, intent(in)  :: element !< value 
   
   select case (method)
      case (ELPA,OMM_DENSE)
         if (mode == REAL_VALUES) then
            if(associated(H_real)) then
               H_real(i_row,i_col) = element
            else 
              call elsi_stop("Hamiltonian not created/linked.",&
                    "elsi_set_real_hamiltonian_element") 
            endif
         else  
            call elsi_stop("Wrong mode: "// &
               "Complex valued hamiltonian to be written in real storage",&
               "elsi_set_real_hamiltonian_element")
         end if
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",&
               "elsi_set_real_hamiltonian_element")
      case DEFAULT
         call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_set_real_hamiltonian_element")
   end select


end subroutine


!>
!!  This routine sets the full local real hamiltonian matrix 
!!
subroutine elsi_set_real_hamiltonian(h,n_rows,n_cols)


   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   real*8, target, intent(in) :: h(n_rows,n_cols) !< hamiltonian 
   
   select case (method)
      case (ELPA)
         if (mode == REAL_VALUES) then
            H_real => h      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Complex valued hamiltonian to be written in real storage"
            stop
         end if
      case (OMM_DENSE)
         call m_register_pdbc(OMM_H_matrix,h,sc_desc)
         H_real => OMM_H_matrix%dval
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",&
               "elsi_set_real_hamiltonian")
      case DEFAULT
          call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_set_real_hamiltonian")
   end select


end subroutine

!>
!!  This routine sets the element (i_row, i_col) in the complex 
!!  hamiltonian matrix 
!!
subroutine elsi_set_complex_hamiltonian_element(element,i_row,i_col)


   implicit none

   integer, intent(in)     :: i_row   !<  row position
   integer, intent(in)     :: i_col   !<  col position
   complex*16, intent(in)  :: element !< value 
   
   select case (method)
      case (ELPA,OMM_DENSE)
         if (mode == COMPLEX_VALUES) then
            H_complex(i_row,i_col) = element      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Complex valued hamiltonian to be written in real storage"
            stop
         end if
      case (PEXSI)
        call elsi_stop("PEXSI not implemented yet!",&
               "elsi_set_complex_hamiltonian_element")
      case DEFAULT
          call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_set_complex_hamiltonian_element")
   end select


end subroutine

!>
!!  This routine sets the full complex hamiltonian matrix 
!!
subroutine elsi_set_complex_hamiltonian(h,n_rows,n_cols)

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   complex*16, target, intent(in) :: h(n_rows,n_cols) !< hamiltonian 
   
   select case (method)
      case (ELPA)
         if (mode == COMPLEX_VALUES) then
            H_complex => h      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Real valued hamiltonian to be written in complex storage"
            stop
         end if
      case (OMM_DENSE)
         call m_register_pdbc(OMM_H_matrix,h,sc_desc)
         H_complex => OMM_H_matrix%zval
         write(*,'(a)') "OMM not implemented yet!"
         stop
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",&
               "elsi_set_complex_hamiltonian")
      case DEFAULT
          call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_set_complex_hamiltonian")
   end select


end subroutine

!>
!!  This routine sets the element (i_row, i_col) in the real 
!!  overlap matrix 
!!
subroutine elsi_set_real_overlap_element(element,i_row,i_col)

   implicit none

   integer, intent(in) :: i_row   !<  row position
   integer, intent(in) :: i_col   !<  col position
   real*8, intent(in)  :: element !< value 
   
   select case (method)
      case (ELPA,OMM_DENSE)
         if (mode == REAL_VALUES) then
            S_real(i_row,i_col) = element      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Complex valued overlap to be written in real storage"
            stop
         end if
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",&
               "elsi_set_real_overlap")
      case DEFAULT
         call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_set_real_overlap_element")
   end select

   overlap_is_unity = .False.


end subroutine

!>
!!  This routine symmetrizes an upper or lower triangle hamiltonian
!!
subroutine elsi_symmetrize_hamiltonian()


  implicit none

  real*8,     allocatable :: buffer_real   (:,:)
  complex*16, allocatable :: buffer_complex(:,:)
  complex*16, parameter   :: CONE = (1d0,0d0)

  integer :: l_row, l_col  !< local matrix indices
  integer :: g_row, g_col  !< global matrix indices

  select case (method)
      case (ELPA,OMM_DENSE)
         if (mode == REAL_VALUES) then
            allocate(buffer_real (n_l_rows, n_l_cols))
            buffer_real(:,:) = H_real(:,:)
            call pdtran(n_g_rank, n_g_rank, &
                  1.d0, buffer_real, 1, 1, sc_desc, &
                  1.d0, H_real, 1, 1, sc_desc)
            deallocate(buffer_real)
            
            do l_row = 1, n_l_rows
               call elsi_get_global_row(g_row, l_row)
               do l_col = l_row, n_l_cols
                  call elsi_get_global_col(g_col, l_col)
                  if (g_row == g_col) then
                     H_real(l_row,l_col) = 0.5d0 * H_real(l_row,l_col)
                  end if
               end do
            end do

         else  
            allocate(buffer_complex (n_l_rows, n_l_cols))
            buffer_complex(:,:) = H_complex(:,:)
            call pztranc(n_g_rank, n_g_rank, &
                  CONE, buffer_complex, 1, 1, sc_desc, &
                  CONE, H_complex, 1, 1, sc_desc)
            deallocate(buffer_complex)

            do l_row = 1, n_l_rows
               call elsi_get_global_row(g_row, l_row)
               do l_col = l_row, n_l_cols
                  call elsi_get_global_col(g_col, l_col)
                  if (g_row == g_col) then
                     H_complex(l_row,l_col) = 0.5d0 * H_complex(l_row,l_col)
                  end if
               end do
            end do

         end if
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",&
               "elsi_symmetrize_hamiltonian")
      case DEFAULT
         call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_symetrize_hamiltonian")
   end select

   if (allocated(buffer_real))    deallocate (buffer_real)
   if (allocated(buffer_complex)) deallocate (buffer_complex)

end subroutine 

!>
!!  This routine symmetrizes an upper or lower triangle overlap
!!
subroutine elsi_symmetrize_overlap()

   !>
  !!  This routine symmetrizes an upper or lower triangle overlap

  implicit none

  real*8,     allocatable :: buffer_real   (:,:)
  complex*16, allocatable :: buffer_complex(:,:)
  complex*16, parameter   :: CONE = (1d0,0d0)

  integer :: l_row, l_col  !< local matrix indices
  integer :: g_row, g_col  !< global matrix indices
  
  select case (method)
      case (ELPA,OMM_DENSE)
         if (mode == REAL_VALUES) then
            allocate(buffer_real (n_l_rows, n_l_cols))
            buffer_real(:,:) = S_real(:,:)
            call pdtran(n_g_rank, n_g_rank, &
                  1.d0, buffer_real, 1, 1, sc_desc, &
                  1.d0, S_real, 1, 1, sc_desc)
            deallocate(buffer_real)

            do l_row = 1, n_l_rows
               call elsi_get_global_row(g_row, l_row)
               do l_col = l_row, n_l_cols
                  call elsi_get_global_col(g_col, l_col)
                  if (g_row == g_col) then
                     S_real(l_row,l_col) = 0.5d0 * S_real(l_row,l_col)
                  end if
               end do
            end do

         else  
            allocate(buffer_complex (n_l_rows, n_l_cols))
            buffer_complex(:,:) = S_complex(:,:)
            call pztranc(n_g_rank, n_g_rank, &
                  CONE, buffer_complex, 1, 1, sc_desc, &
                  CONE, S_complex, 1, 1, sc_desc)
            deallocate(buffer_complex)

            do l_row = 1, n_l_rows
               call elsi_get_global_row(g_row, l_row)
               do l_col = l_row, n_l_cols
                  call elsi_get_global_col(g_col, l_col)
                  if (g_row == g_col) then
                     S_complex(l_row,l_col) = 0.5d0 * S_complex(l_row,l_col)
                  end if
               end do
            end do

         end if
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",&
               "elsi_symmetrize_overlap")
      case DEFAULT
         call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_symmetrize_overlap")
   end select
   
   if (allocated(buffer_real))    deallocate (buffer_real)
   if (allocated(buffer_complex)) deallocate (buffer_complex)

end subroutine 



!>
!!  This routine sets the real overlap matrix
!!
subroutine elsi_set_real_overlap(s,n_rows,n_cols)

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   real*8, target, intent(in) :: s(n_rows,n_cols) !< overlap 
   
   select case (method)
      case (ELPA)
         if (mode == REAL_VALUES) then
            S_real => s      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Complex valued overlap to be written in real storage"
            stop
         end if
      case (OMM_DENSE)
         call m_register_pdbc(OMM_S_matrix,s,sc_desc)
         S_real => OMM_S_matrix%dval
      case (PEXSI)
        call elsi_stop("PEXSI not implemented yet!",&
               "elsi_set_real_overlap")
      case DEFAULT
         call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_set_real_overlap")
   end select

   overlap_is_unity = .False.

end subroutine

!>
!!  This routine sets one element (i_row, i_col) of the complex overlap matrix
!!
subroutine elsi_set_complex_overlap_element(element,i_row,i_col)

   implicit none

   integer, intent(in)     :: i_row   !<  row position
   integer, intent(in)     :: i_col   !<  col position
   complex*16, intent(in)  :: element !< value 
   
   select case (method)
      case (ELPA,OMM_DENSE)
         if (mode == COMPLEX_VALUES) then
            S_complex(i_row,i_col) = element      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Real valued overlap to be written in complex storage"
            stop
         end if
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",&
               "elsi_set_complex_overlap_element")
      case DEFAULT
         call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_set_complex_overlap_element")
   end select

   overlap_is_unity = .False.


end subroutine


!>
!!  This routine sets the full complex overlap matrix
!!
subroutine elsi_set_complex_overlap(s,n_rows,n_cols)

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   complex*16, target, intent(in) :: s(n_rows,n_cols) !< overlap 
   
   select case (method)
      case (ELPA)
         if (mode == REAL_VALUES) then
            S_complex => s      
         else  
            write(*,'(2a)') "Wrong mode:", &
               " Real valued overlap to be written in complex storage"
            stop
         end if
      case (OMM_DENSE)
         call m_register_pdbc(OMM_S_matrix,s,sc_desc)
         S_complex => OMM_S_matrix%zval
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",&
               "elsi_set_complex_overlap")
      case DEFAULT
         call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_set_complex_overlap")
   end select

   overlap_is_unity = .False.

end subroutine

!>
!!  This routine gets the real hamiltonian matrix
!!
subroutine elsi_get_real_hamiltonian(h,n_rows,n_cols)


   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   real*8, intent(out) :: h(n_rows,n_cols) !< hamiltonian 
   
   select case (method)
      case (ELPA,OMM_DENSE)
         if (mode == REAL_VALUES) then
            h(:,:) = H_real(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Real valued hamiltonian to be written in complex storage"
            stop
         end if
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",&
               "elsi_get_real_hamiltonian")
      case DEFAULT
         call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_get_real_hamiltonian")
   end select

end subroutine

!>
!!  This routine gets the complex hamiltonian matrix 
!!
subroutine elsi_get_complex_hamiltonian(h,n_rows,n_cols)

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   complex*16, intent(out) :: h(n_rows,n_cols) !< hamiltonian 
   
   select case (method)
      case (ELPA,OMM_DENSE)
         if (mode == COMPLEX_VALUES) then
            h(:,:) = H_complex(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Complex valued hamiltonian to be written in real storage"
            stop
         end if
      case (PEXSI)
        call elsi_stop("PEXSI not implemented yet!",&
               "elsi_get_complex_hamiltonian")
      case DEFAULT
         call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_get_complex_hamiltonian")
   end select

end subroutine

!>
!!  This routine gets the real overlap matrix
!!
subroutine elsi_get_real_overlap(s,n_rows,n_cols)

   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   real*8, intent(out) :: s(n_rows,n_cols) !< overlap 
   
   select case (method)
      case (ELPA,OMM_DENSE)
         if (mode == REAL_VALUES) then
            s(:,:) = S_real(:,:)    
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Real valued overlap to be written in complex storage"
            stop
         end if
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",&
               "elsi_get_real_overlap")
      case DEFAULT
         call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_get_real_overlap")
   end select

end subroutine

!>
!!  This routine gets the complex overlap matrix 
!!
subroutine elsi_get_complex_overlap(s,n_rows,n_cols)


   implicit none

   integer, intent(in) :: n_rows !< Number of rows
   integer, intent(in) :: n_cols !< Number of cols
   complex*16, intent(out) :: s(n_rows,n_cols) !< overlap 
   
   select case (method)
      case (ELPA,OMM_DENSE)
         if (mode == COMPLEX_VALUES) then
            s(:,:) = S_complex(:,:)      
         else  
            write(*,'(2a)') "Wrong mode:", &
               "Complex valued overlap to be written in real storage"
            stop
         end if
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!",&
               "elsi_get_complex_overlap")
      case DEFAULT
         write(*,'(2a)') "No method has been chosen. ", &
            "Please choose method ELPA, OMM_DENSE, or PEXSI"
         stop
   end select

end subroutine

!>
!!  This routine writes the complete eigenvalue problem into a file
!!
subroutine elsi_write_ev_problem(file_name)


   implicit none
   include "mpif.h"

   character(len=*), intent(in) :: file_name !< File name where to write

   integer :: file_id   !< HDF5 File identifier
   integer :: group_id  !< HDF5 Group identifier
   real*8, allocatable :: buffer(:,:) !< Read buffer for PEXSI

   call elsi_start_write_evp_time()

   if (method == PEXSI) then
      allocate(buffer(n_l_rows, n_l_cols)) 
   end if



   call hdf5_create_file (file_name, mpi_comm_global, mpi_info_null, file_id)

   ! The Hamiltonian
   call hdf5_create_group (file_id, "hamiltonian", group_id)
   
   ! Matrix dimension
   call hdf5_write_attribute (group_id, "n_matrix_rows", n_g_rank)
   call hdf5_write_attribute (group_id, "n_matrix_cols", n_g_rank)

   ! Hamiltonian Write
   call hdf5_get_scalapack_pattern()
   
   select case (method)
      case (ELPA,OMM_DENSE) 
         call hdf5_write_matrix_parallel (group_id, "matrix", H_real)
      case (PEXSI)
         call elsi_ccs_to_dense(buffer, H_real_sparse, sparse_index,&
               sparse_pointer)
         call hdf5_write_matrix_parallel (group_id, "matrix", buffer)
      case DEFAULT
         call elsi_stop("No supported method has been chosen. "&
            // "Please choose method ELPA, OMM_DENSE, or PEXSI",&
               "elsi_write_ev_problem")
   end select

   call hdf5_close_group (group_id)

   ! The Overlap Matrix
   call hdf5_create_group (file_id, "overlap", group_id)
   
   ! Matrix dimension
   call hdf5_write_attribute (group_id, "n_matrix_rows", n_g_rank)
   call hdf5_write_attribute (group_id, "n_matrix_cols", n_g_rank)

   ! Overlap Write
   select case (method)
      case (ELPA,OMM_DENSE) 
         call hdf5_write_matrix_parallel (group_id, "matrix", S_real)
      case (PEXSI)
         call elsi_ccs_to_dense(buffer, S_real_sparse, sparse_index,&
               sparse_pointer)
         call hdf5_write_matrix_parallel (group_id, "matrix", buffer)
      case DEFAULT
         call elsi_stop("No supported method has been chosen. "&
            // "Please choose method ELPA, OMM_DENSE, or PEXSI",&
               "elsi_write_ev_problem")
   end select

   call hdf5_close_group (group_id)

   call hdf5_close_file (file_id)

   call elsi_stop_write_evp_time()

end subroutine

!>
!!  This routine reads the eigenvalue problem from a file
!!
subroutine elsi_read_ev_problem(file_name)


   implicit none
   include "mpif.h"
   
   character(len=*), intent(in) :: file_name !< File to open

   integer :: file_id  !< HDF5 File identifier
   integer :: group_id !< HDF5 Group identifier
   real*8, allocatable :: buffer(:,:) !< Read buffer for PEXSI

   ! For Pexsi we need to create a buffer 
   ! we convert it directly to the CCS format

   call elsi_start_read_evp_time()

   if (method == PEXSI) then
      allocate(buffer(n_l_rows, n_l_cols)) 
   end if

   call hdf5_open_file (file_name, mpi_comm_global, mpi_info_null, file_id)

   ! The Hamiltonian
   call hdf5_open_group (file_id, "hamiltonian", group_id)
   
   ! Hamiltonian Read
   call hdf5_get_scalapack_pattern()
   select case (method)
      case (ELPA,OMM_DENSE)
         call hdf5_read_matrix_parallel (group_id, "matrix", H_real)
      case (PEXSI)
         call hdf5_read_matrix_parallel (group_id, "matrix", buffer)
         call elsi_compute_N_nonzero(buffer)
         allocate(H_real_sparse(n_l_nonzero))
         allocate(sparse_index(n_l_nonzero))
         allocate(sparse_pointer(n_l_cols + 1))
         call elsi_dense_to_ccs(buffer, H_real_sparse, sparse_index,&
               sparse_pointer)
      case DEFAULT
         call elsi_stop("No supported method has been chosen. "&
            // "Please choose method ELPA, OMM_DENSE, or PEXSI",&
               "elsi_read_ev_problem")
    end select

   
   call hdf5_close_group (group_id)

   ! The Overlap Matrix
   call hdf5_open_group (file_id, "overlap", group_id)
   
   ! Overlap Read
   select case (method)
      case (ELPA,OMM_DENSE)
         call hdf5_read_matrix_parallel (group_id, "matrix", S_real)
      case (PEXSI)
         call hdf5_read_matrix_parallel (group_id, "matrix", buffer)
         allocate(S_real_sparse(n_l_nonzero))
         call elsi_dense_to_ccs_by_pattern(buffer, S_real_sparse,&
               sparse_index, sparse_pointer)
      case DEFAULT
         call elsi_stop("No supported method has been chosen. "&
            // "Please choose method ELPA, OMM_DENSE, or PEXSI",&
               "elsi_read_ev_problem")
         stop
   end select

   ! TODO Check if overlap is unity
   overlap_is_unity = .False.
   
   call hdf5_close_group (group_id)

   call hdf5_close_file (file_id)

   if (allocated(buffer)) deallocate (buffer)

   call elsi_stop_read_evp_time()

end subroutine

!>
!!  This routine sets the method of choice for solving the eigenvalue problem
!!
subroutine elsi_initialize_problem_from_file(file_name, block_rows, block_cols)


   implicit none
   include "mpif.h"

   character(len=*), intent(in) :: file_name !< File to open
   integer, intent(in) :: block_rows  !< block rows of matrix
   integer, intent(in) :: block_cols  !< block cols of matrix

   integer :: file_id  !< HDF5 File identifier 
   integer :: group_id !< HDF5 Group identifier

   call elsi_start_read_evp_time()

   n_b_rows = block_rows
   n_b_cols = block_cols

   if (.not. mpi_is_setup) call elsi_stop("MPI needs to be setup first!", &
         "elsi_initialize_problem_from_file")

   call hdf5_open_file (file_name, mpi_comm_global, mpi_info_null, file_id)

   ! The Hamiltonian
   call hdf5_open_group (file_id, "hamiltonian", group_id)

   ! Matrix dimension
   call hdf5_read_attribute (group_id, "n_matrix_rows", n_g_rank)
   call hdf5_read_attribute (group_id, "n_matrix_cols", n_g_rank)

   call hdf5_close_group (group_id)

   call hdf5_close_file (file_id)

   
   call elsi_stop_read_evp_time()

end subroutine

!>
!!  This routine interfaces to the eigenvalue solvers
!!
subroutine elsi_solve_ev_problem(number_of_electrons)


   implicit none
   include "mpif.h"

   !< Number of electrons in the system
   real*8, intent(in) :: number_of_electrons 

   logical :: success
   logical :: two_step_solver
   real*8  :: val
   integer :: i,j
   integer :: i_task
   character(len=4096) :: string_message

   call elsi_start_solve_evp_time()

   if (method == ELPA .or. method == OMM_DENSE) then
      if (mode == REAL_VALUES .and. .not. associated(H_real)) then
         call elsi_stop("Hamiltonian not created/linked.",&
                       "elsi_solve_ev_problem")
      end if
   
      if (mode == COMPLEX_VALUES .and. .not. associated(H_complex)) then
         call elsi_stop("Hamiltonian not created/linked.",&
                       "elsi_solve_ev_problem")
      end if
   end if

   if (method == PEXSI) then
      if (mode == REAL_VALUES .and. .not. allocated(H_real_sparse)) then
         call elsi_stop("Hamiltonian not created/linked.",&
                       "elsi_solve_ev_problem")
      end if
   end if

   n_electrons = number_of_electrons
   n_eigenvectors = ceiling(n_electrons/2d0)

   ! Debug
   !call elsi_print_setup()
   !call elsi_variable_status()

   select case (method)
      case (ELPA)
         two_step_solver = .True.
         if (.not. overlap_is_unity) then
            call elsi_to_standard_eigenvalue_problem ()
         end if

         if (1d0 * n_eigenvectors/n_g_rank > elpa_step_switch) then
           two_step_solver = .False.
         end if


         if (two_step_solver) then
           select case (mode)
             case (COMPLEX_VALUES)
               success = solve_evp_complex_2stage( &
                     n_g_rank, n_eigenvectors, H_complex, &
                     n_l_rows, eigenvalues, vectors_complex, &
                     n_l_rows, n_b_rows, n_l_cols,&
                     mpi_comm_row, mpi_comm_col, mpi_comm_global)
             case (REAL_VALUES)
               success = solve_evp_real_2stage( &
                     n_g_rank, n_eigenvectors, H_real, &
                     n_l_rows, eigenvalues, vectors_real, &
                     n_l_rows, n_b_rows, n_l_cols,&
                     mpi_comm_row, mpi_comm_col, mpi_comm_global)
           end select
         else
           select case (mode)
             case (COMPLEX_VALUES)
               success = solve_evp_complex( &
                     n_g_rank, n_eigenvectors, H_complex, &
                     n_l_rows, eigenvalues, vectors_complex, &
                     n_l_rows, n_b_rows, n_l_cols,&
                     mpi_comm_row, mpi_comm_col)
             case (REAL_VALUES)
               success = solve_evp_real( &
                     n_g_rank, n_eigenvectors, H_real, &
                     n_l_rows, eigenvalues, vectors_real, &
                     n_l_rows, n_b_rows, n_l_cols,&
                     mpi_comm_row, mpi_comm_col)
           end select

         end if
       
         if (.not.success) then
            call elsi_stop("ELPA failed in solving eigenvalue problem.",&
                  "elsi_solve_ev_problem")
         end if

      case (OMM_DENSE)

        call elsi_set_omm_default_options()

        ! Prepare Coefficient matrix
        select case (mode)
            case (COMPLEX_VALUES)
               call m_allocate (OMM_C_matrix, n_eigenvectors, n_g_rank, "pzdbc")
               C_matrix_initialized = .False.
            case (REAL_VALUES)
               call m_allocate (OMM_C_matrix, n_eigenvectors, n_g_rank, "pddbc")
               C_matrix_initialized = .False.
         end select

        ! Shift eigenvalue spectrum
        call m_add(OMM_S_matrix,'N',OMM_H_matrix,-eta,1d0,"lap")

        call omm(n_g_rank, n_eigenvectors, OMM_H_matrix, OMM_S_matrix, &
              new_overlap, total_energy, OMM_D_matrix, calc_ED, eta, &
              OMM_C_matrix, C_matrix_initialized, OMM_T_matrix, &
              scale_kinetic, omm_flavour, nk_times_nspin, i_k_spin,&
              min_tol, omm_verbose, do_dealloc, "pddbc", "lap", myid+1)

      case (PEXSI)

         if (.not. allocated(D_real_sparse))  &
            allocate (D_real_sparse(n_l_nonzero))
         if (.not. allocated(ED_real_sparse)) &
            allocate (ED_real_sparse(n_l_nonzero))
         if (.not. allocated(FD_real_sparse)) &
            allocate (FD_real_sparse(n_l_nonzero))

         ! Set the default options
         !TODO: User interface is missing
         call elsi_set_pexsi_default_options()

         ! Load sparse matrices for PEXSI
         if (overlap_is_unity) then
           call f_ppexsi_load_real_symmetric_hs_matrix( pexsi_plan,&
               pexsi_options, n_g_rank, n_g_nonzero, n_l_nonzero, n_l_cols,&
               sparse_pointer, sparse_index, H_real_sparse, 1, S_real_sparse,&
               pexsi_info)
         else
           call f_ppexsi_load_real_symmetric_hs_matrix( pexsi_plan,&
               pexsi_options, n_g_rank, n_g_nonzero, n_l_nonzero, n_l_cols,&
               sparse_pointer, sparse_index, H_real_sparse, 0, S_real_sparse,&
               pexsi_info)
         end if

         if (pexsi_info /= 0) then
            call elsi_stop("PEXSI not able to load H/S matrix.",&
                  "elsi_solve_ev_problem")
         end if

         ! Solve the eigenvalue problem
         call f_ppexsi_dft_driver(pexsi_plan, pexsi_options, n_electrons,&
               mu_Pexsi, n_electrons_pexsi, mu_min_inertia, mu_max_inertia,&
               n_total_inertia_iter, n_total_pexsi_iter, pexsi_info)
         
         if (pexsi_info /= 0) then
            call elsi_stop("PEXSI DFT Driver not able to solve problem.",&
                  "elsi_solve_ev_problem")
         end if

         ! Get the results
         call f_ppexsi_retrieve_real_symmetric_dft_matrix( pexsi_plan,&
            D_real_sparse, ED_real_sparse, FD_real_sparse, e_tot_H, &
            e_tot_S, f_tot, pexsi_info)

         if (pexsi_info /= 0) then
            call elsi_stop("PEXSI not able to retrieve solution.",&
                  "elsi_solve_ev_problem")
         end if

         if( myid == 0 ) then
           write(*,*) "Total energy (H*DM)         = ", e_tot_H
           write(*,*) "Total energy (S*EDM)        = ", e_tot_S
           write(*,*) "Total free energy           = ", f_tot
         endif

      case DEFAULT
         call elsi_stop("No method has been chosen. "//&
            "Please choose method ELPA, OMM, or PEXSI","elsi_solve_ev_problem")
   end select

   call MPI_BARRIER(mpi_comm_global, mpierr)
   
   call elsi_stop_solve_evp_time()

end subroutine

!>
!!  This routine gets the eigenvalues
!!
subroutine elsi_get_eigenvalues(eigenvalues_out,n_eigenvalues)


   implicit none

   integer, intent(in) :: n_eigenvalues !< Number of eigenvalues
   real*8, intent(out) :: eigenvalues_out(n_eigenvalues) !< eigenvalues 
   
   select case (method)
      case (ELPA)
            eigenvalues_out(1:n_eigenvalues) = &
               eigenvalues(1:n_eigenvalues)      
      case (OMM_DENSE)
         call elsi_stop("OMM_DENSE not implemented yet!","elsi_get_eigenvalues")
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!","elsi_get_eigenvalues")
      case DEFAULT
         call elsi_stop("No method has been chosen. "// &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_get_eigenvalues")
   end select

end subroutine

!>
!!  This routine gets the sum of eigenvalues a.k.a total energy
!!
subroutine elsi_get_total_energy(e_tot)


   implicit none

   real*8, intent(out) :: e_tot !< eigenvalues 
   
   select case (method)
      case (ELPA)

         e_tot = 2d0 * SUM(eigenvalues(1:n_eigenvectors))      
      case (OMM_DENSE)
         if (myid == 0) then
         end if
         e_tot = 2d0 * total_energy + 2d0 * n_eigenvectors * eta
      case (PEXSI)
         e_tot = e_tot_h     
      case DEFAULT
         call elsi_stop("No method has been chosen. "//&
               "Please choose method ELPA, OMM_DENSE, or PEXSI",&
               "elsi_get_total_energy")
   end select

end subroutine



!>
!!  This routine prepares the matrices based on the dimension and method
!!  specified
!!
subroutine elsi_deallocate_matrices()

   implicit none

   ! Nullify pointers
   if (associated(H_real))    nullify (H_real)
   if (associated(H_complex)) nullify (H_complex)
   if (associated(S_real))    nullify (S_real)
   if (associated(S_complex)) nullify (S_complex)
 
   ! Free Memory
   if (allocated(H_real_target))    deallocate(H_real_target)
   if (allocated(H_complex_target)) deallocate(H_complex_target)
   if (allocated(S_real_target))    deallocate(S_real_target)
   if (allocated(S_complex_target)) deallocate(S_complex_target)
   if (allocated(vectors_real))     deallocate(vectors_real)
   if (allocated(vectors_complex))  deallocate(vectors_complex)
   if (allocated(eigenvalues))      deallocate(eigenvalues)
   if (allocated(H_real_sparse))    deallocate(H_real_sparse)
   if (allocated(S_real_sparse))    deallocate(S_real_sparse)
   if (allocated(D_real_sparse))    deallocate(D_real_sparse)
   if (allocated(ED_real_sparse))   deallocate(ED_real_sparse)
   if (allocated(FD_real_sparse))   deallocate(FD_real_sparse)
   if (allocated(sparse_index))     deallocate(sparse_index)
   if (allocated(sparse_pointer))   deallocate(sparse_pointer)

   if (method == OMM_DENSE) then
      call m_deallocate (OMM_H_matrix)
      call m_deallocate (OMM_S_matrix)
      call m_deallocate (OMM_C_matrix)
      call m_deallocate (OMM_D_matrix)
      call m_deallocate (OMM_T_matrix)
   end if

   
end subroutine

!>
!!  This routine prepares the matrices based on the dimension and method
!!  specified
subroutine elsi_finalize()

   implicit none
   include "mpif.h"

   call MPI_BARRIER(mpi_comm_global,mpierr)

   call elsi_deallocate_matrices()

   if (mode == PEXSI) call f_ppexsi_plan_finalize(pexsi_plan, pexsi_info)
   
   call hdf5_finalize()

   if (.not.external_blacs) call elsi_finalize_blacs()

   if (.not.external_mpi)   call elsi_finalize_mpi()

   call elsi_stop_total_time()

   call elsi_print_timers()

end subroutine

!> 
!! This routine transforms a general eigenvalue problem to standart form 
!!
!! Starting from Ha = eSa, we first perform a Cholesky decomposition of S
!! S = U^T U, resulting in Ha = eU^TUa
!!
!! Using 1=U^-1U we define a new standard eigenvalue problem by
!! HU^-1(Ua) = eU^T(Ua) => U^-THU^-1 (Ua) = e(Ua)
!!
subroutine elsi_to_standard_eigenvalue_problem()

   logical :: success  !< Success flag of eigensolver
   real*8,     allocatable :: buffer_real (:,:) !< real valued matrix buffer
   complex*16, allocatable :: buffer_complex (:,:) !< complex valued matrix buffer

   select case (method)
      case (ELPA)
         select case (mode)
            case (COMPLEX_VALUES)
                call elsi_stop ("ELPA general to standard complex EVP not yet "&
                     // " implemented.",& 
                  "elsi_to_standard_eigenvalue_problem")
            case (REAL_VALUES)
               allocate (buffer_real(n_l_rows, n_l_cols))
               buffer_real = 0d0
               
               ! Compute U
               ! S contains then U
               call cholesky_real( &
                     n_g_rank, S_real, &
                     n_l_rows, n_b_rows, n_l_cols,&
                     mpi_comm_row, mpi_comm_col, .False., success)

               ! compute U^-1
               ! S contains U^-1
               call invert_trm_real( &
                     n_g_rank, S_real, &
                     n_l_rows, n_b_rows, n_l_cols,&
                     mpi_comm_row, mpi_comm_col, .False., success)

               ! compute H U^-1 -> buffer
               call pdgemm('N','N', n_g_rank, n_g_rank, n_g_rank,&
                     1.0d0, H_real, 1, 1, sc_desc,&
                     S_real, 1, 1, sc_desc,&
                     0.0d0, buffer_real, 1, 1, sc_desc)

               ! compute U^-TH by (HU^-1)^T 
               call pdtran(n_g_rank, n_g_rank, &
                  1.d0, buffer_real, 1, 1, sc_desc, &
                  0.d0, H_real, 1, 1, sc_desc)

               ! compute (U^-TH)U^-1
               buffer_real = H_real
               call pdgemm('N','N', n_g_rank, n_g_rank, n_g_rank,&
                     1.0d0, buffer_real, 1, 1, sc_desc,&
                     S_real, 1, 1, sc_desc,&
                     0.0d0, H_real, 1, 1,sc_desc)

         end select
         if (.not.success) then
            call elsi_stop ("ELPA general to standard EVP failed",& 
                  "elsi_to_standard_eigenvalue_problem")
         end if

      case (OMM_DENSE)
         ! Nothing to be done, supported by OMM
      case (PEXSI)
         call elsi_stop ("PEXSI not yet implemented",& 
                  "elsi_to_standard_eigenvalue_problem")
      case DEFAULT
         call elsi_stop ("No method has been chosen. " // &
            "Please choose method ELPA, OMM_DENSE, or PEXSI",& 
            "elsi_to_standard_eigenvalue_problem")
   end select

   if (allocated(buffer_real))    deallocate (buffer_real)
   if (allocated(buffer_complex)) deallocate (buffer_complex)

end subroutine

!>
!! This routine checks if the obtained eigenvalues and eigenvectors are indeed 
!! solution of the eigenvalue problem HC = eSC. To use this function set the 
!! hamiltonian and the overlap as in the beginning, the eigenvectors are taken
!! from the elsi module
!!
subroutine elsi_check_solution(success)

   implicit none

   logical, intent(out)    :: success  !< Success flag of eigensolver
   real*8,     allocatable :: buffer1_real (:,:) !< real valued matrix buffer
   real*8,     allocatable :: buffer2_real (:,:) !< real valued matrix buffer
   real*8,     allocatable :: buffer3_real (:,:) !< real valued matrix buffer
   complex*16, allocatable :: buffer1_complex (:,:) !< complex valued matrix buffer
   complex*16, allocatable :: buffer2_complex (:,:) !< complex valued matrix buffer
   real*8                  :: norm !< norm of Hc - eSc
   integer                 :: i_val !< eigenvalue index

   success = .True.

   select case (method)
      case (ELPA)
         select case (mode)
            case (COMPLEX_VALUES)
                write(*,'(a)') "COMPLEX check not implemented yet!"
                stop
            case (REAL_VALUES)
               allocate (buffer1_real(n_l_rows, n_l_cols))
               allocate (buffer2_real(n_l_rows, n_l_cols))
               allocate (buffer3_real(n_l_rows, n_l_cols))
               buffer1_real = 0d0
               ! HC
               call pdgemm('N','N', n_g_rank, n_eigenvectors, 1.0d0, H_real, &
                     1, 1, sc_desc, vectors_real, 1, 1, sc_desc, &
                     0.0d0, buffer1_real, 1, 1, sc_desc)

               buffer2_real = vectors_real
               do i_val = 1, n_eigenvectors
                  call pdscal(n_g_rank, eigenvalues(i_val), &
                        buffer2_real, 1, i_val, sc_desc, 1)
               end do

               call pdgemm('N','N', n_g_rank, buffer2_real, 1.0d0, S_real, &
                     1, 1, sc_desc, vectors_real, 1, 1, sc_desc, &
                     0.0d0, buffer3_real, 1, 1, sc_desc)
               
               buffer3_real = buffer3_real - buffer1_real

               do i_val = 1, n_eigenvectors
                  call pdnrm2(n_g_rank, norm, buffer3_real, 1, &
                        i_val, sc_desc, 1)
                  if (norm > 1d-10) then
                    success = .False.
                    if (myid == 0) print *, "Large Residuum: ", norm
                  end if
               end do
               deallocate(buffer1_real)
               deallocate(buffer2_real)
               deallocate(buffer3_real)
         end select
         if (.not.success) then
            call elsi_stop("ELPA failed.","elsi_check_solution")
         end if

      case (OMM_DENSE)
         call elsi_stop("OMM_DENSE not implemented yet!","elsi_check_solution")
      case (PEXSI)
         call elsi_stop("PEXSI not implemented yet!","elsi_check_solution")
      case DEFAULT
         call elsi_stop("No method has been chosen. "//&
            "Please choose method ELPA, OMM_DENSE, or PEXSI",&
            "elsi_check_solution")
   end select

end subroutine


end module ELSI
