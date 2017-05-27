! Copyright (c) 2015-2017, the ELSI team. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!  * Neither the name of the "ELectronic Structure Infrastructure" project nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
! OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!>
!! This module contains utilities for ELSI, including array allocations and
!! printings for debugging.
!!
module ELSI_UTILS

   use iso_c_binding
   use ELSI_CONSTANTS, only: AUTO,ELPA,LIBOMM,PEXSI,CHESS,SIPS,BLACS_DENSE,MULTI_PROC,N_SOLVERS,&
                             N_MATRIX_DATA_TYPES,N_MATRIX_STORAGE_FORMATS,N_PARALLEL_MODES,UNSET
   use ELSI_DIMENSIONS, only: elsi_handle,print_info
   use ELSI_PRECISION, only: r8,i4
   use f_ppexsi_interface
   use MatrixSwitch, only: matrix,m_register_pdbc,m_deallocate
   use m_qetsc

   implicit none
   private

   public :: elsi_set_hamiltonian
   public :: elsi_set_overlap
   public :: elsi_set_eigenvalue
   public :: elsi_set_eigenvector
   public :: elsi_set_density_matrix
   public :: elsi_set_row_ind
   public :: elsi_set_col_ptr
   public :: elsi_set_sparse_hamiltonian
   public :: elsi_set_sparse_overlap
   public :: elsi_set_sparse_density_matrix
   public :: elsi_statement_print
   public :: elsi_allocate
   public :: elsi_cleanup
   public :: elsi_stop
   public :: elsi_check
   public :: elsi_check_handle
   public :: elsi_get_global_row
   public :: elsi_get_global_col
   public :: elsi_get_local_nnz

   interface elsi_set_hamiltonian
      module procedure elsi_set_real_hamiltonian,&
                       elsi_set_complex_hamiltonian
   end interface

   interface elsi_set_sparse_hamiltonian
      module procedure elsi_set_sparse_real_hamiltonian,&
                       elsi_set_sparse_complex_hamiltonian
   end interface

   interface elsi_set_overlap
      module procedure elsi_set_real_overlap,&
                       elsi_set_complex_overlap
   end interface

   interface elsi_set_sparse_overlap
      module procedure elsi_set_sparse_real_overlap,&
                       elsi_set_sparse_complex_overlap
   end interface

   interface elsi_set_eigenvector
      module procedure elsi_set_real_eigenvector,&
                       elsi_set_complex_eigenvector
   end interface

   interface elsi_allocate
      module procedure elsi_allocate_integer_1d,&
                       elsi_allocate_integer_2d,&
                       elsi_allocate_integer_3d,&
                       elsi_allocate_real_1d,&
                       elsi_allocate_real_2d,&
                       elsi_allocate_real_3d,&
                       elsi_allocate_complex_1d,&
                       elsi_allocate_complex_2d,&
                       elsi_allocate_complex_3d
   end interface

contains

!>
!! This routine prints a statement.
!!
subroutine elsi_statement_print(message,elsi_h)

   implicit none

   character(len=*),  intent(in) :: message
   type(elsi_handle), intent(in) :: elsi_h

   if(print_info) then
      if(elsi_h%myid == 0) then
         write(*,'(A)') trim(message)
      endif
   endif

end subroutine

!>
!! This routine allocates a 1D array with real(kind=r8).
!!
subroutine elsi_allocate_real_1d(elsi_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle),          intent(in)    :: elsi_h
   real(kind=r8), allocatable, intent(inout) :: array(:)  !< Data
   integer(kind=i4),           intent(in)    :: dim1      !< Size
   character(len=*),           intent(in)    :: arrayname !< Name
   character(len=*),           intent(in)    :: caller    !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1),stat=error)

   if(error > 0) then 
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,elsi_h,caller)
   endif

   array = 0.0_r8

end subroutine

!>
!! This routine allocates a 1D array with integer.
!!
subroutine elsi_allocate_integer_1d(elsi_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h
   integer(kind=i4), allocatable, intent(inout) :: array(:)  !< Data
   integer(kind=i4),              intent(in)    :: dim1      !< Size
   character(len=*),              intent(in)    :: arrayname !< Name
   character(len=*),              intent(in)    :: caller    !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,elsi_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 1D array with complex(kind=r8).
!!
subroutine elsi_allocate_complex_1d(elsi_h,array,dim1,arrayname,caller)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h
   complex(kind=r8), allocatable, intent(inout) :: array(:)  !< Data
   integer(kind=i4),              intent(in)    :: dim1      !< Size
   character(len=*),              intent(in)    :: arrayname !< Name
   character(len=*),              intent(in)    :: caller    !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,elsi_h,caller)
   endif

   array = cmplx(0.0_r8,0.0_r8)

end subroutine

!>
!! This routine allocates a 2D array of real(kind=r8).
!!
subroutine elsi_allocate_real_2d(elsi_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle),          intent(in)    :: elsi_h
   real(kind=r8), allocatable, intent(inout) :: array(:,:) !< Data
   integer(kind=i4),           intent(in)    :: dim1       !< Size
   integer(kind=i4),           intent(in)    :: dim2       !< Size
   character(len=*),           intent(in)    :: arrayname  !< Name
   character(len=*),           intent(in)    :: caller     !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,elsi_h,caller)
   endif

   array = 0.0_r8

end subroutine

!>
!! This routine allocates a 2D array of integer.
!!
subroutine elsi_allocate_integer_2d(elsi_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h
   integer(kind=i4), allocatable, intent(inout) :: array(:,:) !< Data
   integer(kind=i4),              intent(in)    :: dim1       !< Size
   integer(kind=i4),              intent(in)    :: dim2       !< Size
   character(len=*),              intent(in)    :: arrayname  !< Name
   character(len=*),              intent(in)    :: caller     !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,elsi_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 2D array of complex(kind=r8).
!!
subroutine elsi_allocate_complex_2d(elsi_h,array,dim1,dim2,arrayname,caller)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h
   complex(kind=r8), allocatable, intent(inout) :: array(:,:) !< Data
   integer(kind=i4),              intent(in)    :: dim1       !< Size
   integer(kind=i4),              intent(in)    :: dim2       !< Size
   character(len=*),              intent(in)    :: arrayname  !< Name
   character(len=*),              intent(in)    :: caller     !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,elsi_h,caller)
   endif

   array = cmplx(0.0_r8,0.0_r8)

end subroutine

!>
!! This routine allocates a 3D array of real(kind=r8).
!!
subroutine elsi_allocate_real_3d(elsi_h,array,dim1,dim2,dim3,arrayname,caller)

   implicit none

   type(elsi_handle),          intent(in)    :: elsi_h
   real(kind=r8), allocatable, intent(inout) :: array(:,:,:) !< Data
   integer(kind=i4),           intent(in)    :: dim1         !< Size
   integer(kind=i4),           intent(in)    :: dim2         !< Size
   integer(kind=i4),           intent(in)    :: dim3         !< Size
   character(len=*),           intent(in)    :: arrayname    !< Name
   character(len=*),           intent(in)    :: caller       !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1,dim2,dim3),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,elsi_h,caller)
   endif

   array = 0.0_r8

end subroutine

!>
!! This routine allocates a 3D array of integer.
!!
subroutine elsi_allocate_integer_3d(elsi_h,array,dim1,dim2,dim3,arrayname,caller)

   implicit none

   type(elsi_handle),             intent(in)    :: elsi_h
   integer(kind=i4), allocatable, intent(inout) :: array(:,:,:) !< Data
   integer(kind=i4),              intent(in)    :: dim1         !< Size
   integer(kind=i4),              intent(in)    :: dim2         !< Size
   integer(kind=i4),              intent(in)    :: dim3         !< Size
   character(len=*),              intent(in)    :: arrayname    !< Name
   character(len=*),              intent(in)    :: caller       !< Caller
   
   integer :: error
   character*200 :: message
   
   allocate(array(dim1,dim2,dim3),stat=error)
         
   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,elsi_h,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 3D array of complex(kind=r8).
!!
subroutine elsi_allocate_complex_3d(elsi_h,array,dim1,dim2,dim3,arrayname,caller)

   implicit none
   

   type(elsi_handle),             intent(in)    :: elsi_h
   complex(kind=r8), allocatable, intent(inout) :: array(:,:,:) !< Data
   integer(kind=i4),              intent(in)    :: dim1         !< Size
   integer(kind=i4),              intent(in)    :: dim2         !< Size
   integer(kind=i4),              intent(in)    :: dim3         !< Size
   character(len=*),              intent(in)    :: arrayname    !< Name
   character(len=*),              intent(in)    :: caller       !< Caller
   
   integer :: error
   character*200 :: message
   
   allocate(array(dim1,dim2,dim3),stat=error)
         
   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,elsi_h,caller)
   endif

   array = cmplx(0.0_r8,0.0_r8)

end subroutine

!>
!! Clean shutdown in case of error.
!!
subroutine elsi_stop(message,elsi_h,caller)

   implicit none
   include "mpif.h"

   character(len=*),  intent(in) :: message
   type(elsi_handle), intent(in) :: elsi_h
   character(len=*),  intent(in) :: caller

   character(len=4096) :: string_message
   integer :: i_task
   integer :: mpierr

   if(elsi_h%mpi_is_setup) then
      do i_task = 0,elsi_h%n_procs-1
         if(elsi_h%myid == i_task) then
            write(string_message,"(1X,'*** Proc',I5,' in ',A,': ',A)")&
                  elsi_h%myid,trim(caller),trim(message)
            write(*,'(A)') trim(string_message)
         endif
         call MPI_Barrier(elsi_h%mpi_comm,mpierr)
      enddo

      if(elsi_h%n_procs > 1) then
         call MPI_Abort(elsi_h%mpi_comm,0,mpierr)
      endif
   else
      write(string_message,"(1X,'*** Proc  N/A',' in ',A,': ',A)")&
            trim(caller),trim(message)
      write(*,'(A)') trim(string_message)
   endif

   stop

end subroutine

!>
!! This routine sets the real hamiltonian matrix.
!!
subroutine elsi_set_real_hamiltonian(elsi_h,H_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   real(kind=r8),     target        :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols)

   character*40, parameter :: caller = "elsi_set_real_hamiltonian"

   select case (elsi_h%solver)
      case (ELPA)
         elsi_h%ham_real => H_in
      case (LIBOMM)
         call m_register_pdbc(elsi_h%ham_omm,H_in,elsi_h%sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine sets the complex hamiltonian matrix.
!!
subroutine elsi_set_complex_hamiltonian(elsi_h,H_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   complex(kind=r8),  target        :: H_in(elsi_h%n_l_rows,elsi_h%n_l_cols)

   character*40, parameter :: caller = "elsi_set_complex_hamiltonian"

   select case (elsi_h%solver)
      case (ELPA)
         elsi_h%ham_complex => H_in
      case (LIBOMM)
         call m_register_pdbc(elsi_h%ham_omm,H_in,elsi_h%sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine sets the sparse real Hamiltonian matrix.
!!
subroutine elsi_set_sparse_real_hamiltonian(elsi_h,H_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   real(kind=r8),     target        :: H_in(elsi_h%nnz_l_pexsi)

   character*40, parameter :: caller = "elsi_set_sparse_real_hamiltonian"

   select case (elsi_h%solver)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         elsi_h%ham_real_ccs => H_in
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         elsi_h%ham_real_ccs => H_in
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine sets the sparse complex Hamiltonian matrix.
!!
subroutine elsi_set_sparse_complex_hamiltonian(elsi_h,H_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   complex(kind=r8),  target        :: H_in(elsi_h%nnz_l_pexsi)

   character*40, parameter :: caller = "elsi_set_sparse_complex_hamiltonian"

   select case (elsi_h%solver)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         elsi_h%ham_complex_ccs => H_in         
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine sets the real overlap matrix.
!!
subroutine elsi_set_real_overlap(elsi_h,S_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   real(kind=r8),     target        :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols)

   character*40, parameter :: caller = "elsi_set_real_overlap"

   if(.not. elsi_h%overlap_is_unit) then
      select case (elsi_h%solver)
         case (ELPA)
            elsi_h%ovlp_real => S_in
         case (LIBOMM)
            call m_register_pdbc(elsi_h%ovlp_omm,S_in,elsi_h%sc_desc)
         case (PEXSI)
            ! Nothing to be done here
         case (CHESS)
            call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
         case (SIPS)
            ! Nothing to be done here
         case DEFAULT
            call elsi_stop(" No supported solver has been chosen."//&
                           " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                           " Exiting...",elsi_h,caller)
      end select
   endif

end subroutine

!>
!! This routine sets the complex overlap matrix.
!!
subroutine elsi_set_complex_overlap(elsi_h,S_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   complex(kind=r8),  target        :: S_in(elsi_h%n_l_rows,elsi_h%n_l_cols)

   character*40, parameter :: caller = "elsi_set_complex_overlap"

   if(.not. elsi_h%overlap_is_unit) then
      select case (elsi_h%solver)
         case (ELPA)
            elsi_h%ovlp_complex => S_in
         case (LIBOMM)
            call m_register_pdbc(elsi_h%ovlp_omm,S_in,elsi_h%sc_desc)
         case (PEXSI)
            ! Nothing to be done here
         case (CHESS)
            call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
         case (SIPS)
            call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
         case DEFAULT
            call elsi_stop(" No supported solver has been chosen."//&
                           " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                           " Exiting...",elsi_h,caller)
      end select
   endif
end subroutine

!>
!! This routine sets the sparse real overlap matrix.
!!
subroutine elsi_set_sparse_real_overlap(elsi_h,S_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   real(kind=r8),     target        :: S_in(elsi_h%nnz_l_pexsi)

   character*40, parameter :: caller = "elsi_set_sparse_real_overlap"

   if(.not. elsi_h%overlap_is_unit) then
      select case (elsi_h%solver)
         case (ELPA)
            ! Nothing to be done here
         case (LIBOMM)
            ! Nothing to be done here
         case (PEXSI)
            elsi_h%ovlp_real_ccs => S_in
         case (CHESS)
            call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
         case (SIPS)
            elsi_h%ovlp_real_ccs => S_in
         case DEFAULT
            call elsi_stop(" No supported solver has been chosen."//&
                           " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                           " Exiting...",elsi_h,caller)
      end select
   endif

end subroutine

!>
!! This routine sets the sparse complex overlap matrix.
!!
subroutine elsi_set_sparse_complex_overlap(elsi_h,S_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   complex(kind=r8),  target        :: S_in(elsi_h%nnz_l_pexsi)

   character*40, parameter :: caller = "elsi_set_sparse_complex_overlap"

   if(.not. elsi_h%overlap_is_unit) then
      select case (elsi_h%solver)
         case (ELPA)
            ! Nothing to be done here
         case (LIBOMM)
            ! Nothing to be done here
         case (PEXSI)
            elsi_h%ovlp_complex_ccs => S_in
         case (CHESS)
            call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
         case (SIPS)
            call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
         case DEFAULT
            call elsi_stop(" No supported solver has been chosen."//&
                           " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                           " Exiting...",elsi_h,caller)
      end select
   endif

end subroutine

!>
!! This routine sets the eigenvalues.
!!
subroutine elsi_set_eigenvalue(elsi_h,e_val_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   real(kind=r8),     target        :: e_val_in(elsi_h%n_g_size)

   character*40, parameter :: caller = "elsi_set_eigenvalue"

   select case (elsi_h%solver)
      case (ELPA)
         elsi_h%eval => e_val_in
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         elsi_h%eval => e_val_in
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine sets the real eigenvectors.
!!
subroutine elsi_set_real_eigenvector(elsi_h,e_vec_in)

   implicit none
   
   type(elsi_handle), intent(inout) :: elsi_h
   real(kind=r8),     target        :: e_vec_in(elsi_h%n_l_rows,elsi_h%n_l_cols)

   character*40, parameter :: caller = "elsi_set_real_eigenvector"

   select case (elsi_h%solver)
      case (ELPA)
         elsi_h%evec_real => e_vec_in
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         elsi_h%evec_real => e_vec_in
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine sets the complex eigenvectors.
!!
subroutine elsi_set_complex_eigenvector(elsi_h,e_vec_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   complex(kind=r8),  target        :: e_vec_in(elsi_h%n_l_rows,elsi_h%n_l_cols)

   character*40, parameter :: caller = "elsi_set_complex_eigenvector"

   select case (elsi_h%solver)
      case (ELPA)
         elsi_h%evec_complex => e_vec_in
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine sets the density matrix.
!!
subroutine elsi_set_density_matrix(elsi_h,D_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   real(kind=r8),     target        :: D_in(elsi_h%n_l_rows,elsi_h%n_l_cols)

   character*40, parameter :: caller = "elsi_set_density_matrix"

   select case (elsi_h%solver)
      case (ELPA)
         elsi_h%den_mat => D_in
      case (LIBOMM)
         call m_register_pdbc(elsi_h%den_mat_omm,D_in,elsi_h%sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine sets the sparse density matrix.
!!
subroutine elsi_set_sparse_density_matrix(elsi_h,D_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   real(kind=r8),     target        :: D_in(elsi_h%nnz_l_pexsi)

   character*40, parameter :: caller = "elsi_set_sparse_density_matrix"

   select case (elsi_h%solver)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         elsi_h%den_mat_ccs => D_in
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",elsi_h,caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine sets the row index.
!!
subroutine elsi_set_row_ind(elsi_h,row_ind_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   integer(kind=i4),  target        :: row_ind_in(elsi_h%nnz_l_pexsi)

   character*40, parameter :: caller = "elsi_set_row_ind"

   select case (elsi_h%solver)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         elsi_h%row_ind_ccs => row_ind_in
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         elsi_h%row_ind_ccs => row_ind_in
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine sets the column pointer.
!!
subroutine elsi_set_col_ptr(elsi_h,col_ptr_in)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h
   integer(kind=i4),  target        :: col_ptr_in(elsi_h%n_l_cols_pexsi+1)

   character*40, parameter :: caller = "elsi_set_col_ptr"

   select case (elsi_h%solver)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         elsi_h%col_ptr_ccs => col_ptr_in
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",elsi_h,caller)
      case (SIPS)
         elsi_h%col_ptr_ccs => col_ptr_in
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",elsi_h,caller)
   end select

end subroutine

!>
!! This routine frees memory.
!!
subroutine elsi_cleanup(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h

   character*40, parameter :: caller = "elsi_cleanup"

   ! Nullify pointers
   if(associated(elsi_h%ham_real))         nullify(elsi_h%ham_real)
   if(associated(elsi_h%ham_complex))      nullify(elsi_h%ham_complex)
   if(associated(elsi_h%ovlp_real))        nullify(elsi_h%ovlp_real)
   if(associated(elsi_h%ovlp_complex))     nullify(elsi_h%ovlp_complex)
   if(associated(elsi_h%evec_real))        nullify(elsi_h%evec_real)
   if(associated(elsi_h%evec_complex))     nullify(elsi_h%evec_complex)
   if(associated(elsi_h%eval))             nullify(elsi_h%eval)
   if(associated(elsi_h%den_mat))          nullify(elsi_h%den_mat)
   if(associated(elsi_h%ham_real_ccs))     nullify(elsi_h%ham_real_ccs)
   if(associated(elsi_h%ham_complex_ccs))  nullify(elsi_h%ham_complex_ccs)
   if(associated(elsi_h%ovlp_real_ccs))    nullify(elsi_h%ovlp_real_ccs)
   if(associated(elsi_h%ovlp_complex_ccs)) nullify(elsi_h%ovlp_complex_ccs)
   if(associated(elsi_h%den_mat_ccs))      nullify(elsi_h%den_mat_ccs)
   if(associated(elsi_h%row_ind_ccs))      nullify(elsi_h%row_ind_ccs)
   if(associated(elsi_h%col_ptr_ccs))      nullify(elsi_h%col_ptr_ccs)

   ! ELPA
   if(allocated(elsi_h%ham_real_elpa))     deallocate(elsi_h%ham_real_elpa)
   if(allocated(elsi_h%ham_complex_elpa))  deallocate(elsi_h%ham_complex_elpa)
   if(allocated(elsi_h%ovlp_real_elpa))    deallocate(elsi_h%ovlp_real_elpa)
   if(allocated(elsi_h%ovlp_complex_elpa)) deallocate(elsi_h%ovlp_complex_elpa)
   if(allocated(elsi_h%evec_real_elpa))    deallocate(elsi_h%evec_real_elpa)
   if(allocated(elsi_h%evec_complex_elpa)) deallocate(elsi_h%evec_complex_elpa)
   if(allocated(elsi_h%eval_elpa))         deallocate(elsi_h%eval_elpa)
   if(allocated(elsi_h%den_mat_elpa))      deallocate(elsi_h%den_mat_elpa)
   if(allocated(elsi_h%occ_elpa))          deallocate(elsi_h%occ_elpa)

   ! libOMM
   if(elsi_h%ham_omm%is_initialized)       call m_deallocate(elsi_h%ham_omm)
   if(elsi_h%ovlp_omm%is_initialized)      call m_deallocate(elsi_h%ovlp_omm)
   if(elsi_h%den_mat_omm%is_initialized)   call m_deallocate(elsi_h%den_mat_omm)
   if(elsi_h%coeff_omm%is_initialized)     call m_deallocate(elsi_h%coeff_omm)
   if(elsi_h%t_den_mat_omm%is_initialized) call m_deallocate(elsi_h%t_den_mat_omm)

   ! PEXSI
   if(allocated(elsi_h%ham_real_pexsi))    deallocate(elsi_h%ham_real_pexsi)
   if(allocated(elsi_h%ovlp_real_pexsi))   deallocate(elsi_h%ovlp_real_pexsi)
   if(allocated(elsi_h%den_mat_pexsi))     deallocate(elsi_h%den_mat_pexsi)
   if(allocated(elsi_h%e_den_mat_pexsi))   deallocate(elsi_h%e_den_mat_pexsi)
   if(allocated(elsi_h%f_den_mat_pexsi))   deallocate(elsi_h%f_den_mat_pexsi)
   if(allocated(elsi_h%row_ind_pexsi))     deallocate(elsi_h%row_ind_pexsi)
   if(allocated(elsi_h%col_ptr_pexsi))     deallocate(elsi_h%col_ptr_pexsi)

   ! SIPs
   if(allocated(elsi_h%inertias))          deallocate(elsi_h%inertias)
   if(allocated(elsi_h%shifts))            deallocate(elsi_h%shifts)
   if(allocated(elsi_h%slices))            deallocate(elsi_h%slices)

   if(allocated(elsi_h%local_row))         deallocate(elsi_h%local_row)
   if(allocated(elsi_h%local_col))         deallocate(elsi_h%local_col)

   ! Finalize PEXSI plan
   if(elsi_h%solver == PEXSI) then
      call f_ppexsi_plan_finalize(elsi_h%pexsi_plan,elsi_h%pexsi_info)
   endif

   ! Finalize QETSC
   if(elsi_h%solver == SIPS) then
      call finalize_sips()
   endif

   ! Reset elsi_h
   elsi_h%solver                  = UNSET
   elsi_h%matrix_data_type        = UNSET
   elsi_h%matrix_storage_format   = UNSET
   elsi_h%parallel_mode           = UNSET
   elsi_h%n_elsi_calls            = 0
   elsi_h%n_g_size                = UNSET
   elsi_h%n_b_rows                = UNSET
   elsi_h%n_b_cols                = UNSET
   elsi_h%n_p_rows                = UNSET
   elsi_h%n_p_cols                = UNSET
   elsi_h%n_l_rows                = UNSET
   elsi_h%n_l_cols                = UNSET
   elsi_h%myid                    = UNSET
   elsi_h%n_procs                 = UNSET
   elsi_h%mpi_comm                = UNSET
   elsi_h%mpi_is_setup            = .false.
   elsi_h%blacs_ctxt              = UNSET
   elsi_h%sc_desc                 = UNSET
   elsi_h%mpi_comm_row            = UNSET
   elsi_h%mpi_comm_col            = UNSET
   elsi_h%my_p_row                = UNSET
   elsi_h%my_p_col                = UNSET
   elsi_h%blacs_is_setup          = .false.
   elsi_h%nnz_g                   = UNSET
   elsi_h%nnz_l                   = UNSET
   elsi_h%zero_threshold          = 1.0e-15_r8
   elsi_h%overlap_is_unit         = .false.
   elsi_h%overlap_is_singular     = .false.
   elsi_h%no_singularity_check    = .false.
   elsi_h%singularity_tolerance   = 1.0e-5_r8
   elsi_h%stop_singularity        = .false.
   elsi_h%n_nonsingular           = UNSET
   elsi_h%n_electrons             = 0.0_r8
   elsi_h%n_states                = UNSET
   elsi_h%n_occupied_states       = UNSET
   elsi_h%broadening_scheme       = 0
   elsi_h%broadening_width        = 1.0e-2_r8
   elsi_h%occ_tolerance           = 1.0e-13_r8
   elsi_h%max_mu_steps            = 100
   elsi_h%spin_degen              = 0.0_r8
   elsi_h%elpa_one_always         = .false.
   elsi_h%elpa_two_always         = .false.
   elsi_h%n_elpa_steps            = UNSET
   elsi_h%new_overlap             = .true.
   elsi_h%coeff_initialized       = .false.
   elsi_h%total_energy            = 0.0_r8
   elsi_h%omm_flavor              = UNSET
   elsi_h%scale_kinetic           = 0.0_r8
   elsi_h%calc_ed                 = .false.
   elsi_h%eta                     = 0.0_r8
   elsi_h%min_tol                 = 1.0e-12_r8
   elsi_h%nk_times_nspin          = UNSET
   elsi_h%i_k_spin                = UNSET
   elsi_h%omm_verbose             = .false.
   elsi_h%do_dealloc              = .false.
   elsi_h%use_psp                 = .false.
   elsi_h%my_p_row_pexsi          = UNSET
   elsi_h%my_p_col_pexsi          = UNSET
   elsi_h%n_b_rows_pexsi          = UNSET
   elsi_h%n_b_cols_pexsi          = UNSET
   elsi_h%n_p_rows_pexsi          = UNSET
   elsi_h%n_p_cols_pexsi          = UNSET
   elsi_h%n_l_rows_pexsi          = UNSET
   elsi_h%n_l_cols_pexsi          = UNSET
   elsi_h%n_p_per_pole_pexsi      = UNSET
   elsi_h%nnz_l_pexsi             = UNSET
   elsi_h%sparsity_pattern_ready  = .false.
   elsi_h%n_p_per_pole_ready      = .false.
   elsi_h%small_pexsi_tol         = .false.
   elsi_h%final_pexsi_tol         = 1.0e-2_r8
   elsi_h%pexsi_info              = UNSET
   elsi_h%pexsi_output_file_index = UNSET
   elsi_h%mu_pexsi                = 0.0_r8
   elsi_h%n_electrons_pexsi       = 0.0_r8
   elsi_h%mu_min_inertia          = 0.0_r8
   elsi_h%mu_max_inertia          = 0.0_r8
   elsi_h%n_total_inertia_iter    = UNSET
   elsi_h%n_total_pexsi_iter      = UNSET
   elsi_h%e_tot_h                 = 0.0_r8
   elsi_h%e_tot_s                 = 0.0_r8
   elsi_h%f_tot                   = 0.0_r8
   elsi_h%n_b_rows_sips           = UNSET
   elsi_h%n_b_cols_sips           = UNSET
   elsi_h%n_l_rows_sips           = UNSET
   elsi_h%n_l_cols_sips           = UNSET
   elsi_h%nnz_l_sips              = UNSET
   elsi_h%n_p_per_slice_sips      = UNSET
   elsi_h%n_inertia_steps         = UNSET
   elsi_h%n_solve_steps           = UNSET
   elsi_h%slicing_method          = UNSET
   elsi_h%inertia_option          = UNSET
   elsi_h%unbound                 = UNSET
   elsi_h%n_slices                = UNSET
   elsi_h%interval                = 0.0_r8
   elsi_h%slice_buffer            = 0.0_r8

end subroutine

!>
!! This routine guarantees that there are no mutually conflicting parameters.
!!
subroutine elsi_check(elsi_h,caller)

   implicit none

   type(elsi_handle), intent(in) :: elsi_h
   character(len=*),  intent(in) :: caller

   ! General check of solver, parallel mode, matrix data type, matrix storage format
   if(elsi_h%solver == UNSET) then
      call elsi_stop(" The solver has not been set."//&
                     " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                     " Exiting...",elsi_h,caller)
   else if((elsi_h%solver < 0) .or. (elsi_h%solver .ge. N_SOLVERS)) then
      call elsi_stop(" An unsupported solver has been chosen."//&
                     " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                     " Exiting...",elsi_h,caller)
   endif

   if(elsi_h%parallel_mode == UNSET) then
      call elsi_stop(" The parallel mode has not been set."//&
                     " Please choose either SINGLE_PROC or MULTI_PROC."//&
                     " Exiting...",elsi_h,caller)
   else if((elsi_h%parallel_mode < 0) .or. &
      (elsi_h%parallel_mode .ge. N_PARALLEL_MODES)) then
      call elsi_stop(" An unsupported parallel mode has been chosen."//&
                     " Please choose either SINGLE_PROC or MULTI_PROC."//&
                     " Exiting...",elsi_h,caller)
   endif

   if(elsi_h%matrix_data_type == UNSET) then
      call elsi_stop(" The matrix data type has not been set."//& 
                     " Please choose either REAL_VALUES or COMPLEX_VALUES."//&
                     " Exiting...",elsi_h,caller)
   else if((elsi_h%matrix_data_type < 0) .or. &
      (elsi_h%matrix_data_type .ge. N_MATRIX_DATA_TYPES)) then
      call elsi_stop(" An unsupported matirx data type has been chosen."//&
                     " Please choose either REAL_VALUES or COMPLEX_VALUES."//&
                     " Exiting...",elsi_h,caller)
   endif

   if(elsi_h%matrix_storage_format == UNSET) then
      call elsi_stop(" The matrix storage format has not been set."//&             
                     " Please choose either BLACS_DENSE or PEXSI_CSC."//&
                     " Exiting...",elsi_h,caller)
   else if((elsi_h%matrix_storage_format < 0) .or. &
      (elsi_h%matrix_storage_format .ge. N_MATRIX_STORAGE_FORMATS)) then
      call elsi_stop(" An unsupported matirx storage format has been chosen."//&
                     " Please choose either BLACS_DENSE or PEXSI_CSC."//&
                     " Exiting...",elsi_h,caller)
   endif

   ! Specific check for each solver
   if(elsi_h%solver == AUTO) then
      call elsi_stop(" AUTO not yet available."//&
                     " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                     " Exiting...",elsi_h,caller)

   else if(elsi_h%solver == ELPA) then
      if(elsi_h%parallel_mode == MULTI_PROC) then
         if(.not. elsi_h%mpi_is_setup) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI being"//&
                           " set up before calling the solver."//&
                           " Exiting...",elsi_h,caller)
         endif

         if(.not. elsi_h%blacs_is_setup) then
            call elsi_stop(" MULTI_PROC parallel mode requires BLACS being"//&
                           " set up before calling the solver."//&
                           " Exiting...",elsi_h,caller)
         endif
      endif

      if(elsi_h%matrix_storage_format == PEXSI_CSC) then
         if(.not. elsi_h%sparsity_pattern_ready) then
            call elsi_stop(" The PEXSI_CSC format has been chosen."//&
                           " Please set the sparsity pattern before"//&
                           " calling the solver. Exiting...",elsi_h,caller)
         endif
      endif

   else if(elsi_h%solver == LIBOMM) then
      if(elsi_h%parallel_mode == MULTI_PROC) then
         if(.not. elsi_h%mpi_is_setup) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI being"//&
                           " set up before calling the solver."//&
                           " Exiting...",elsi_h,caller)
         endif

         if(.not. elsi_h%blacs_is_setup) then
            call elsi_stop(" MULTI_PROC parallel mode requires BLACS being"//&
                           " set up before calling the solver."//&
                           " Exiting...",elsi_h,caller)
         endif
      else
         call elsi_stop(" libOMM has been chosen as the solver."//&
                        " Please choose MULTI_PROC parallel mode."//&
                        " Exiting...",elsi_h,caller)
      endif

      if(elsi_h%matrix_storage_format /= BLACS_DENSE) then
         call elsi_stop(" libOMM has been chosen as the solver."//&
                        " Please choose BLACS_DENSE matrix storage format."//&
                        " Exiting...",elsi_h,caller)
      endif

   else if(elsi_h%solver == PEXSI) then
      if(elsi_h%parallel_mode == MULTI_PROC) then
         if(.not. elsi_h%mpi_is_setup) then
            call elsi_stop(" MULTI_PROC parallel mode requires MPI being"//&
                           " set up before calling the solver."//&
                           " Exiting...",elsi_h,caller)
         endif
      else
         call elsi_stop(" PEXSI has been chosen as the solver."//&
                        " Please choose MULTI_PROC parallel mode."//&
                        " Exiting...",elsi_h,caller)
      endif

      if(elsi_h%matrix_storage_format == BLACS_DENSE) then
         if(.not. elsi_h%blacs_is_setup) then
            call elsi_stop(" The BLACS_DENSE format has been chosen."//&
                           " Please set up BLACS before calling the"//&
                           " solver. Exiting...",elsi_h,caller)
         endif
      else
         if(.not. elsi_h%sparsity_pattern_ready) then
            call elsi_stop(" The PEXSI_CSC format has been chosen."//&
                           " Please set the sparsity pattern before"//&
                           " calling the solver. Exiting...",elsi_h,caller)
         endif
      endif

      if(elsi_h%n_g_size < elsi_h%n_p_per_pole_pexsi) then
         call elsi_stop(" PEXSI has been chosen as the solver."//&
                        " The size of matrix is too small for this"//&
                        " number of processes. Exiting...",elsi_h,caller)
      endif

   else if(elsi_h%solver == CHESS) then
      call elsi_stop(" CHESS not yet available."//&
                     " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                     " Exiting...",elsi_h,caller)

   else if(elsi_h%solver == SIPS) then
      call elsi_statement_print("  ATTENTION! SIPS is EXPERIMENTAL.",elsi_h)

      if(elsi_h%n_g_size < elsi_h%n_procs) then
         call elsi_stop(" SIPs has been chosen as the solver."//&
                        " The size of matrix is too small for this"//&
                        " number of processes. Exiting...",elsi_h,caller)
      endif

   else
      call elsi_stop(" No supported solver has been chosen."//&
                     " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                     " Exiting...",elsi_h,caller)
   endif

end subroutine

!>
!! This routine checks whether the input elsi_h has been properly
!! initialized by the user or not.
!!
!! Do NOT call this subroutine outside the (main) ELSI module.
!!
subroutine elsi_check_handle(elsi_h,caller)

   implicit none

   type(elsi_handle), intent(in) :: elsi_h
   character(len=*),  intent(in) :: caller

   if(.not. elsi_h%handle_initialized) then
      call elsi_stop(" Invalid handle! An ELSI handle must be properly"//&
                     " initialized with 'elsi_init' before being used."//&
                     " Exiting...",elsi_h,caller)
   endif

end subroutine

!> 
!! This routine computes global row index based on local row index.
!!
subroutine elsi_get_global_row(elsi_h,global_idx,local_idx)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h
   integer(kind=i4),  intent(in)  :: local_idx  !< Local index
   integer(kind=i4),  intent(out) :: global_idx !< Global index

   integer(kind=i4) :: block !< Local block
   integer(kind=i4) :: idx   !< Local index in block

   character*40, parameter :: caller = "elsi_get_global_row"

   block = (local_idx-1)/elsi_h%n_b_rows
   idx = local_idx-block*elsi_h%n_b_rows

   global_idx = elsi_h%my_p_row*elsi_h%n_b_rows+&
                block*elsi_h%n_b_rows*elsi_h%n_p_rows+idx

end subroutine

!> 
!! This routine computes global column index based on local column index.
!!
subroutine elsi_get_global_col(elsi_h,global_idx,local_idx)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h
   integer(kind=i4),  intent(in)  :: local_idx  !< Local index
   integer(kind=i4),  intent(out) :: global_idx !< Global index

   integer(kind=i4) :: block !< Local block
   integer(kind=i4) :: idx   !< Local index in block

   character*40, parameter :: caller = "elsi_get_global_col"

   block = (local_idx-1)/elsi_h%n_b_cols
   idx = local_idx-block*elsi_h%n_b_cols

   global_idx = elsi_h%my_p_col*elsi_h%n_b_cols+&
                block*elsi_h%n_b_cols*elsi_h%n_p_cols+idx

end subroutine

!>
!! This routine counts the local number of non_zero elements.
!!
subroutine elsi_get_local_nnz(elsi_h,matrix,n_rows,n_cols,nnz)

   implicit none

   type(elsi_handle), intent(in)  :: elsi_h
   real(kind=r8),     intent(in)  :: matrix(n_rows,n_cols) !< Local matrix
   integer(kind=i4),  intent(in)  :: n_rows                !< Local rows
   integer(kind=i4),  intent(in)  :: n_cols                !< Local cols
   integer(kind=i4),  intent(out) :: nnz                   !< Number of non-zero

   integer(kind=i4) :: i_row !< Row counter
   integer(kind=i4) :: i_col !< Column counter
   integer(kind=i4) :: i_nz  !< Non-zero element counter

   character*40, parameter :: caller = "elsi_get_local_nnz"

   nnz = 0

   do i_col = 1,n_cols
      do i_row = 1,n_rows
         if(abs(matrix(i_row,i_col)) > elsi_h%zero_threshold) then
            nnz = nnz+1
         endif
      enddo
   enddo

end subroutine

end module ELSI_UTILS
