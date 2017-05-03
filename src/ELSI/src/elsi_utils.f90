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
   use ELSI_PRECISION, only : dp
   use ELSI_DIMENSIONS
   use MatrixSwitch
   use f_ppexsi_interface
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
subroutine elsi_statement_print(message)

   implicit none

   character(len=*), intent(in) :: message !< Message to print

   if(print_info) then
      if(myid == 0) then
         write(*,'(A)') trim(message)
      endif
   endif

end subroutine

!>
!! This routine allocates a 1D array with real(kind=dp).
!!
subroutine elsi_allocate_real_1d(array,dim1,arrayname,caller)

   implicit none

   real(kind=dp), allocatable, intent(inout) :: array(:)  !< Data
   integer,                    intent(in)    :: dim1      !< Size
   character(len=*),           intent(in)    :: arrayname !< Name
   character(len=*),           intent(in)    :: caller    !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1),stat=error)

   if(error > 0) then 
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,caller)
   endif

   array = 0.0_dp

end subroutine

!>
!! This routine allocates a 1D array with integer.
!!
subroutine elsi_allocate_integer_1d(array,dim1,arrayname,caller)

   implicit none

   integer, allocatable, intent(inout) :: array(:)  !< Data
   integer,              intent(in)    :: dim1      !< Size
   character(len=*),     intent(in)    :: arrayname !< Name
   character(len=*),     intent(in)    :: caller    !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 1D array with complex(kind=dp).
!!
subroutine elsi_allocate_complex_1d(array,dim1,arrayname,caller)

   implicit none

   complex(kind=dp), allocatable, intent(inout) :: array(:)  !< Data
   integer,                       intent(in)    :: dim1      !< Size
   character(len=*),              intent(in)    :: arrayname !< Name
   character(len=*),              intent(in)    :: caller    !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,caller)
   endif

   array = cmplx(0.0_dp,0.0_dp)

end subroutine

!>
!! This routine allocates a 2D array of real(kind=dp).
!!
subroutine elsi_allocate_real_2d(array,dim1,dim2,arrayname,caller)

   implicit none

   real(kind=dp), allocatable, intent(inout) :: array(:,:) !< Data
   integer,                    intent(in)    :: dim1       !< Size
   integer,                    intent(in)    :: dim2       !< Size
   character(len=*),           intent(in)    :: arrayname  !< Name
   character(len=*),           intent(in)    :: caller     !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,caller)
   endif

   array = 0.0_dp

end subroutine

!>
!! This routine allocates a 2D array of integer.
!!
subroutine elsi_allocate_integer_2d(array,dim1,dim2,arrayname,caller)

   implicit none

   integer, allocatable, intent(inout) :: array(:,:) !< Data
   integer,              intent(in)    :: dim1       !< Size
   integer,              intent(in)    :: dim2       !< Size
   character(len=*),     intent(in)    :: arrayname  !< Name
   character(len=*),     intent(in)    :: caller     !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 2D array of complex(kind=dp).
!!
subroutine elsi_allocate_complex_2d(array,dim1,dim2,arrayname,caller)

   implicit none

   complex(kind=dp), allocatable, intent(inout) :: array(:,:) !< Data
   integer,                       intent(in)    :: dim1       !< Size
   integer,                       intent(in)    :: dim2       !< Size
   character(len=*),              intent(in)    :: arrayname  !< Name
   character(len=*),              intent(in)    :: caller     !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1,dim2),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,caller)
   endif

   array = cmplx(0.0_dp,0.0_dp)

end subroutine

!>
!! This routine allocates a 3D array of real(kind=dp).
!!
subroutine elsi_allocate_real_3d(array,dim1,dim2,dim3,arrayname,caller)

   implicit none

   real(kind=dp), allocatable, intent(inout) :: array(:,:,:) !< Data
   integer,                    intent(in)    :: dim1         !< Size
   integer,                    intent(in)    :: dim2         !< Size
   integer,                    intent(in)    :: dim3         !< Size
   character(len=*),           intent(in)    :: arrayname    !< Name
   character(len=*),           intent(in)    :: caller       !< Caller

   integer :: error
   character*200 :: message

   allocate(array(dim1,dim2,dim3),stat=error)

   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,caller)
   endif

   array = 0.0_dp

end subroutine

!>
!! This routine allocates a 3D array of integer.
!!
subroutine elsi_allocate_integer_3d(array,dim1,dim2,dim3,arrayname,caller)

   implicit none
   
   integer, allocatable, intent(inout) :: array(:,:,:) !< Data
   integer,              intent(in)    :: dim1         !< Size
   integer,              intent(in)    :: dim2         !< Size
   integer,              intent(in)    :: dim3         !< Size
   character(len=*),     intent(in)    :: arrayname    !< Name
   character(len=*),     intent(in)    :: caller       !< Caller
   
   integer :: error
   character*200 :: message
   
   allocate(array(dim1,dim2,dim3),stat=error)
         
   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,caller)
   endif

   array = 0

end subroutine

!>
!! This routine allocates a 3D array of complex(kind=dp).
!!
subroutine elsi_allocate_complex_3d(array,dim1,dim2,dim3,arrayname,caller)

   implicit none
   
   complex(kind=dp), allocatable, intent(inout) :: array(:,:,:) !< Data
   integer,                       intent(in)    :: dim1         !< Size
   integer,                       intent(in)    :: dim2         !< Size
   integer,                       intent(in)    :: dim3         !< Size
   character(len=*),              intent(in)    :: arrayname    !< Name
   character(len=*),              intent(in)    :: caller       !< Caller
   
   integer :: error
   character*200 :: message
   
   allocate(array(dim1,dim2,dim3),stat=error)
         
   if(error > 0) then
      write(message,"(A,A)") "Error in allocating ",trim(arrayname)
      call elsi_stop(message,caller)
   endif

   array = cmplx(0.0_dp,0.0_dp)

end subroutine

!>
!! Clean shutdown in case of error.
!!
subroutine elsi_stop(message,caller)

   implicit none
   include "mpif.h"

   character(len=*), intent(in) :: message !< Error message
   character(len=*), intent(in) :: caller  !< Caller

   character(len=4096) :: string_message
   integer :: i_task

   if (mpi_is_setup) then
      do i_task = 0,n_procs-1
         if(myid == i_task) then
            write(string_message,"(1X,'*** Proc',I5,' in ',A,': ',A)")&
                  myid,trim(caller),trim(message)
            write(*,'(A)') trim(string_message)
         endif
         call MPI_Barrier(mpi_comm_global,mpierr)
      enddo

      if(n_procs > 1) then
         call MPI_Abort(mpi_comm_global,0,mpierr)
      endif
   else
      write(string_message,"(1X,'*** Proc  N/A',' in ',A,': ',A)")&
            trim(caller),trim(message)
      write(*,'(A)') trim(string_message)
   end if

   stop

end subroutine

!>
!! This routine sets the real hamiltonian matrix.
!!
subroutine elsi_set_real_hamiltonian(H_in)

   implicit none

   real(kind=dp), target :: H_in(n_l_rows,n_l_cols) !< Real Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_real_hamiltonian"

   select case (method)
      case (ELPA)
         ham_real => H_in
      case (LIBOMM)
         call m_register_pdbc(ham_omm,H_in,sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the complex hamiltonian matrix.
!!
subroutine elsi_set_complex_hamiltonian(H_in)

   implicit none

   complex(kind=dp), target :: H_in(n_l_rows,n_l_cols) !< Complex Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_complex_hamiltonian"

   select case (method)
      case (ELPA)
         ham_complex => H_in
      case (LIBOMM)
         call m_register_pdbc(ham_omm,H_in,sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the sparse real Hamiltonian matrix.
!!
subroutine elsi_set_sparse_real_hamiltonian(H_in)

   implicit none

   real(kind=dp), target :: H_in(nnz_l_pexsi) !< Sparse real Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_sparse_real_hamiltonian"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ham_real_ccs => H_in
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         ham_real_ccs => H_in
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the sparse complex Hamiltonian matrix.
!!
subroutine elsi_set_sparse_complex_hamiltonian(H_in)

   implicit none

   complex(kind=dp), target :: H_in(nnz_l_pexsi) !< Sparse complex Hamiltonian matrix

   character*40, parameter :: caller = "elsi_set_sparse_complex_hamiltonian"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ham_complex_ccs => H_in         
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the real overlap matrix.
!!
subroutine elsi_set_real_overlap(S_in)

   implicit none

   real(kind=dp), target :: S_in(n_l_rows,n_l_cols) !< Real overlap matrix

   character*40, parameter :: caller = "elsi_set_real_overlap"

   select case (method)
      case (ELPA)
         ovlp_real => S_in
      case (LIBOMM)
         call m_register_pdbc(ovlp_omm,S_in,sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         ! Nothing to be done here
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the complex overlap matrix.
!!
subroutine elsi_set_complex_overlap(S_in)

   implicit none

   complex(kind=dp), target :: S_in(n_l_rows,n_l_cols) !< Complex overlap matrix

   character*40, parameter :: caller = "elsi_set_complex_overlap"

   select case (method)
      case (ELPA)
         ovlp_complex => S_in
      case (LIBOMM)
         call m_register_pdbc(ovlp_omm,S_in,sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the sparse real overlap matrix.
!!
subroutine elsi_set_sparse_real_overlap(S_in)

   implicit none

   real(kind=dp), target :: S_in(nnz_l_pexsi) !< Sparse real overlap matrix

   character*40, parameter :: caller = "elsi_set_sparse_real_overlap"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ovlp_real_ccs => S_in
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         ovlp_real_ccs => S_in
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the sparse complex overlap matrix.
!!
subroutine elsi_set_sparse_complex_overlap(S_in)

   implicit none

   complex(kind=dp), target :: S_in(nnz_l_pexsi) !< Sparse complex overlap matrix

   character*40, parameter :: caller = "elsi_set_sparse_complex_overlap"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ovlp_complex_ccs => S_in
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the eigenvalues.
!!
subroutine elsi_set_eigenvalue(e_val_in)

   implicit none

   real(kind=dp), target :: e_val_in(n_g_size) !< Eigenvalues

   character*40, parameter :: caller = "elsi_set_eigenvalue"

   select case (method)
      case (ELPA)
         eval => e_val_in
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         eval => e_val_in
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the real eigenvectors.
!!
subroutine elsi_set_real_eigenvector(e_vec_in)

   implicit none
   
   real(kind=dp), target :: e_vec_in(n_l_rows,n_l_cols) !< Real eigenvectors

   character*40, parameter :: caller = "elsi_set_real_eigenvector"

   select case (method)
      case (ELPA)
         evec_real => e_vec_in
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         evec_real => e_vec_in
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the complex eigenvectors.
!!
subroutine elsi_set_complex_eigenvector(e_vec_in)

   implicit none

   complex(kind=dp), target :: e_vec_in(n_l_rows,n_l_cols) !< Complex eigenvectors

   character*40, parameter :: caller = "elsi_set_complex_eigenvector"

   select case (method)
      case (ELPA)
         evec_complex => e_vec_in
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the density matrix.
!!
subroutine elsi_set_density_matrix(D_in)

   implicit none

   real(kind=dp), target :: D_in(n_l_rows,n_l_cols) !< Density matrix

   character*40, parameter :: caller = "elsi_set_density_matrix"

   select case (method)
      case (ELPA)
         den_mat => D_in
      case (LIBOMM)
         call m_register_pdbc(den_mat_omm,D_in,sc_desc)
      case (PEXSI)
         ! Nothing to be done here
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the sparse density matrix.
!!
subroutine elsi_set_sparse_density_matrix(D_in)

   implicit none

   real(kind=dp), target :: D_in(nnz_l_pexsi) !< Sparse density matrix

   character*40, parameter :: caller = "elsi_set_sparse_density_matrix"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         den_mat_ccs => D_in
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         call elsi_stop(" SIPS not yet implemented. Exiting...",caller)
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the row index.
!!
subroutine elsi_set_row_ind(row_ind_in)

   implicit none

   integer, target :: row_ind_in(nnz_l_pexsi) !< Row index

   character*40, parameter :: caller = "elsi_set_row_ind"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         row_ind_ccs => row_ind_in
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         row_ind_ccs => row_ind_in
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine sets the column pointer.
!!
subroutine elsi_set_col_ptr(col_ptr_in)

   implicit none

   integer, target :: col_ptr_in(n_l_cols_pexsi+1) !< Column pointer

   character*40, parameter :: caller = "elsi_set_cplumn_pointer"

   select case (method)
      case (ELPA)
         ! Nothing to be done here
      case (LIBOMM)
         ! Nothing to be done here
      case (PEXSI)
         col_ptr_ccs => col_ptr_in
      case (CHESS)
         call elsi_stop(" CHESS not yet implemented. Exiting...",caller)
      case (SIPS)
         col_ptr_ccs => col_ptr_in
      case DEFAULT
         call elsi_stop(" No supported solver has been chosen."//&
                        " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                        " Exiting...",caller)
   end select

end subroutine

!>
!! This routine frees memory.
!!
subroutine elsi_cleanup()

   implicit none

   ! Nullify pointers
   if(associated(ham_real))         nullify(ham_real)
   if(associated(ham_complex))      nullify(ham_complex)
   if(associated(ovlp_real))        nullify(ovlp_real)
   if(associated(ovlp_complex))     nullify(ovlp_complex)
   if(associated(evec_real))        nullify(evec_real)
   if(associated(evec_complex))     nullify(evec_complex)
   if(associated(eval))             nullify(eval)
   if(associated(den_mat))          nullify(den_mat)
   if(associated(ham_real_ccs))     nullify(ham_real_ccs)
   if(associated(ham_complex_ccs))  nullify(ham_complex_ccs)
   if(associated(ovlp_real_ccs))    nullify(ovlp_real_ccs)
   if(associated(ovlp_complex_ccs)) nullify(ovlp_complex_ccs)
   if(associated(den_mat_ccs))      nullify(den_mat_ccs)
   if(associated(row_ind_ccs))      nullify(row_ind_ccs)
   if(associated(col_ptr_ccs))      nullify(col_ptr_ccs)

   ! ELPA
   if(allocated(ham_real_elpa))     deallocate(ham_real_elpa)
   if(allocated(ham_complex_elpa))  deallocate(ham_complex_elpa)
   if(allocated(ovlp_real_elpa))    deallocate(ovlp_real_elpa)
   if(allocated(ovlp_complex_elpa)) deallocate(ovlp_complex_elpa)
   if(allocated(evec_real_elpa))    deallocate(evec_real_elpa)
   if(allocated(evec_complex_elpa)) deallocate(evec_complex_elpa)
   if(allocated(eval_elpa))         deallocate(eval_elpa)
   if(allocated(den_mat_elpa))      deallocate(den_mat_elpa)
   if(allocated(occ_elpa))          deallocate(occ_elpa)

   ! libOMM
   if(ham_omm%is_initialized)       call m_deallocate(ham_omm)
   if(ovlp_omm%is_initialized)      call m_deallocate(ovlp_omm)
   if(den_mat_omm%is_initialized)   call m_deallocate(den_mat_omm)
   if(coeff_omm%is_initialized)     call m_deallocate(coeff_omm)
   if(t_den_mat_omm%is_initialized) call m_deallocate(t_den_mat_omm)

   ! PEXSI
   if(allocated(ham_real_pexsi))    deallocate(ham_real_pexsi)
   if(allocated(ovlp_real_pexsi))   deallocate(ovlp_real_pexsi)
   if(allocated(den_mat_pexsi))     deallocate(den_mat_pexsi)
   if(allocated(e_den_mat_pexsi))   deallocate(e_den_mat_pexsi)
   if(allocated(f_den_mat_pexsi))   deallocate(f_den_mat_pexsi)
   if(allocated(row_ind_pexsi))     deallocate(row_ind_pexsi)
   if(allocated(col_ptr_pexsi))     deallocate(col_ptr_pexsi)

   ! SIPs
   if(allocated(inertias))          deallocate(inertias)
   if(allocated(shifts))            deallocate(shifts)
   if(allocated(slices))            deallocate(slices)

   if(allocated(local_row))         deallocate(local_row)
   if(allocated(local_col))         deallocate(local_col)

   ! Finalize PEXSI plan
   if(method == PEXSI) then
      call f_ppexsi_plan_finalize(pexsi_plan,pexsi_info)
   endif

   ! Finalize QETSC
   if(method == SIPS) then
      call finalize_sips()
   endif

end subroutine

!>
!! This routine guarantees that there are no mutually conflicting parameters.
!!
subroutine elsi_check()

   implicit none

   character*40, parameter :: caller = "elsi_check"

   ! General check of method, parallelism, storage
   if(method < 0 .or. method > 5) then
      call elsi_stop(" No supported solver has been chosen."//&
                     " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                     " Exiting...",caller)
   endif

   if(parallelism < 0 .or. parallelism > 1) then
      call elsi_stop(" No supported parallelism has been chosen."//&
                     " Please choose SINGLE_PROC or MULTI_PROC parallelism."//&
                     " Exiting...",caller)
   endif

   if(storage < 0 .or. storage > 1) then
      call elsi_stop(" No supported format has been chosen."//&
                     " Please choose BLACS_DENSE or PEXSI_CSC format."//&
                     " Exiting...",caller)
   endif

   ! Specific check for each solver
   if(method == AUTO) then
      call elsi_stop(" AUTO not yet available."//&
                     " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                     " Exiting...",caller)

   else if(method == ELPA) then
      if(parallelism == MULTI_PROC) then
         if(.not.mpi_is_setup) then
            call elsi_stop(" MULTI_PROC parallelism requires MPI being"//&
                           " set up before calling the solver."//&
                           " Exiting...",caller)
         endif

         if(.not.blacs_is_setup) then
            call elsi_stop(" MULTI_PROC parallelism requires BLACS being"//&
                           " set up before calling the solver."//&
                           " Exiting...",caller)
         endif
      endif

      if(storage /= BLACS_DENSE) then
         call elsi_stop(" ELPA has been chosen as the solver."//&
                        " Please choose BLACS_DENSE matrix format."//&
                        " Exiting...",caller)
      endif

   else if(method == LIBOMM) then
      if(parallelism == MULTI_PROC) then
         if(.not.mpi_is_setup) then
            call elsi_stop(" MULTI_PROC parallelism requires MPI being"//&
                           " set up before calling the solver."//&
                           " Exiting...",caller)
         endif

         if(.not.blacs_is_setup) then
            call elsi_stop(" MULTI_PROC parallelism requires BLACS being"//&
                           " set up before calling the solver."//&
                           " Exiting...",caller)
         endif
      else
         call elsi_stop(" libOMM has been chosen as the solver."//&
                        " Please choose MULTI_PROC parallelism."//&
                        " Exiting...",caller)
      endif

      if(storage /= BLACS_DENSE) then
         call elsi_stop(" libOMM has been chosen as the solver."//&
                        " Please choose BLACS_DENSE matrix format."//&
                        " Exiting...",caller)
      endif

   else if(method == PEXSI) then
      if(parallelism == MULTI_PROC) then
         if(.not.mpi_is_setup) then
            call elsi_stop(" MULTI_PROC parallelism requires MPI being"//&
                           " set up before calling the solver."//&
                           " Exiting...",caller)
         endif
      else
         call elsi_stop(" PEXSI has been chosen as the solver."//&
                        " Please choose MULTI_PROC parallelism."//&
                        " Exiting...",caller)
      endif

      if(storage == BLACS_DENSE) then
         if(.not.blacs_is_setup) then
            call elsi_stop(" The BLACS_DENSE format has been chosen."//&
                           " Please set up BLACS before calling the"//&
                           " Solver. Exiting...",caller)
         endif
      else
         if(.not.sparsity_pattern_ready) then
            call elsi_stop(" The PEXSI_CSC format has been chosen."//&
                           " Please set the sparsity pattern before"//&
                           " calling the solver. Exiting...",caller)
         endif
      endif

      if(n_g_size < n_p_per_pole_pexsi) then
         call elsi_stop(" PEXSI has been chosen as the solver."//&
                        " The size of matrix is too small for this"//&
                        " number of processes. Exiting...",caller)
      endif

   else if(method == CHESS) then
      call elsi_stop(" CHESS not yet available."//&
                     " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                     " Exiting...",caller)

   else if(method == SIPS) then
      call elsi_statement_print("  ATTENTION! SIPS is EXPERIMENTAL.")

      if(n_g_size < n_procs) then
         call elsi_stop(" SIPs has been chosen as the solver."//&
                        " The size of matrix is too small for this"//&
                        " number of processes. Exiting...",caller)
      endif

   else
      call elsi_stop(" No supported solver has been chosen."//&
                     " Please choose ELPA, LIBOMM, or PEXSI solver."//&
                     " Exiting...",caller)
   endif

end subroutine

!> 
!! This routine computes global row index based on local row index.
!!
subroutine elsi_get_global_row(global_idx,local_idx)

   implicit none

   integer, intent(in)  :: local_idx  !< Local index
   integer, intent(out) :: global_idx !< Global index

   integer :: block !< Local block
   integer :: idx !< Local index in block

   block = (local_idx-1)/n_b_rows
   idx = local_idx-block*n_b_rows

   global_idx = my_p_row*n_b_rows+block*n_b_rows*n_p_rows+idx

end subroutine

!> 
!! This routine computes global column index based on local column index.
!!
subroutine elsi_get_global_col(global_idx,local_idx)

   implicit none

   integer, intent(in)  :: local_idx  !< Local index
   integer, intent(out) :: global_idx !< Global index

   integer :: block !< Local block
   integer :: idx !< Local index in block

   block = (local_idx-1)/n_b_cols
   idx = local_idx-block*n_b_cols

   global_idx = my_p_col*n_b_cols+block*n_b_cols*n_p_cols+idx

end subroutine

!>
!! This routine counts the local number of non_zero elements.
!!
subroutine elsi_get_local_nnz(matrix,n_rows,n_cols,nnz)

   implicit none

   real(kind=dp), intent(in)  :: matrix(n_rows,n_cols) !< Local matrix
   integer,       intent(in)  :: n_rows                !< Local rows
   integer,       intent(in)  :: n_cols                !< Local cols
   integer,       intent(out) :: nnz                   !< Number of non-zero

   integer :: i_row !< Row counter
   integer :: i_col !< Column counter
   integer :: i_nz !< Non-zero element counter

   nnz = 0

   do i_col = 1,n_cols
      do i_row = 1,n_rows
         if(abs(matrix(i_row,i_col)) > zero_threshold) then
            nnz = nnz+1
         endif
      enddo
   enddo

end subroutine

end module ELSI_UTILS
