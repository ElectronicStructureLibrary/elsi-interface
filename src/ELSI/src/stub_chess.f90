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
!! This module is only compiled when the actual CheSS is not available,
!! to make the CheSS part of ELSI compile.
!!
module FOE_BASE

   use ELSI_PRECISION, only: r8,i4

   implicit none

   public :: foe_data_deallocate
   public :: foe_data_get_real

   type, public :: foe_data
      integer(kind=i4) :: dummy
   end type

contains

subroutine foe_data_deallocate(foe_obj)

   implicit none

   type(foe_data) :: foe_obj

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

function foe_data_get_real(foe_obj,fieldname,ind) result(val)

   implicit none

   type(foe_data)    :: foe_obj
   character(len=*)  :: fieldname
   integer, optional :: ind
   real(kind=r8)     :: val

   val = 0.0_r8

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end function

end module

module FOE_COMMON

   use ELSI_PRECISION, only: r8,i4
   use FOE_BASE, only: foe_data
   
   implicit none

   public :: init_foe

contains

subroutine init_foe(iproc,nproc,nspin,charge,foe_obj,ef,tmprtr,&
              evbounds_nsatur,evboundsshrink_nsatur,evlow,evhigh,&
              fscale,ef_interpol_det,ef_interpol_chargediff,&
              fscale_lowerbound,fscale_upperbound,eval_multiplicator,&
              npl_min,npl_max,npl_stride,betax,ntemp,accuracy_function,&
              accuracy_penalty)

   implicit none

   integer(kind=i4)           :: iproc
   integer(kind=i4)           :: nproc
   integer(kind=i4)           :: nspin
   real(kind=r8)              :: charge(nspin)
   type(foe_data)             :: foe_obj
   integer(kind=i4), optional :: evbounds_nsatur
   integer(kind=i4), optional :: evboundsshrink_nsatur
   real(kind=r8),    optional :: ef
   real(kind=r8),    optional :: evlow
   real(kind=r8),    optional :: evhigh
   real(kind=r8),    optional :: fscale
   real(kind=r8),    optional :: ef_interpol_det
   real(kind=r8),    optional :: ef_interpol_chargediff
   real(kind=r8),    optional :: fscale_lowerbound
   real(kind=r8),    optional :: fscale_upperbound
   real(kind=r8),    optional :: tmprtr
   real(kind=r8),    optional :: eval_multiplicator
   integer(kind=i4), optional :: npl_min
   integer(kind=i4), optional :: npl_max
   integer(kind=i4), optional :: npl_stride
   real(kind=r8),    optional :: betax
   integer(kind=i4), optional :: ntemp
   real(kind=r8),    optional :: accuracy_function
   real(kind=r8),    optional :: accuracy_penalty

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

end module

module SPARSEMATRIX_BASE

   use ELSI_PRECISION, only: r8,i4
   
   implicit none

   public :: sparsematrix_init_errors
   public :: sparsematrix_initialize_timing_categories
   public :: deallocate_sparse_matrix
   public :: deallocate_matrices

   integer(kind=i4), parameter, public :: SPARSE_TASKGROUP = 50

   type, public :: matrices
      real(kind=r8), pointer :: matrix_compr(:)
   end type

   type, public :: sparse_matrix
      integer(kind=i4) :: nspin = 1
   end type

contains

subroutine sparsematrix_init_errors()

   implicit none

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sparsematrix_initialize_timing_categories()

   implicit none

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine deallocate_sparse_matrix(sparsemat)

   implicit none

   type(sparse_matrix) :: sparsemat

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine deallocate_matrices(mat)

   implicit none

   type(matrices) :: mat

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

end module

module SPARSEMATRIX_HIGHLEVEL

   use ELSI_PRECISION, only: r8,i4
   use FOE_BASE, only: foe_data
   use SPARSEMATRIX_BASE, only: matrices,sparse_matrix
   
   implicit none

   public :: matrices_init
   public :: matrix_fermi_operator_expansion
   public :: sparse_matrix_init_from_data_ccs

contains

subroutine matrices_init(smat,mat,matsize)

   implicit none

   type(sparse_matrix) :: smat
   type(matrices)      :: mat
   integer, optional   :: matsize

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine matrix_fermi_operator_expansion(iproc,nproc,comm,foe_obj,ice_obj,&
              smat_s,smat_h,smat_k,overlap,ham,overlap_minus_one_half,kernel,&
              ebs,calculate_minusonehalf,foe_verbosity,symmetrize_kernel,&
              calculate_energy_density_kernel,calculate_spin_channels,&
              energy_kernel,inversion_method,pexsi_np_sym_fact)

   implicit none

   integer(kind=i4)           :: iproc
   integer(kind=i4)           :: nproc
   integer(kind=i4)           :: comm
   type(foe_data)             :: foe_obj
   type(foe_data)             :: ice_obj
   type(sparse_matrix)        :: smat_s
   type(sparse_matrix)        :: smat_h
   type(sparse_matrix)        :: smat_k
   type(matrices)             :: overlap
   type(matrices)             :: ham
   type(matrices)             :: overlap_minus_one_half(1)
   type(matrices)             :: kernel
   real(kind=r8)              :: ebs
   logical,          optional :: calculate_minusonehalf
   logical,          optional :: symmetrize_kernel
   logical,          optional :: calculate_energy_density_kernel
   integer(kind=i4), optional :: foe_verbosity
   type(matrices),   optional :: energy_kernel
   logical,          optional :: calculate_spin_channels(smat_k%nspin)
   character(len=*), optional :: inversion_method
   integer(kind=i4), optional :: pexsi_np_sym_fact

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine sparse_matrix_init_from_data_ccs(iproc,nproc,comm,nfvctr,&
              nvctr,row_ind,col_ptr,smat,init_matmul,nvctr_mult,&
              row_ind_mult,col_ptr_mult)

   implicit none

   integer(kind=i4)           :: iproc
   integer(kind=i4)           :: nproc
   integer(kind=i4)           :: comm
   integer(kind=i4)           :: nfvctr
   integer(kind=i4)           :: nvctr
   integer(kind=i4)           :: row_ind(nvctr)
   integer(kind=i4)           :: col_ptr(nfvctr)
   type(sparse_matrix)        :: smat
   logical,          optional :: init_matmul
   integer(kind=i4), optional :: nvctr_mult
   integer(kind=i4), optional :: row_ind_mult(:)
   integer(kind=i4), optional :: col_ptr_mult(:)

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

end module

module SPARSEMATRIX_INIT

   use ELSI_PRECISION, only: r8,i4
   use SPARSEMATRIX_BASE, only: sparse_matrix
   
   implicit none

   public :: init_matrix_taskgroups_wrapper

contains

subroutine init_matrix_taskgroups_wrapper(iproc,nproc,comm,&
              enable_matrix_taskgroups,nmat,smat,ind_minmax)

   implicit none

   integer(kind=i4)    :: iproc
   integer(kind=i4)    :: nproc
   integer(kind=i4)    :: comm
   logical             :: enable_matrix_taskgroups
   integer(kind=i4)    :: nmat
   type(sparse_matrix) :: smat(nmat)
   integer, optional   :: ind_minmax(2,nmat)

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

end module

subroutine f_lib_initialize()

   implicit none

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine

subroutine f_lib_finalize()

   implicit none

   write(*,"(A)") " A CheSS stub routine was called. Check ELSI installation."
   write(*,"(A)") " Exiting..."
   stop

end subroutine
