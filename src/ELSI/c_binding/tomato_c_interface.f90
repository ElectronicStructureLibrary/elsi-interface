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
!! This module contains the C interfaces of TOMATO. Only used in the C test
!! programs of ELSI.
!!
module TOMATO_C_INTERFACE

   use, intrinsic :: iso_c_binding
   use ELSI_C2F
   use MATRIXSWITCH

   implicit none

   private

   integer(kind=c_int) :: n_l_rows_c
   integer(kind=c_int) :: n_l_cols_c

contains

subroutine m_deallocate_c_wrapper(handle_c)&
   bind(C,name="c_m_deallocate")

   implicit none

   type(c_ptr), value :: handle_c

   type(matrix), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   call m_deallocate(handle_f)

end subroutine

subroutine ms_scalapack_setup_c_wrapper(mpi_comm,nprow,order,blk,ctxt)&
   bind(C,name="c_ms_scalapack_setup")

   implicit none

   integer(kind=c_int), value,                 intent(in) :: mpi_comm
   integer(kind=c_int), value,                 intent(in) :: nprow
   character(kind=c_char,len=1), dimension(2), intent(in) :: order
   integer(kind=c_int), value,                 intent(in) :: blk
   integer(kind=c_int), value,                 intent(in) :: ctxt

   character(2) :: order_f

   order_f = convert_c_string_to_f_string(order)

   call ms_scalapack_setup(mpi_comm,nprow,order_f,blk,icontxt=ctxt)

end subroutine

subroutine tomato_tb_c_wrapper(seed_dir,system,switch1,frac_occ,&
              n_basis_per_atom,switch2,n_basis,supercell,switch3,&
              sparsity,r_cut,n_states,gamma_only,k_point,defect,&
              perturbation,h_c,s_c,ms_storage,build_matrix)&
   bind(C,name="c_tomato_tb")

   implicit none

   character(kind=c_char,len=1), dimension(128) :: seed_dir
   character(kind=c_char,len=1), dimension(128) :: system
   integer(c_int), value,                        intent(in) :: switch1
   real(c_double)       ,                        intent(in) :: frac_occ
   integer(c_int), value,                        intent(in) :: n_basis_per_atom
   integer(c_int), value,                        intent(in) :: switch2
   integer(c_int)       ,                        intent(in) :: n_basis
   integer(c_int)       ,                        intent(in) :: supercell(3)
   integer(c_int), value,                        intent(in) :: switch3
   real(c_double)       ,                        intent(in) :: sparsity
   real(c_double)       ,                        intent(in) :: r_cut
   integer(c_int)       ,                        intent(in) :: n_states
   integer(c_int), value,                        intent(in) :: gamma_only
   real(c_double)       ,                        intent(in) :: k_point(3)
   integer(c_int), value,                        intent(in) :: defect
   real(c_double), value,                        intent(in) :: perturbation
   type(c_ptr)                                              :: h_c
   type(c_ptr)                                              :: s_c
   character(kind=c_char,len=1), dimension(6) :: ms_storage
   integer(c_int), value,                        intent(in) :: build_matrix

   type(matrix), pointer :: h_f
   type(matrix), pointer :: s_f

   logical :: switch1_f
   logical :: switch2_f
   logical :: switch3_f
   logical :: gamma_only_f
   logical :: defect_f
   logical :: build_matrix_f

   character(6)                  :: ms_storage_f
   character(len=:), allocatable :: seed_dir_f
   character(len=:), allocatable :: system_f

   allocate(h_f)
   allocate(s_f)

   switch1_f      = convert_c_int_to_f_logical(switch1)
   switch2_f      = convert_c_int_to_f_logical(switch2)
   switch3_f      = convert_c_int_to_f_logical(switch3)
   gamma_only_f   = convert_c_int_to_f_logical(gamma_only)
   defect_f       = convert_c_int_to_f_logical(defect)
   build_matrix_f = convert_c_int_to_f_logical(build_matrix)

   seed_dir_f     = convert_c_string_to_f_string(seed_dir)
   system_f       = convert_c_string_to_f_string(system)
   ms_storage_f   = convert_c_string_to_f_string(ms_storage)

   call tomato_TB(seed_dir_f,system_f,switch1_f,frac_occ,n_basis_per_atom,&
           switch2_f,n_basis,supercell,switch3_f,sparsity,r_cut,n_states,&
           gamma_only_f,k_point,defect_f,perturbation,h_f,s_f,ms_storage_f,&
           build_matrix_f)

   h_c = c_loc(h_f)
   s_c = c_loc(s_f)

   n_l_rows_c = h_f%iaux2(1)
   n_l_cols_c = h_f%iaux2(2)

end subroutine

subroutine get_mat_info_c_wrapper(handle_c,is_real,l_row,l_col)&
   bind(C,name="c_get_mat_info")

   implicit none

   type(c_ptr), value  :: handle_c
   integer(kind=c_int) :: is_real
   integer(kind=c_int) :: l_row
   integer(kind=c_int) :: l_col

   type(matrix), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   if(handle_f%is_real) then
      is_real = 1
   else
      is_real = 0
   endif

   l_row = handle_f%iaux2(1)
   l_col = handle_f%iaux2(2)

end subroutine

subroutine get_mat_real_c_wrapper(handle_c,mat_real)&
   bind(C,name="c_get_mat_real")

   implicit none

   type(c_ptr), value  :: handle_c
   real(kind=c_double) :: mat_real(n_l_rows_c,n_l_cols_c)

   type(matrix), pointer :: handle_f

   call c_f_pointer(handle_c,handle_f)

   mat_real = handle_f%dval

end subroutine

end module TOMATO_C_INTERFACE
