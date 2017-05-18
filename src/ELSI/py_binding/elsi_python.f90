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
!! This file contains the Python interfaces of ELSI.
!!

!  Compiling and Using the Python Interface
!  - You will need f2py, packaged as part of numpy as of version XXX.
!  - You will also need an MPI-aware Python shell.  The mpi4py module seems to be one 
!    of the more popular choices.  Note that you will need to have Python be MPI-aware 
!    even if you are running ELPA in single-proc mode where you never need to set up
!    the MPI communicators (or BLACS desciptors); you must import the MPI module (for
!    mpi4py, via "from mpi4py import MPI") anyway!
!  - f2py creates a shared library, so all source (including dependences!) must
!    be position-independent code.  (This in practice means compiled with -fPIC.)
!    This will introduce a performance cost.
!  - f2py automatically converts all subroutines to lower-case, and Python is a case
!    sensitive language.
!  - Python uses "our" standard lengths for primitives:
!    - 64-bit floating point (numpy.float64)
!    - 32-bit integers (numpy.int32)
!    - False = 0 and True = 1
!    but, since Python loves implicit type casting, manually specifying this for arrays
!    is a good idea.
!
!  Implementing the Python Interface
!  - f2py does not support (Fortran) pointers or derived types.  This does not affect 
!    us right now (other than in the tomato interface, but I have a temporary hack around
!    that) but it's good to keep in mind for future work on the interface.
!  - f2py assumes that by default that pass-by-reference should be done.  Which
!    is normal behavior when writing wrappers.  However, if you then specify the 
!    "value" keyword, it will honor that in its own way and will make a copy of 
!    what was passed.  Which was the address.  So you'll essentially get a random
!    number that is a pain-in-the-ass to debug.  What I'm saying is don't
!    use the "value" keyword.
!  - f2py does not support assumed size/shape arrays.  It also does not support array
!    dimensions taken from modules.  This means that array dimensions must be coded 
!    explicitly as arguments, but this is easy to implement (see next.)
!  - f2py will create a Pythonic subroutine.  That means that it will hide arguments
!    that it has determined are array dimensions and output intent(out) variables as a
!    tuple, implying a memory copy and a fork in the argument list (i.e. the Python
!    function will have a different number of arguments as the C/Fortran functions.)  
!    Accordingly, a subroutine that in ordinary Fortran would have the header
!      subroutine A(I,O,a_number)
!        real*8,  intent(in)  :: I(:)
!        real*8,  intent(out) :: O(:)
!        real*8,  intent(out) :: a_number
!        ...
!    should have a header in the f2py wrapper subroutine of
!      subroutine A(O,I,a_number,size_I,size_O)
!        integer, intent(in)    :: size_I
!        integer, intent(in)    :: size_O
!        real*8,  intent(inout) :: I(size_I)
!        real*8,  intent(inout) :: O(size_O)
!        real*8,  intent(out)   :: a_number
!        ...
!    (To make sure f2py has really hidden the array indices properly, look at
!    the signature file or the info string of the resulting Python function.)
!    The compiled subroutine would be called from within Python as
!      a_number = A(O, I)
!    where O and I are arrays of size size_I and size_O, respectively.
!
!  TODO
!  - Finalize naming scheme

subroutine elsi_init_py_wrapper(solver,parallel_mode,matrix_format,&
                               matrix_size,n_electrons_in,n_states_in)

   use ELSI, only: elsi_init

   implicit none

   ! f2py arguments
   integer, intent(in) :: solver
   integer, intent(in) :: parallel_mode
   integer, intent(in) :: matrix_format
   integer, intent(in) :: matrix_size
   real*8,  intent(in) :: n_electrons_in
   integer, intent(in) :: n_states_in

   call elsi_init(solver,parallel_mode,matrix_format,&
                  matrix_size,n_electrons_in,n_states_in)

end subroutine

subroutine elsi_set_mpi_py_wrapper(mpi_comm_global_in)

   use ELSI, only: elsi_set_mpi

   implicit none

   ! f2py arguments
   integer, intent(in) :: mpi_comm_global_in

   call elsi_set_mpi(mpi_comm_global_in)

end subroutine

subroutine elsi_set_blacs_py_wrapper(icontext,block_size)

   use ELSI, only: elsi_set_blacs

   implicit none

   ! f2py arguments
   integer, intent(in) :: icontext
   integer, intent(in) :: block_size

   call elsi_set_blacs(icontext,block_size)

end subroutine

subroutine elsi_set_csc_py_wrapper(nnz_g_in,nnz_l_in,nnz_l_cols_in,row_ind_in,&
                                  col_ptr_in,size_row_ind_in,size_col_ptr_in)

   use ELSI, only: elsi_set_csc

   implicit none

   ! Array dimensions hidden by f2py
   integer, intent(in) :: size_row_ind_in
   integer, intent(in) :: size_col_ptr_in
   ! f2py arguments
   integer, intent(in) :: nnz_g_in
   integer, intent(in) :: nnz_l_in
   integer, intent(in) :: nnz_l_cols_in
   integer, intent(in) :: row_ind_in(size_row_ind_in)
   integer, intent(in) :: col_ptr_in(size_col_ptr_in)

   call elsi_set_csc(nnz_g_in,nnz_l_in,nnz_l_cols_in,row_ind_in,col_ptr_in)

end subroutine

subroutine elsi_finalize_py_wrapper()

   use ELSI, only: elsi_finalize

   implicit none

   call elsi_finalize()

end subroutine

subroutine elsi_customize_py_wrapper(print_detail,unit_overlap,numerical_zero,&
                                    no_check_singularity,singularity_threshold,&
                                    stop_singularity)

   use ELSI, only: elsi_customize

   implicit none

   ! f2py arguments
   logical, intent(in) :: print_detail           ! Default: .false.
   logical, intent(in) :: unit_overlap           ! Default: .false.
   real*8,  intent(in) :: numerical_zero         ! Default:  1d-13
   logical, intent(in) :: no_check_singularity   ! Default:  .false.
   real*8,  intent(in) :: singularity_threshold  ! Default: 1d-5
   logical, intent(in) :: stop_singularity ! Default: .false.

!   logical :: print_detail_f
!   logical :: unit_overlap_f
!   logical :: no_check_singularity_f
!   logical :: stop_singularity_f

!   if(print_detail == 0) then
!      print_detail_f = .false.
!   else
!      print_detail_f = .true.
!   endif

!   if(unit_overlap == 0) then
!      unit_overlap_f = .false.
!   else
!      unit_overlap_f = .true.
!   endif

!   if(no_check_singularity == 0) then
!      no_check_singularity_f = .false.
!   else
!      no_check_singularity_f = .true.
!   endif

!   if(stop_singularity == 0) then
!      stop_singularity_f = .false.
!   else
!      stop_singularity_f = .true.
!   endif

   call elsi_customize(print_detail,unit_overlap,numerical_zero,&
                       no_check_singularity,singularity_threshold,&
                       stop_singularity)

end subroutine

subroutine elsi_customize_mu_py_wrapper(broadening_scheme,broadening_width,&
                                       mu_accuracy,mu_max_steps,spin_degeneracy)

   use ELSI, only: elsi_customize_mu

   implicit none

   ! f2py arguments
   integer, intent(in) :: broadening_scheme
   real*8,  intent(in) :: broadening_width
   real*8,  intent(in) :: mu_accuracy
   integer, intent(in) :: mu_max_steps
   real*8,  intent(in) :: spin_degeneracy

   call elsi_customize_mu(broadening_scheme,broadening_width,&
                          mu_accuracy,mu_max_steps,spin_degeneracy)

end subroutine

subroutine elsi_customize_omm_py_wrapper(n_elpa_steps_omm,omm_method,eigen_shift,&
                                        omm_tolerance,use_pspblas,omm_output)

   use ELSI, only: elsi_customize_omm

   implicit none

   ! f2py arguments
   integer, intent(in) :: n_elpa_steps_omm
   integer, intent(in) :: omm_method
   real*8,  intent(in) :: eigen_shift
   real*8,  intent(in) :: omm_tolerance
   logical, intent(in) :: use_pspblas
   logical, intent(in) :: omm_output

!   logical :: use_pspblas_f
!   logical :: omm_output_f

!   if(use_pspblas == 0) then
!      use_pspblas_f = .false.
!   else
!      use_pspblas_f = .true.
!   endif

!   if(omm_output == 0) then
!      omm_output_f = .false.
!   else
!      omm_output_f = .true.
!   endif

   call elsi_customize_omm(n_elpa_steps_omm,omm_method,eigen_shift,&
                           omm_tolerance,use_pspblas,omm_output)

end subroutine

subroutine elsi_customize_pexsi_py_wrapper(temperature,gap,delta_E,n_poles,&
                                          n_procs_per_pole,max_iteration,mu_min,&
                                          mu_max,mu0,mu_inertia_tolerance,&
                                          mu_inertia_expansion,mu_safeguard,&
                                          n_electron_accuracy,matrix_type,&
                                          is_symbolic_factorize,ordering,&
                                          np_symbolic_factorize,verbosity)

   use ELSI, only: elsi_customize_pexsi

   implicit none

   ! f2py arguments
   real*8,  intent(in) :: temperature
   real*8,  intent(in) :: gap
   real*8,  intent(in) :: delta_E
   integer, intent(in) :: n_poles
   integer, intent(in) :: n_procs_per_pole
   integer, intent(in) :: max_iteration
   real*8,  intent(in) :: mu_min
   real*8,  intent(in) :: mu_max
   real*8,  intent(in) :: mu0
   real*8,  intent(in) :: mu_inertia_tolerance
   real*8,  intent(in) :: mu_inertia_expansion
   real*8,  intent(in) :: mu_safeguard
   real*8,  intent(in) :: n_electron_accuracy
   integer, intent(in) :: matrix_type
   integer, intent(in) :: is_symbolic_factorize
   integer, intent(in) :: ordering
   integer, intent(in) :: np_symbolic_factorize
   integer, intent(in) :: verbosity

   call elsi_customize_pexsi(temperature,gap,delta_E,n_poles,n_procs_per_pole,&
           max_iteration,mu_min,mu_max,mu0,mu_inertia_tolerance,mu_inertia_expansion,&
           mu_safeguard,n_electron_accuracy,matrix_type,is_symbolic_factorize,&
           ordering,np_symbolic_factorize,verbosity)

end subroutine

subroutine elsi_customize_elpa_py_wrapper(elpa_solver)

   use ELSI, only: elsi_customize_elpa

   implicit none

   ! f2py arguments
   integer, intent(in) :: elpa_solver

   call elsi_customize_elpa(elpa_solver)

end subroutine

subroutine elsi_ev_real_py_wrapper(H_in,S_in,e_val_out,e_vec_out,n_l_rows,n_g_size,n_cols_H_in,n_cols_S_in,n_cols_e_vec_out)

   use ELSI, only: elsi_ev_real

   implicit none

   ! Array dimensions hidden by f2py
   integer, intent(in)    :: n_l_rows
   integer, intent(in)    :: n_g_size
   integer, intent(in)    :: n_cols_H_in
   integer, intent(in)    :: n_cols_S_in
   integer, intent(in)    :: n_cols_e_vec_out
   ! f2py arguments
   real*8,  intent(in)    :: H_in(0:n_l_rows-1,0:n_cols_H_in-1)
   real*8,  intent(in)    :: S_in(0:n_l_rows-1,0:n_cols_S_in-1)
   real*8,  intent(inout) :: e_val_out(0:n_g_size-1)
   real*8,  intent(inout) :: e_vec_out(0:n_l_rows-1,0:n_cols_e_vec_out-1)

   call elsi_ev_real(H_in,S_in,e_val_out,e_vec_out)

end subroutine

subroutine elsi_ev_complex_py_wrapper(H_in,S_in,e_val_out,e_vec_out,n_l_rows,n_g_size,n_cols_H_in,n_cols_S_in,n_cols_e_vec_out)

   use ELSI, only: elsi_ev_complex

   implicit none

   ! Array dimensions hidden by f2py
   integer,    intent(in)    :: n_l_rows
   integer,    intent(in)    :: n_g_size
   integer,    intent(in)    :: n_cols_H_in
   integer,    intent(in)    :: n_cols_S_in
   integer,    intent(in)    :: n_cols_e_vec_out
   ! f2py arguments
   complex*16, intent(in)    :: H_in(0:n_l_rows-1,0:n_cols_H_in-1)
   complex*16, intent(in)    :: S_in(0:n_l_rows-1,0:n_cols_S_in-1)
   real*8,     intent(inout) :: e_val_out(0:n_g_size-1)
   complex*16, intent(inout) :: e_vec_out(0:n_l_rows-1,0:n_cols_e_vec_out-1)

   call elsi_ev_complex(H_in,S_in,e_val_out,e_vec_out)

end subroutine

function elsi_dm_real_py_wrapper(H_in,S_in,D_out,n_l_rows, n_cols_H_in, n_cols_S_in, n_cols_D_out) result(energy_out)

   use ELSI, only: elsi_dm_real

   implicit none

   ! Array dimensions hidden by f2py
   integer, intent(in)    :: n_l_rows
   integer, intent(in)    :: n_cols_H_in
   integer, intent(in)    :: n_cols_S_in
   integer, intent(in)    :: n_cols_D_out
   ! f2py arguments
   real*8,  intent(in)    :: H_in(n_l_rows,n_cols_H_in)
   real*8,  intent(in)    :: S_in(n_l_rows,n_cols_S_in)
   real*8,  intent(inout) :: D_out(n_l_rows,n_cols_D_out)
   real*8                 :: energy_out

   call elsi_dm_real(H_in,S_in,D_out,energy_out)

end function

function elsi_dm_real_sparse_py_wrapper(H_in,S_in,D_out,size_H_in,size_S_in,size_D_out) result(energy_out)

   use ELSI, only: elsi_dm_real_sparse

   implicit none

   ! Array dimensions hidden by f2py
   integer, intent(in)    :: size_H_in
   integer, intent(in)    :: size_S_in
   integer, intent(in)    :: size_D_out
   ! f2py arguments
   real*8,  intent(in)    :: H_in(size_H_in)
   real*8,  intent(in)    :: S_in(size_S_in)
   real*8,  intent(inout) :: D_out(size_D_out)
   real*8                 :: energy_out

   call elsi_dm_real_sparse(H_in,S_in,D_out,energy_out)

end function

subroutine elsi_collect_pexsi_py_wrapper(mu_out,edm_out,fdm_out, size_edm_out, size_fdm_out)

   use ELSI, only: elsi_collect_pexsi

   implicit none

   ! Array dimensions hidden by f2py
   integer, intent(in)    :: size_edm_out
   integer, intent(in)    :: size_fdm_out
   ! f2py arguments
   real*8,  intent(inout) :: mu_out
   real*8,  intent(inout) :: edm_out(size_edm_out)
   real*8,  intent(inout) :: fdm_out(size_fdm_out)

   call elsi_collect_pexsi(mu_out,edm_out,fdm_out)

end subroutine
