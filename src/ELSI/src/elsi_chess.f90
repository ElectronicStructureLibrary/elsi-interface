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
!! This module provides interfaces to CheSS.
!!
module ELSI_CHESS

   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS
   use FOE_BASE, only: foe_data_get_real
   use FOE_COMMON, only: init_foe
   use SPARSEMATRIX_BASE, only: sparsematrix_init_errors,&
                                sparsematrix_initialize_timing_categories,&
                                SPARSE_TASKGROUP
   use SPARSEMATRIX_HIGHLEVEL, only: matrices_init,&
                                     matrix_fermi_operator_expansion,&
                                     sparse_matrix_init_from_data_ccs
   use SPARSEMATRIX_INIT, only: init_matrix_taskgroups_wrapper

   implicit none
   private

   public :: elsi_init_chess
   public :: elsi_solve_evp_chess
   public :: elsi_compute_edm_chess
   public :: elsi_set_chess_default
   public :: elsi_print_chess_options

contains

!>
!! This routine initializes CheSS.
!!
subroutine elsi_init_chess(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   real(kind=r8) :: n_electron(1)

   character*40, parameter :: caller = "elsi_init_chess"

   if(elsi_h%n_elsi_calls == 1) then
      ! Initialize f_lib
      call f_lib_initialize()

      ! Initialize sparsematrix error handling and timing
      call sparsematrix_init_errors()
      call sparsematrix_initialize_timing_categories()

      ! Initialize sparse matrices
      call sparse_matrix_init_from_data_ccs(elsi_h%myid,elsi_h%n_procs,&
              elsi_h%mpi_comm,elsi_h%n_basis,elsi_h%nnz_g,elsi_h%row_ind_chess,&
              elsi_h%col_ptr_chess,elsi_h%sparse_mat(1),init_matmul=.false.)

      ! TODO: Sparsity buffer
      call sparse_matrix_init_from_data_ccs(elsi_h%myid,elsi_h%n_procs,&
              elsi_h%mpi_comm,elsi_h%n_basis,elsi_h%nnz_g,elsi_h%row_ind_chess,&
              elsi_h%col_ptr_chess,elsi_h%sparse_mat(2),init_matmul=.true.,&
              nvctr_mult=elsi_h%nnz_g,row_ind_mult=elsi_h%row_ind_chess,&
              col_ptr_mult=elsi_h%col_ptr_chess(1:elsi_h%n_basis))

      ! Initialize task groups
      call init_matrix_taskgroups_wrapper(elsi_h%myid,elsi_h%n_procs,&
              elsi_h%mpi_comm,.false.,2,elsi_h%sparse_mat)

      n_electron = elsi_h%n_electrons

      ! Initialize FOE objects
      call init_foe(elsi_h%myid,elsi_h%n_procs,1,n_electron,&
              elsi_h%foe_obj,fscale=elsi_h%erf_decay,&
              fscale_lowerbound=elsi_h%erf_decay_min,&
              fscale_upperbound=elsi_h%erf_decay_max,&
              evlow=elsi_h%ev_ham_min,evhigh=elsi_h%ev_ham_max,&
              betax=elsi_h%beta)

      call init_foe(elsi_h%myid,elsi_h%n_procs,1,n_electron,&
              elsi_h%ice_obj,evlow=elsi_h%ev_ovlp_min,&
              evhigh=elsi_h%ev_ovlp_max,betax=elsi_h%beta)

      ! Allocate CheSS matrices
      call matrices_init(elsi_h%sparse_mat(1),elsi_h%ham_chess,&
              matsize=SPARSE_TASKGROUP)
      call matrices_init(elsi_h%sparse_mat(1),elsi_h%ovlp_chess,&
              matsize=SPARSE_TASKGROUP)
      call matrices_init(elsi_h%sparse_mat(2),elsi_h%dm_chess,&
              matsize=SPARSE_TASKGROUP)
      call matrices_init(elsi_h%sparse_mat(2),elsi_h%edm_chess,&
              matsize=SPARSE_TASKGROUP)
      call matrices_init(elsi_h%sparse_mat(2),elsi_h%ovlp_inv_sqrt(1),&
              matsize=SPARSE_TASKGROUP)

      elsi_h%ovlp_chess%matrix_compr = elsi_h%ovlp_real_chess

      elsi_h%chess_started = .true.
   endif

   elsi_h%ham_chess%matrix_compr = elsi_h%ham_real_chess

end subroutine

!>
!! This routine interfaces to CheSS.
!!
subroutine elsi_solve_evp_chess(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   logical          :: calc_ovlp_inv_sqrt
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   integer(kind=i4) :: mpierr
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_solve_evp_chess"

   call elsi_get_time(elsi_h,t0)

   if(elsi_h%n_elsi_calls == 1) then
      calc_ovlp_inv_sqrt = .true.
   else
      calc_ovlp_inv_sqrt = .false.
   endif

   call elsi_statement_print("  Starting CheSS density matrix solver",elsi_h)

   call matrix_fermi_operator_expansion(elsi_h%myid,elsi_h%n_procs,&
           elsi_h%mpi_comm,elsi_h%foe_obj,elsi_h%ice_obj,&
           elsi_h%sparse_mat(1),elsi_h%sparse_mat(1),&
           elsi_h%sparse_mat(2),elsi_h%ovlp_chess,elsi_h%ham_chess,&
           elsi_h%ovlp_inv_sqrt,elsi_h%dm_chess,elsi_h%energy_hdm,&
           calculate_minusonehalf=calc_ovlp_inv_sqrt,&
           foe_verbosity=0,symmetrize_kernel=.true.,&
           calculate_energy_density_kernel=.false.,&
           energy_kernel=elsi_h%edm_chess)

   elsi_h%mu = foe_data_get_real(elsi_h%foe_obj,"ef",1)

   call MPI_Barrier(elsi_h%mpi_comm,mpierr)

   call elsi_get_time(elsi_h,t1)

   write(info_str,"('  Finished density matrix calculation')")
   call elsi_statement_print(info_str,elsi_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,elsi_h)

end subroutine

!> 
!! This routine computes the energy-weighted density matrix.
!! 
subroutine elsi_compute_edm_chess(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   character*40, parameter :: caller = "elsi_compute_edm_chess"

   ! TODO

end subroutine

!> 
!! This routine sets default CheSS parameters.
!! 
subroutine elsi_set_chess_default(elsi_h)
   
   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   character*40, parameter :: caller = "elsi_set_chess_default"

   ! Initial guess of the decay length of the error function
   elsi_h%erf_decay = 5.0e-2_r8

   ! Lower bound of the decay length
   elsi_h%erf_decay_min = 5.0e-3_r8

   ! Upper bound of the decay length
   elsi_h%erf_decay_max = 5.0e-2_r8

   ! Lower bound of the eigenvalues of H
   elsi_h%ev_ham_min = -2.0_r8

   ! Upper bound of the eigenvalues of H
   elsi_h%ev_ham_max = 2.0_r8

   ! Lower bound of the eigenvalues of S
   elsi_h%ev_ovlp_min = 1.0e-4_r8

   ! Upper bound of the eigenvalues of S
   elsi_h%ev_ovlp_max = 2.0_r8

   ! A patameter used to estimate eigenspectrum
   elsi_h%beta = -1.0e3_r8

end subroutine

!>
!! This routine prints CheSS settings.
!!          
subroutine elsi_print_chess_options(elsi_h)

   implicit none

   type(elsi_handle), intent(in) :: elsi_h !< Handle

   character*200 :: info_str

   character*40, parameter :: caller = "elsi_print_chess_options"

   write(info_str,"(A)") "  CheSS settings (in the same unit of Hamiltonian):"
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Initial guess of error function decay length ',E10.1)") &
      elsi_h%erf_decay
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Lower bound of decay length ',E10.1)") &
      elsi_h%erf_decay_min
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Upper bound of decay length ',E10.1)") &
      elsi_h%erf_decay_max
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Lower bound of H eigenvalue ',E10.1)") &
      elsi_h%ev_ham_min
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Upper bound of H eigenvalue ',E10.1)") &
      elsi_h%ev_ham_max
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Lower bound of S eigenvalue ',E10.1)") &
      elsi_h%ev_ovlp_min
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Upper bound of S eigenvalue ',E10.1)") &
      elsi_h%ev_ovlp_max
   call elsi_statement_print(info_str,elsi_h)

end subroutine

end module ELSI_CHESS
