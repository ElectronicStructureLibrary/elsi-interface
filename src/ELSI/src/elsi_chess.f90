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

   public :: elsi_set_chess_default
   public :: elsi_init_chess
   public :: elsi_solve_evp_chess
   public :: elsi_compute_edm_chess

contains

!>
!! This routine initializes CheSS.
!!
subroutine elsi_init_chess(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   real(kind=r8) :: n_electron(1)

   character*40, parameter :: caller = "elsi_init_chess"

   if(e_h%n_elsi_calls == 1) then
      ! Initialize f_lib
      call f_lib_initialize()

      ! Initialize sparsematrix error handling and timing
      call sparsematrix_init_errors()
      call sparsematrix_initialize_timing_categories()

      ! Initialize sparse matrices
      call sparse_matrix_init_from_data_ccs(e_h%myid,e_h%n_procs,e_h%mpi_comm,&
              e_h%n_basis,e_h%nnz_g,e_h%row_ind_chess,e_h%col_ptr_chess,&
              e_h%sparse_mat(1),init_matmul=.false.)

      ! TODO: Sparsity buffer
      call sparse_matrix_init_from_data_ccs(e_h%myid,e_h%n_procs,e_h%mpi_comm,&
              e_h%n_basis,e_h%nnz_g,e_h%row_ind_chess,e_h%col_ptr_chess,&
              e_h%sparse_mat(2),init_matmul=.true.,nvctr_mult=e_h%nnz_g,&
              row_ind_mult=e_h%row_ind_chess,&
              col_ptr_mult=e_h%col_ptr_chess(1:e_h%n_basis))

      ! Initialize task groups
      call init_matrix_taskgroups_wrapper(e_h%myid,e_h%n_procs,e_h%mpi_comm,&
              .false.,2,e_h%sparse_mat)

      n_electron = e_h%n_electrons

      ! Initialize FOE objects
      call init_foe(e_h%myid,e_h%n_procs,1,n_electron,e_h%foe_obj,&
              fscale=e_h%erf_decay,fscale_lowerbound=e_h%erf_decay_min,&
              fscale_upperbound=e_h%erf_decay_max,evlow=e_h%ev_ham_min,&
              evhigh=e_h%ev_ham_max,betax=e_h%beta)

      call init_foe(e_h%myid,e_h%n_procs,1,n_electron,e_h%ice_obj,&
              evlow=e_h%ev_ovlp_min,evhigh=e_h%ev_ovlp_max,betax=e_h%beta)

      ! Allocate CheSS matrices
      call matrices_init(e_h%sparse_mat(1),e_h%ham_chess,&
              matsize=SPARSE_TASKGROUP)
      call matrices_init(e_h%sparse_mat(1),e_h%ovlp_chess,&
              matsize=SPARSE_TASKGROUP)
      call matrices_init(e_h%sparse_mat(2),e_h%dm_chess,&
              matsize=SPARSE_TASKGROUP)
      call matrices_init(e_h%sparse_mat(2),e_h%edm_chess,&
              matsize=SPARSE_TASKGROUP)
      call matrices_init(e_h%sparse_mat(2),e_h%ovlp_inv_sqrt(1),&
              matsize=SPARSE_TASKGROUP)

      e_h%ovlp_chess%matrix_compr = e_h%ovlp_real_chess

      e_h%chess_started = .true.
   endif

   e_h%ham_chess%matrix_compr = e_h%ham_real_chess

end subroutine

!>
!! This routine interfaces to CheSS.
!!
subroutine elsi_solve_evp_chess(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   logical          :: calc_ovlp_inv_sqrt
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   integer(kind=i4) :: mpierr
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_solve_evp_chess"

   call elsi_get_time(e_h,t0)

   if(e_h%n_elsi_calls == 1) then
      calc_ovlp_inv_sqrt = .true.
   else
      calc_ovlp_inv_sqrt = .false.
   endif

   call elsi_statement_print("  Starting CheSS density matrix solver",e_h)

   call matrix_fermi_operator_expansion(e_h%myid,e_h%n_procs,e_h%mpi_comm,&
           e_h%foe_obj,e_h%ice_obj,e_h%sparse_mat(1),e_h%sparse_mat(1),&
           e_h%sparse_mat(2),e_h%ovlp_chess,e_h%ham_chess,e_h%ovlp_inv_sqrt,&
           e_h%dm_chess,e_h%energy_hdm,&
           calculate_minusonehalf=calc_ovlp_inv_sqrt,foe_verbosity=0,&
           symmetrize_kernel=.true.,calculate_energy_density_kernel=.false.,&
           energy_kernel=e_h%edm_chess)

   e_h%mu = foe_data_get_real(e_h%foe_obj,"ef",1)

   call MPI_Barrier(e_h%mpi_comm,mpierr)

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished density matrix calculation')")
   call elsi_statement_print(info_str,e_h)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_statement_print(info_str,e_h)

end subroutine

!>
!! This routine computes the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_chess(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character*40, parameter :: caller = "elsi_compute_edm_chess"

   ! TODO: Compute edm

end subroutine

!>
!! This routine sets default CheSS parameters.
!!
subroutine elsi_set_chess_default(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character*40, parameter :: caller = "elsi_set_chess_default"

   ! Initial guess of the decay length of the error function
   e_h%erf_decay = 5.0e-2_r8

   ! Lower bound of the decay length
   e_h%erf_decay_min = 5.0e-3_r8

   ! Upper bound of the decay length
   e_h%erf_decay_max = 5.0e-2_r8

   ! Lower bound of the eigenvalues of H
   e_h%ev_ham_min = -2.0_r8

   ! Upper bound of the eigenvalues of H
   e_h%ev_ham_max = 2.0_r8

   ! Lower bound of the eigenvalues of S
   e_h%ev_ovlp_min = 1.0e-4_r8

   ! Upper bound of the eigenvalues of S
   e_h%ev_ovlp_max = 2.0_r8

   ! A patameter used to estimate eigenspectrum
   e_h%beta = -1.0e3_r8

end subroutine

end module ELSI_CHESS
