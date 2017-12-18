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
!! This module provides interfaces to libOMM.
!!
module ELSI_OMM

   use ELSI_CONSTANTS, only: REAL_VALUES,COMPLEX_VALUES,BLACS_DENSE
   use ELSI_DATATYPE
   use ELSI_PRECISION, only: r8,i4
   use ELSI_TIMINGS,   only: elsi_get_time
   use ELSI_UTILS
   use ELPA1,          only: elpa_cholesky_real_double,&
                             elpa_cholesky_complex_double,&
                             elpa_invert_trm_real_double,&
                             elpa_invert_trm_complex_double
   use MATRIXSWITCH,   only: m_add

   implicit none

   private

   public :: elsi_set_omm_default
   public :: elsi_solve_evp_omm
   public :: elsi_compute_edm_omm

contains

!>
!! This routine interfaces to libOMM.
!!
subroutine elsi_solve_evp_omm(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   logical          :: success
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   integer(kind=i4) :: mpierr
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_solve_evp_omm"

   ! Compute sparsity
   if(e_h%n_elsi_calls == 1 .and. e_h%matrix_format == BLACS_DENSE) then
      select case(e_h%data_type)
      case(COMPLEX_VALUES)
         call elsi_get_local_nnz_complex(e_h,e_h%ham_omm%zval,e_h%n_lrow,&
                 e_h%n_lcol,e_h%nnz_l)

         call MPI_Allreduce(e_h%nnz_l,e_h%nnz_g,1,mpi_integer4,mpi_sum,&
                 e_h%mpi_comm,mpierr)
      case(REAL_VALUES)
         call elsi_get_local_nnz_real(e_h,e_h%ham_omm%dval,e_h%n_lrow,&
                 e_h%n_lcol,e_h%nnz_l)

         call MPI_Allreduce(e_h%nnz_l,e_h%nnz_g,1,mpi_integer4,mpi_sum,&
                 e_h%mpi_comm,mpierr)
      end select
   endif

   if(.not. e_h%ovlp_is_unit) then
      if(e_h%omm_flavor == 2) then
         if(e_h%n_elsi_calls == 1) then
            call elsi_get_time(e_h,t0)

            ! Cholesky factorization
            select case(e_h%data_type)
            case(COMPLEX_VALUES)
               ! Compute S = (U^T)U, U -> S
               success = elpa_cholesky_complex_double(e_h%n_basis,&
                            e_h%ovlp_omm%zval,e_h%n_lrow,e_h%blk_row,&
                            e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                            .false.)

               success = elpa_invert_trm_complex_double(e_h%n_basis,&
                            e_h%ovlp_omm%zval,e_h%n_lrow,e_h%blk_row,&
                            e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                            .false.)
            case(REAL_VALUES)
               ! Compute S = (U^T)U, U -> S
               success = elpa_cholesky_real_double(e_h%n_basis,&
                            e_h%ovlp_omm%dval,e_h%n_lrow,e_h%blk_row,&
                            e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                            .false.)

               success = elpa_invert_trm_real_double(e_h%n_basis,&
                            e_h%ovlp_omm%dval,e_h%n_lrow,e_h%blk_row,&
                            e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                            .false.)
            end select

            call elsi_get_time(e_h,t1)

            write(info_str,"('  Finished Cholesky decomposition')")
            call elsi_say(e_h,info_str)
            write(info_str,"('  | Time :',F10.3,' s')") t1-t0
            call elsi_say(e_h,info_str)
         endif

         if(e_h%n_elsi_calls > e_h%omm_n_elpa+1) then
            ! Invert one more time
            select case(e_h%data_type)
            case(COMPLEX_VALUES)
               success = elpa_invert_trm_complex_double(e_h%n_basis,&
                            e_h%ovlp_omm%zval,e_h%n_lrow,e_h%blk_row,&
                            e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                            .false.)
            case(REAL_VALUES)
               success = elpa_invert_trm_real_double(e_h%n_basis,&
                            e_h%ovlp_omm%dval,e_h%n_lrow,e_h%blk_row,&
                            e_h%n_lcol,e_h%mpi_comm_row,e_h%mpi_comm_col,&
                            .false.)
            end select
         endif
      endif ! omm_flavor == 2
   endif ! ovlp_is_unit

   if(e_h%n_elsi_calls == 1) then
      e_h%coeff_ready = .false.
   else
      e_h%coeff_ready = .true.
   endif

   if(e_h%n_elsi_calls == e_h%omm_n_elpa+1) then
      e_h%new_overlap = .true.
   else
      e_h%new_overlap = .false.
   endif

   call elsi_get_time(e_h,t0)

   call elsi_say(e_h,"  Starting OMM density matrix solver")

   select case(e_h%data_type)
   case(COMPLEX_VALUES)
      call omm(e_h%n_basis,e_h%n_states_omm,e_h%ham_omm,e_h%ovlp_omm,&
              e_h%new_overlap,e_h%energy_hdm,e_h%dm_omm,e_h%calc_ed,e_h%eta,&
              e_h%coeff,e_h%coeff_ready,e_h%tdm_omm,e_h%scale_kinetic,&
              e_h%omm_flavor,1,1,e_h%min_tol,e_h%omm_output,e_h%do_dealloc,&
              "pzdbc","lap")

   case(REAL_VALUES)
      call omm(e_h%n_basis,e_h%n_states_omm,e_h%ham_omm,e_h%ovlp_omm,&
              e_h%new_overlap,e_h%energy_hdm,e_h%dm_omm,e_h%calc_ed,e_h%eta,&
              e_h%coeff,e_h%coeff_ready,e_h%tdm_omm,e_h%scale_kinetic,&
              e_h%omm_flavor,1,1,e_h%min_tol,e_h%omm_output,e_h%do_dealloc,&
              "pddbc","lap")
   end select

   call MPI_Barrier(e_h%mpi_comm,mpierr)

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished density matrix calculation')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine computes the energy-weighted density matrix.
!!
subroutine elsi_compute_edm_omm(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character*200 :: info_str

   character*40, parameter :: caller = "elsi_compute_edm_omm"

   call elsi_get_time(e_h,t0)

   e_h%calc_ed = .true.

   select case(e_h%data_type)
   case(COMPLEX_VALUES)
      call omm(e_h%n_basis,e_h%n_states_omm,e_h%ham_omm,e_h%ovlp_omm,&
              e_h%new_overlap,e_h%energy_hdm,e_h%dm_omm,e_h%calc_ed,e_h%eta,&
              e_h%coeff,e_h%coeff_ready,e_h%tdm_omm,e_h%scale_kinetic,&
              e_h%omm_flavor,1,1,e_h%min_tol,e_h%omm_output,e_h%do_dealloc,&
              "pzdbc","lap")
   case(REAL_VALUES)
      call omm(e_h%n_basis,e_h%n_states_omm,e_h%ham_omm,e_h%ovlp_omm,&
              e_h%new_overlap,e_h%energy_hdm,e_h%dm_omm,e_h%calc_ed,e_h%eta,&
              e_h%coeff,e_h%coeff_ready,e_h%tdm_omm,e_h%scale_kinetic,&
              e_h%omm_flavor,1,1,e_h%min_tol,e_h%omm_output,e_h%do_dealloc,&
              "pddbc","lap")
   end select

   e_h%calc_ed = .false.

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished energy density matrix calculation')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine sets default libOMM parameters.
!!
subroutine elsi_set_omm_default(e_h)

   implicit none

   type(elsi_handle), intent(inout) :: e_h !< Handle

   character*40, parameter :: caller = "elsi_set_omm_default"

   ! How many steps of ELPA to run before OMM
   e_h%omm_n_elpa = 6

   ! How do we perform the calculation
   ! 0 = Basic
   ! 2 = Cholesky already performed, U is provided in S
   e_h%omm_flavor = 0

   ! How to scale the kinetic energy matrix
   e_h%scale_kinetic = 5.0_r8

   ! Calculate the energy density matrix
   e_h%calc_ed = .false.

   ! Eigenspectrum shift parameter
   e_h%eta = 0.0_r8

   ! Tolerance for minimization
   e_h%min_tol = 1.0e-12_r8

   ! Output level?
   e_h%omm_output = .false.

   ! Deallocate temporary arrays?
   e_h%do_dealloc = .false.

   ! Use pspBLAS sparse linear algebra?
   e_h%use_psp = .false.

end subroutine

end module ELSI_OMM
