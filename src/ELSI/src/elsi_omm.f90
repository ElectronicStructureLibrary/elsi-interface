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

   use ELSI_CONSTANTS, only: BLACS_DENSE
   use ELSI_DATATYPE
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS
   use ELPA1,          only: elpa_cholesky_real_double,&
                             elpa_cholesky_complex_double,&
                             elpa_invert_trm_real_double,&
                             elpa_invert_trm_complex_double
   use MATRIXSWITCH,   only: m_add,m_register_pdbc

   implicit none

   private

   public :: elsi_set_omm_default
   public :: elsi_solve_evp_omm_real
   public :: elsi_compute_edm_omm_real
   public :: elsi_solve_evp_omm_cmplx
   public :: elsi_compute_edm_omm_cmplx

contains

!>
!! This routine interfaces to libOMM.
!!
subroutine elsi_solve_evp_omm_real(e_h,ham,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: dm(e_h%n_lrow,e_h%n_lcol)

   logical          :: success
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   integer(kind=i4) :: mpierr
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_solve_evp_omm_real"

   call m_register_pdbc(e_h%ham_omm,ham,e_h%sc_desc)
   call m_register_pdbc(e_h%ovlp_omm,ovlp,e_h%sc_desc)
   call m_register_pdbc(e_h%dm_omm,dm,e_h%sc_desc)

   ! Compute sparsity
   if(e_h%n_elsi_calls == 1 .and. e_h%matrix_format == BLACS_DENSE) then
      call elsi_get_local_nnz_real(e_h,ham,e_h%n_lrow,e_h%n_lcol,e_h%nnz_l)

      call MPI_Allreduce(e_h%nnz_l,e_h%nnz_g,1,mpi_integer4,mpi_sum,&
              e_h%mpi_comm,mpierr)
   endif

   if(.not. e_h%ovlp_is_unit) then
      if(e_h%omm_flavor == 2) then
         if(e_h%n_elsi_calls == 1) then
            call elsi_get_time(e_h,t0)

            ! Cholesky factorization
            success = elpa_cholesky_real_double(e_h%n_basis,ovlp,e_h%n_lrow,&
                         e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,&
                         e_h%mpi_comm_col,.false.)

            success = elpa_invert_trm_real_double(e_h%n_basis,ovlp,e_h%n_lrow,&
                         e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,&
                         e_h%mpi_comm_col,.false.)

            call elsi_get_time(e_h,t1)

            write(info_str,"('  Finished Cholesky decomposition')")
            call elsi_say(e_h,info_str)
            write(info_str,"('  | Time :',F10.3,' s')") t1-t0
            call elsi_say(e_h,info_str)
         endif

         if(e_h%n_elsi_calls > e_h%omm_n_elpa+1) then
            success = elpa_invert_trm_real_double(e_h%n_basis,ovlp,e_h%n_lrow,&
                         e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,&
                         e_h%mpi_comm_col,.false.)
         endif
      endif ! omm_flavor == 2
   endif ! ovlp_is_unit

   if(e_h%n_elsi_calls == 1) then
      e_h%coeff_ready = .false.
   else
      e_h%coeff_ready = .true.
   endif

   if(e_h%n_elsi_calls == e_h%omm_n_elpa+1) then
      e_h%new_ovlp = .true.
   else
      e_h%new_ovlp = .false.
   endif

   call elsi_get_time(e_h,t0)

   call elsi_say(e_h,"  Starting OMM density matrix solver")

   call omm(e_h%n_basis,e_h%n_states_omm,e_h%ham_omm,e_h%ovlp_omm,e_h%new_ovlp,&
           e_h%energy_hdm,e_h%dm_omm,e_h%calc_ed,e_h%eta,e_h%coeff,&
           e_h%coeff_ready,e_h%tdm_omm,e_h%scale_kinetic,e_h%omm_flavor,1,1,&
           e_h%min_tol,e_h%omm_output,e_h%do_dealloc,"pddbc","lap")

   dm = e_h%spin_degen*dm

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
subroutine elsi_compute_edm_omm_real(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   real(kind=r8),     intent(inout) :: edm(e_h%n_lrow,e_h%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character*200 :: info_str

   character*40, parameter :: caller = "elsi_compute_edm_omm_real"

   call elsi_get_time(e_h,t0)

   call m_register_pdbc(e_h%dm_omm,edm,e_h%sc_desc)

   e_h%calc_ed = .true.

   call omm(e_h%n_basis,e_h%n_states_omm,e_h%ham_omm,e_h%ovlp_omm,e_h%new_ovlp,&
           e_h%energy_hdm,e_h%dm_omm,e_h%calc_ed,e_h%eta,e_h%coeff,&
           e_h%coeff_ready,e_h%tdm_omm,e_h%scale_kinetic,e_h%omm_flavor,1,1,&
           e_h%min_tol,e_h%omm_output,e_h%do_dealloc,"pddbc","lap")

   e_h%calc_ed = .false.

   edm = e_h%spin_degen*edm

   call elsi_get_time(e_h,t1)

   write(info_str,"('  Finished energy density matrix calculation')")
   call elsi_say(e_h,info_str)
   write(info_str,"('  | Time :',F10.3,' s')") t1-t0
   call elsi_say(e_h,info_str)

end subroutine

!>
!! This routine interfaces to libOMM.
!!
subroutine elsi_solve_evp_omm_cmplx(e_h,ham,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(inout) :: ham(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: ovlp(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),  intent(inout) :: dm(e_h%n_lrow,e_h%n_lcol)

   logical          :: success
   real(kind=r8)    :: t0
   real(kind=r8)    :: t1
   integer(kind=i4) :: mpierr
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_solve_evp_omm_cmplx"

   call m_register_pdbc(e_h%ham_omm,ham,e_h%sc_desc)
   call m_register_pdbc(e_h%ovlp_omm,ovlp,e_h%sc_desc)
   call m_register_pdbc(e_h%dm_omm,dm,e_h%sc_desc)

   ! Compute sparsity
   if(e_h%n_elsi_calls == 1 .and. e_h%matrix_format == BLACS_DENSE) then
      call elsi_get_local_nnz_cmplx(e_h,ham,e_h%n_lrow,e_h%n_lcol,e_h%nnz_l)

      call MPI_Allreduce(e_h%nnz_l,e_h%nnz_g,1,mpi_integer4,mpi_sum,&
              e_h%mpi_comm,mpierr)
   endif

   if(.not. e_h%ovlp_is_unit) then
      if(e_h%omm_flavor == 2) then
         if(e_h%n_elsi_calls == 1) then
            call elsi_get_time(e_h,t0)

            ! Cholesky factorization
            success = elpa_cholesky_complex_double(e_h%n_basis,ovlp,e_h%n_lrow,&
                         e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,&
                         e_h%mpi_comm_col,.false.)

            success = elpa_invert_trm_complex_double(e_h%n_basis,ovlp,&
                         e_h%n_lrow,e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,&
                         e_h%mpi_comm_col,.false.)

            call elsi_get_time(e_h,t1)

            write(info_str,"('  Finished Cholesky decomposition')")
            call elsi_say(e_h,info_str)
            write(info_str,"('  | Time :',F10.3,' s')") t1-t0
            call elsi_say(e_h,info_str)
         endif

         if(e_h%n_elsi_calls > e_h%omm_n_elpa+1) then
            success = elpa_invert_trm_complex_double(e_h%n_basis,ovlp,&
                         e_h%n_lrow,e_h%blk_row,e_h%n_lcol,e_h%mpi_comm_row,&
                         e_h%mpi_comm_col,.false.)
         endif
      endif ! omm_flavor == 2
   endif ! ovlp_is_unit

   if(e_h%n_elsi_calls == 1) then
      e_h%coeff_ready = .false.
   else
      e_h%coeff_ready = .true.
   endif

   if(e_h%n_elsi_calls == e_h%omm_n_elpa+1) then
      e_h%new_ovlp = .true.
   else
      e_h%new_ovlp = .false.
   endif

   call elsi_get_time(e_h,t0)

   call elsi_say(e_h,"  Starting OMM density matrix solver")

   call omm(e_h%n_basis,e_h%n_states_omm,e_h%ham_omm,e_h%ovlp_omm,e_h%new_ovlp,&
           e_h%energy_hdm,e_h%dm_omm,e_h%calc_ed,e_h%eta,e_h%coeff,&
           e_h%coeff_ready,e_h%tdm_omm,e_h%scale_kinetic,e_h%omm_flavor,1,1,&
           e_h%min_tol,e_h%omm_output,e_h%do_dealloc,"pzdbc","lap")

   dm = e_h%spin_degen*dm

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
subroutine elsi_compute_edm_omm_cmplx(e_h,edm)

   implicit none

   type(elsi_handle), intent(inout) :: e_h
   complex(kind=r8),  intent(inout) :: edm(e_h%n_lrow,e_h%n_lcol)

   real(kind=r8) :: t0
   real(kind=r8) :: t1
   character*200 :: info_str

   character*40, parameter :: caller = "elsi_compute_edm_omm_cmplx"

   call elsi_get_time(e_h,t0)

   call m_register_pdbc(e_h%dm_omm,edm,e_h%sc_desc)

   e_h%calc_ed = .true.

   call omm(e_h%n_basis,e_h%n_states_omm,e_h%ham_omm,e_h%ovlp_omm,e_h%new_ovlp,&
           e_h%energy_hdm,e_h%dm_omm,e_h%calc_ed,e_h%eta,e_h%coeff,&
           e_h%coeff_ready,e_h%tdm_omm,e_h%scale_kinetic,e_h%omm_flavor,1,1,&
           e_h%min_tol,e_h%omm_output,e_h%do_dealloc,"pzdbc","lap")

   e_h%calc_ed = .false.

   edm = e_h%spin_degen*edm

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

   type(elsi_handle), intent(inout) :: e_h

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
