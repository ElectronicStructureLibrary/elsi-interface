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

   use ISO_C_BINDING
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4
   use ELSI_TIMERS
   use ELSI_UTILS
!TODO: use CHESS

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

   character*40, parameter :: caller = "elsi_init_chess"

   ! Initialize Futile
   call f_lib_initialize()

   ! Initialize sparsematrix error handling and timing
   call sparsematrix_init_errors()
   call sparsematrix_initialize_timing_categories()

   ! Initialize sparse matrices
   call sparse_matrix_init_from_data_ccs()

   ! Initialize task groups
   call init_matrix_taskgroups_wrapper(elsi_h%myid,elsi_h%n_procs,&
           elsi_h%mpi_comm,.false.,2,elsi_h%sparse_mat)

   ! Initialize FOE objects
   call init_foe()

   ! Allocate CheSS matrices
   call matrices_init()

end subroutine

!>
!! This routine interfaces to CheSS.
!!
subroutine elsi_solve_evp_chess(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   integer(kind=i4) :: mpierr

   character*40, parameter :: caller = "elsi_solve_evp_chess"

   call elsi_start_density_matrix_time(elsi_h)

   call matrix_fermi_operator_expansion()

   call MPI_Barrier(elsi_h%mpi_comm,mpierr)
   call elsi_stop_density_matrix_time(elsi_h)

end subroutine

!> 
!! This routine computes the energy-weighted density matrix.
!! 
subroutine elsi_compute_edm_chess(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   character*40, parameter :: caller = "elsi_compute_edm_chess"

end subroutine

!> 
!! This routine sets default CheSS parameters.
!! 
subroutine elsi_set_chess_default(elsi_h)
   
   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

   character*40, parameter :: caller = "elsi_set_chess_default"

   ! Initial guess of the decay length of the error function
   elsi_h%erf_decay = 0.1_r8

   ! Lower bound of the decay length
   elsi_h%erf_decay_min = 0.01_r8

   ! Upper bound of the decay length
   elsi_h%erf_decay_max = 0.1_r8

   ! Lower bound of the eigenvalues of H
   elsi_h%ev_ham_min = -20.0_r8

   ! Upper bound of the eigenvalues of H
   elsi_h%ev_ham_max = 10.0_r8

   ! Lower bound of the eigenvalues of S
   elsi_h%ev_ovlp_min = 1.0e-6_r8

   ! Upper bound of the eigenvalues of S
   elsi_h%ev_ovlp_max = 5.0_r8

   ! ???
   elsi_h%betax = -1.0e3_r8

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

   write(info_str,"(1X,' | Initial guess of error function decay length ',F10.4)") &
      elsi_h%erf_decay
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Lower bound of decay length ',F10.4)") &
      elsi_h%erf_decay_min
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Upper bound of decay length ',F10.4)") &
      elsi_h%erf_decay_max
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Lower bound of H eigenvalue ',F10.4)") &
      elsi_h%ev_ham_min
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Upper bound of H eigenvalue ',F10.4)") &
      elsi_h%ev_ham_max
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Lower bound of S eigenvalue ',F10.4)") &
      elsi_h%ev_ovlp_min
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Upper bound of S eigenvalue ',F10.4)") &
      elsi_h%ev_ovlp_max
   call elsi_statement_print(info_str,elsi_h)

   write(info_str,"(1X,' | Betax ',F10.4)") &
      elsi_h%betax
   call elsi_statement_print(info_str,elsi_h)

end subroutine

end module ELSI_CHESS
