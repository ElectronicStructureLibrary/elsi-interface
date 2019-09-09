! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Perform wavefunction and density matrix extrapolation between geometry steps.
!!
module ELSI_GEO

   use ELSI_CONSTANT, only: ELPA_SOLVER,EIGENEXA_SOLVER,PEXSI_CSC,SIESTA_CSC,&
       GENERIC_COO
   use ELSI_DATATYPE, only: elsi_handle
   use ELSI_ELPA, only: elsi_update_dm_elpa
   use ELSI_MALLOC, only: elsi_allocate
   use ELSI_MPI, only: elsi_stop
   use ELSI_NTPOLY, only: elsi_update_dm_ntpoly,CopyMatrix
   use ELSI_PRECISION, only: r8
   use ELSI_REDIST, only: elsi_blacs_to_generic_dm,elsi_blacs_to_ntpoly_hs,&
       elsi_blacs_to_siesta_dm,elsi_blacs_to_sips_dm,elsi_generic_to_blacs_hs,&
       elsi_generic_to_ntpoly_hs,elsi_ntpoly_to_blacs_dm,&
       elsi_ntpoly_to_generic_dm,elsi_ntpoly_to_siesta_dm,&
       elsi_ntpoly_to_sips_dm,elsi_siesta_to_blacs_hs,elsi_siesta_to_ntpoly_hs,&
       elsi_sips_to_blacs_hs,elsi_sips_to_ntpoly_hs
   use ELSI_UTIL, only: elsi_check_init,elsi_gram_schmidt

   implicit none

   private

   public :: elsi_orthonormalize_ev_real
   public :: elsi_orthonormalize_ev_complex
   public :: elsi_orthonormalize_ev_real_sparse
   public :: elsi_orthonormalize_ev_complex_sparse
   public :: elsi_extrapolate_dm_real
   public :: elsi_extrapolate_dm_complex
   public :: elsi_extrapolate_dm_real_sparse
   public :: elsi_extrapolate_dm_complex_sparse

contains

!>
!! Orthonormalize eigenvectors with respect to an overlap matrix.
!!
subroutine elsi_orthonormalize_ev_real(eh,ovlp,evec)

   implicit none

   type(elsi_handle), intent(in) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   real(kind=r8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   character(len=*), parameter :: caller = "elsi_orthonormalize_ev_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_gram_schmidt(eh%ph,eh%bh,ovlp,evec)

end subroutine

!>
!! Orthonormalize eigenvectors with respect to an overlap matrix.
!!
subroutine elsi_orthonormalize_ev_complex(eh,ovlp,evec)

   implicit none

   type(elsi_handle), intent(in) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   complex(kind=r8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   character(len=*), parameter :: caller = "elsi_orthonormalize_ev_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_gram_schmidt(eh%ph,eh%bh,ovlp,evec)

end subroutine

!>
!! Orthonormalize eigenvectors with respect to an overlap matrix.
!!
subroutine elsi_orthonormalize_ev_real_sparse(eh,ovlp,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp(eh%bh%nnz_l_sp) !< New overlap
   real(kind=r8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   real(kind=r8) :: dummy1(1)
   real(kind=r8) :: dummy2(1,1)
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_orthonormalize_ev_real_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(.not. allocated(eh%ovlp_real_den)) then
      call elsi_allocate(eh%bh,eh%ovlp_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
           "ovlp_real_den",caller)
   end if

   eh%ph%unit_ovlp = .true.

   select case(eh%ph%matrix_format)
   case(PEXSI_CSC)
      call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp1,&
           eh%col_ptr_sp1,eh%ovlp_real_den,dummy2)
   case(SIESTA_CSC)
      call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp2,&
           eh%col_ptr_sp2,eh%ovlp_real_den,dummy2)
   case(GENERIC_COO)
      if(.not. allocated(eh%map_den)) then
         call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "map_den",caller)
      end if

      call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp3,&
           eh%col_ind_sp3,eh%ovlp_real_den,dummy2,eh%map_den)
   case default
      write(msg,"(A)") "Unsupported matrix format"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%unit_ovlp = .false.

   call elsi_gram_schmidt(eh%ph,eh%bh,eh%ovlp_real_den,evec)

end subroutine

!>
!! Orthonormalize eigenvectors with respect to an overlap matrix.
!!
subroutine elsi_orthonormalize_ev_complex_sparse(eh,ovlp,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp(eh%bh%nnz_l_sp) !< New overlap
   complex(kind=r8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   complex(kind=r8) :: dummy1(1)
   complex(kind=r8) :: dummy2(1,1)
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_orthonormalize_ev_complex_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   if(.not. allocated(eh%ovlp_cmplx_den)) then
      call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
           "ovlp_cmplx_den",caller)
   end if

   eh%ph%unit_ovlp = .true.

   select case(eh%ph%matrix_format)
   case(PEXSI_CSC)
      call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp1,&
           eh%col_ptr_sp1,eh%ovlp_cmplx_den,dummy2)
   case(SIESTA_CSC)
      call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp2,&
           eh%col_ptr_sp2,eh%ovlp_cmplx_den,dummy2)
   case(GENERIC_COO)
      if(.not. allocated(eh%map_den)) then
         call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "map_den",caller)
      end if

      call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp3,&
           eh%col_ind_sp3,eh%ovlp_cmplx_den,dummy2,eh%map_den)
   case default
      write(msg,"(A)") "Unsupported matrix format"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%unit_ovlp = .false.

   call elsi_gram_schmidt(eh%ph,eh%bh,eh%ovlp_cmplx_den,evec)

end subroutine

!>
!! Extrapolate density matrix for a new overlap.
!!
subroutine elsi_extrapolate_dm_real(eh,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< New overlap
   real(kind=r8), intent(inout) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix

   real(kind=r8) :: dummy(1,1)

   character(len=*), parameter :: caller = "elsi_extrapolate_dm_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   select case(eh%ph%solver)
   case(ELPA_SOLVER,EIGENEXA_SOLVER)
      call elsi_update_dm_elpa(eh%ph,eh%bh,eh%ovlp_real_copy,ovlp,dm)

      eh%ovlp_real_copy = ovlp
   case default
      eh%ph%unit_ovlp = .true.

      call elsi_blacs_to_ntpoly_hs(eh%ph,eh%bh,dm,dummy,eh%ph%nt_ham,&
           eh%ph%nt_dm)

      eh%ph%first_blacs_to_ntpoly = .true.

      call elsi_blacs_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy,eh%ph%nt_ovlp,&
           eh%ph%nt_dm)

      eh%ph%unit_ovlp = .false.

      call elsi_update_dm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ovlp_copy,eh%ph%nt_ovlp,&
           eh%ph%nt_ham,eh%ph%nt_dm)
      call elsi_ntpoly_to_blacs_dm(eh%ph,eh%bh,eh%ph%nt_dm,dm)
      call CopyMatrix(eh%ph%nt_ovlp,eh%ph%nt_ovlp_copy)
   end select

end subroutine

!>
!! Extrapolate density matrix for a new overlap.
!!
subroutine elsi_extrapolate_dm_complex(eh,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< New overlap
   complex(kind=r8), intent(inout) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix

   complex(kind=r8) :: dummy(1,1)

   character(len=*), parameter :: caller = "elsi_extrapolate_dm_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_update_dm_elpa(eh%ph,eh%bh,eh%ovlp_cmplx_copy,ovlp,dm)

      eh%ovlp_cmplx_copy = ovlp
   case default
      eh%ph%unit_ovlp = .true.

      call elsi_blacs_to_ntpoly_hs(eh%ph,eh%bh,dm,dummy,eh%ph%nt_ham,&
           eh%ph%nt_dm)

      eh%ph%first_blacs_to_ntpoly = .true.

      call elsi_blacs_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy,eh%ph%nt_ovlp,&
           eh%ph%nt_dm)

      eh%ph%unit_ovlp = .false.

      call elsi_update_dm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ovlp_copy,eh%ph%nt_ovlp,&
           eh%ph%nt_ham,eh%ph%nt_dm)
      call elsi_ntpoly_to_blacs_dm(eh%ph,eh%bh,eh%ph%nt_dm,dm)
      call CopyMatrix(eh%ph%nt_ovlp,eh%ph%nt_ovlp_copy)
   end select

end subroutine

!>
!! Extrapolate density matrix for a new overlap.
!!
subroutine elsi_extrapolate_dm_real_sparse(eh,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(in) :: ovlp(eh%bh%nnz_l_sp) !< New overlap
   real(kind=r8), intent(out) :: dm(eh%bh%nnz_l_sp) !< New density matrix

   real(kind=r8) :: dummy1(1)
   real(kind=r8) :: dummy2(1,1)
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_extrapolate_dm_real_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   select case(eh%ph%solver)
   case(ELPA_SOLVER,EIGENEXA_SOLVER)
      eh%ph%unit_ovlp = .true.

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%ovlp_real_den,dummy2)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ovlp_real_den,dummy2)
      case(GENERIC_COO)
         if(.not. allocated(eh%map_den)) then
            call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "map_den",caller)
         end if

         call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ovlp_real_den,dummy2,eh%map_den)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      eh%ph%unit_ovlp = .false.

      call elsi_update_dm_elpa(eh%ph,eh%bh,eh%ovlp_real_copy,eh%ovlp_real_den,&
           eh%ham_real_den)

      eh%ovlp_real_copy = eh%ovlp_real_den

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%ham_real_den,dm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%ham_real_den,dm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_blacs_to_generic_dm(eh%ph,eh%bh,eh%ham_real_den,eh%map_den,&
              dm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select
   case default
      eh%ph%unit_ovlp = .true.

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%ph%nt_ovlp,eh%ph%nt_dm)
      case(SIESTA_CSC)
         call elsi_siesta_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ph%nt_ovlp,eh%ph%nt_dm)
      case(GENERIC_COO)
         call elsi_generic_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ph%nt_ovlp,eh%ph%nt_dm,eh%ph%nt_map)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      eh%ph%unit_ovlp = .false.

      call elsi_update_dm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ovlp_copy,eh%ph%nt_ovlp,&
           eh%ph%nt_ham,eh%ph%nt_dm)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_ntpoly_to_sips_dm(eh%ph,eh%bh,eh%ph%nt_dm,dm,eh%row_ind_sp1,&
              eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_ntpoly_to_siesta_dm(eh%ph,eh%bh,eh%ph%nt_dm,dm,&
              eh%row_ind_sp2,eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_ntpoly_to_generic_dm(eh%ph,eh%bh,eh%ph%nt_dm,eh%ph%nt_map,&
              dm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call CopyMatrix(eh%ph%nt_ovlp,eh%ph%nt_ovlp_copy)
   end select

end subroutine

!>
!! Extrapolate density matrix for a new overlap.
!!
subroutine elsi_extrapolate_dm_complex_sparse(eh,ovlp,dm)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(in) :: ovlp(eh%bh%nnz_l_sp) !< New overlap
   complex(kind=r8), intent(out) :: dm(eh%bh%nnz_l_sp) !< New density matrix

   complex(kind=r8) :: dummy1(1)
   complex(kind=r8) :: dummy2(1,1)
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_extrapolate_dm_complex_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      eh%ph%unit_ovlp = .true.

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%ovlp_cmplx_den,dummy2)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ovlp_cmplx_den,dummy2)
      case(GENERIC_COO)
         if(.not. allocated(eh%map_den)) then
            call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "map_den",caller)
         end if

         call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ovlp_cmplx_den,dummy2,eh%map_den)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      eh%ph%unit_ovlp = .false.

      call elsi_update_dm_elpa(eh%ph,eh%bh,eh%ovlp_cmplx_copy,&
           eh%ovlp_cmplx_den,eh%ham_cmplx_den)

      eh%ovlp_cmplx_copy = eh%ovlp_cmplx_den

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%ham_cmplx_den,dm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%ham_cmplx_den,dm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_blacs_to_generic_dm(eh%ph,eh%bh,eh%ham_cmplx_den,eh%map_den,&
              dm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select
   case default
      eh%ph%unit_ovlp = .true.

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%ph%nt_ovlp,eh%ph%nt_dm)
      case(SIESTA_CSC)
         call elsi_siesta_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ph%nt_ovlp,eh%ph%nt_dm)
      case(GENERIC_COO)
         call elsi_generic_to_ntpoly_hs(eh%ph,eh%bh,ovlp,dummy1,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ph%nt_ovlp,eh%ph%nt_dm,eh%ph%nt_map)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      eh%ph%unit_ovlp = .false.

      call elsi_update_dm_ntpoly(eh%ph,eh%bh,eh%ph%nt_ovlp_copy,eh%ph%nt_ovlp,&
           eh%ph%nt_ham,eh%ph%nt_dm)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_ntpoly_to_sips_dm(eh%ph,eh%bh,eh%ph%nt_dm,dm,eh%row_ind_sp1,&
              eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_ntpoly_to_siesta_dm(eh%ph,eh%bh,eh%ph%nt_dm,dm,&
              eh%row_ind_sp2,eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_ntpoly_to_generic_dm(eh%ph,eh%bh,eh%ph%nt_dm,eh%ph%nt_map,&
              dm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call CopyMatrix(eh%ph%nt_ovlp,eh%ph%nt_ovlp_copy)
   end select

end subroutine

end module ELSI_GEO
