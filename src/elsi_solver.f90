! Copyright (c) 2015-2021, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Provide interface routines to eigensolvers and density matrix solvers.
!!
module ELSI_SOLVER

   use ELSI_BSEPACK, only: elsi_solve_bsepack
   use ELSI_CONSTANT, only: ELPA_SOLVER,OMM_SOLVER,PEXSI_SOLVER,SIPS_SOLVER,&
       NTPOLY_SOLVER,EIGENEXA_SOLVER,MAGMA_SOLVER,BSEPACK_SOLVER,MULTI_PROC,&
       SINGLE_PROC,PEXSI_CSC,SIESTA_CSC,GENERIC_COO,GET_DM,GET_EDM
   use ELSI_DATATYPE, only: elsi_handle,elsi_param_t,elsi_basic_t
   use ELSI_DECISION, only: elsi_decide_ev,elsi_decide_dm
   use ELSI_EIGENEXA, only: elsi_init_eigenexa,elsi_solve_eigenexa
   use ELSI_ELPA, only: elsi_init_elpa,elsi_solve_elpa,elsi_do_fc_elpa,&
       elsi_undo_fc_elpa
   use ELSI_LAPACK, only: elsi_solve_lapack
   use ELSI_MAGMA, only: elsi_init_magma,elsi_solve_magma
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI
   use ELSI_NTPOLY, only: elsi_init_ntpoly,elsi_solve_ntpoly
   use ELSI_OCC, only: elsi_mu_and_occ,elsi_entropy,elsi_get_occ_for_dm
   use ELSI_OMM, only: elsi_init_omm,elsi_solve_omm
   use ELSI_OUTPUT, only: elsi_add_log,elsi_get_time,fjson_get_datetime_rfc3339
   use ELSI_PEXSI, only: elsi_init_pexsi,elsi_solve_pexsi
   use ELSI_PRECISION, only: r8,i4
   use ELSI_REDIST, only: elsi_blacs_to_generic_dm,elsi_blacs_to_mask,&
       elsi_blacs_to_ntpoly_hs,elsi_blacs_to_pexsi_hs_dim,&
       elsi_blacs_to_pexsi_hs,elsi_blacs_to_siesta_dm,elsi_blacs_to_sips_dm,&
       elsi_blacs_to_sips_hs_dim,elsi_blacs_to_sips_hs,&
       elsi_generic_to_blacs_hs,elsi_generic_to_ntpoly_hs,&
       elsi_generic_to_pexsi_hs_dim,elsi_generic_to_pexsi_hs,&
       elsi_generic_to_sips_hs_dim,elsi_generic_to_sips_hs,&
       elsi_ntpoly_to_blacs_dm,elsi_ntpoly_to_generic_dm,&
       elsi_ntpoly_to_siesta_dm,elsi_ntpoly_to_sips_dm,elsi_pexsi_to_blacs_dm,&
       elsi_pexsi_to_generic_dm,elsi_pexsi_to_siesta_dm,&
       elsi_siesta_to_blacs_hs,elsi_siesta_to_ntpoly_hs,&
       elsi_siesta_to_pexsi_hs_dim,elsi_siesta_to_pexsi_hs,&
       elsi_siesta_to_sips_hs_dim,elsi_siesta_to_sips_hs,elsi_sips_to_blacs_dm,&
       elsi_sips_to_blacs_ev,elsi_sips_to_blacs_hs,elsi_sips_to_generic_dm,&
       elsi_sips_to_ntpoly_hs,elsi_sips_to_siesta_dm
   use ELSI_SETUP, only: elsi_set_blacs
   use ELSI_SIPS, only: elsi_init_sips,elsi_solve_sips,elsi_build_dm_edm_sips
   use ELSI_UTIL, only: elsi_check,elsi_check_init,elsi_reduce_energy,&
       elsi_build_dm_edm

   implicit none

   private

   public :: elsi_ev_real
   public :: elsi_ev_complex
   public :: elsi_ev_real_sparse
   public :: elsi_ev_complex_sparse
   public :: elsi_dm_real
   public :: elsi_dm_complex
   public :: elsi_dm_real_sparse
   public :: elsi_dm_complex_sparse
   public :: elsi_bse_real
   public :: elsi_bse_complex
   public :: elsi_compute_dm_real
   public :: elsi_compute_dm_complex
   public :: elsi_compute_edm_real
   public :: elsi_compute_edm_complex
   public :: elsi_compute_mu_and_occ
   public :: elsi_compute_entropy

contains

!>
!! Get the band structure energy.
!!
subroutine elsi_get_band_energy(ph,bh,ebs,solver)

   implicit none

   type(elsi_param_t), intent(inout) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: solver
   real(kind=r8), intent(out) :: ebs

   character(len=*), parameter :: caller = "elsi_get_band_energy"

   ebs = ph%ebs*ph%i_wt

   ! Handle different conventions
   if(solver == OMM_SOLVER) then
      ebs = ph%spin_degen*ebs
   else if(solver == NTPOLY_SOLVER) then
      ebs = ph%spin_degen*ebs/2.0_r8
   end if

   call elsi_reduce_energy(ph,bh,ebs)

end subroutine

!>
!! Initialize BLACS, in case that the user selects a sparse format and BLACS is
!! still needed internally.
!!
subroutine elsi_init_blacs(eh)

   implicit none

   type(elsi_handle), intent(inout) :: eh

   integer(kind=i4) :: nprow
   integer(kind=i4) :: npcol
   integer(kind=i4) :: blacs_ctxt
   integer(kind=i4) :: block_size
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_init_blacs"

   if(eh%ph%parallel_mode == MULTI_PROC .and. .not. eh%bh%blacs_ready) then
      ! Set square-like process grid
      do nprow = nint(sqrt(real(eh%bh%n_procs,kind=r8)),kind=i4),2,-1
         if(mod(eh%bh%n_procs,nprow) == 0) then
            exit
         end if
      end do

      npcol = eh%bh%n_procs/nprow

      if(max(nprow,npcol) > eh%ph%n_basis) then
         write(msg,"(A)") "Number of MPI tasks too large"
         call elsi_stop(eh%bh,msg,caller)
      end if

      ! Initialize BLACS
      blacs_ctxt = eh%bh%comm

      call BLACS_Gridinit(blacs_ctxt,"r",nprow,npcol)

      ! Find block size
      block_size = 1

      do while(2*block_size*max(nprow,npcol) <= eh%ph%n_basis)
         block_size = 2*block_size
      end do

      ! Maximum allowed value: 256
      block_size = min(256,block_size)

      ! ELPA works better with a small block_size
      if(eh%ph%solver == ELPA_SOLVER) then
         block_size = min(32,block_size)
      end if

      call elsi_set_blacs(eh,blacs_ctxt,block_size)
   end if

end subroutine

!>
!! Compute the eigenvalues and eigenvectors of a Kohn-Sham problem. Note the
!! intent(inout), everything may be reused in the next call.
!!
subroutine elsi_ev_real(eh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(inout) :: ham(eh%bh%n_lrow,eh%bh%n_lcol) !< Hamiltonian
   real(kind=r8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   real(kind=r8), intent(inout) :: eval(eh%ph%n_basis) !< Eigenvalues
   real(kind=r8), intent(out) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   real(kind=r8) :: t0
   integer(kind=i4) :: solver
   character(len=29) :: dt0
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_ev_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)
   call elsi_decide_ev(eh%ph,eh%bh)

   eh%ph%n_calls = eh%ph%n_calls+1
   solver = eh%ph%solver

   if(eh%ph%solver == SIPS_SOLVER .and. eh%ph%n_calls <= eh%ph%sips_n_elpa) then
      solver = ELPA_SOLVER
   end if

   if(eh%ph%solver /= solver) then
      if(.not. allocated(eh%ovlp_real_copy)) then
         call elsi_allocate(eh%bh,eh%ovlp_real_copy,eh%bh%n_lrow,eh%bh%n_lcol,&
              "ovlp_real_copy",caller)
      end if

      if(eh%ph%n_calls == 1) then
         eh%ovlp_real_copy(:,:) = ovlp
      end if
   end if

   select case(solver)
   case(ELPA_SOLVER)
      if(eh%ph%parallel_mode == SINGLE_PROC) then
         call elsi_solve_lapack(eh%ph,eh%bh,ham,ovlp,eval,evec)
      else
         if(eh%ph%n_basis_c > 0) then
            if(.not. allocated(eh%ham_real_v)) then
               call elsi_allocate(eh%bh,eh%ham_real_v,eh%ph%n_lrow_v,&
                    eh%ph%n_lcol_v,"ham_real_v",caller)
               call elsi_allocate(eh%bh,eh%ovlp_real_v,eh%ph%n_lrow_v,&
                    eh%ph%n_lcol_v,"ovlp_real_v",caller)
               call elsi_allocate(eh%bh,eh%evec_real_v,eh%ph%n_lrow_v,&
                    eh%ph%n_lcol_v,"evec_real_v",caller)
            end if

            call elsi_do_fc_elpa(eh%ph,eh%bh,ham,ovlp,evec,eh%perm_fc,&
                 eh%ham_real_v,eh%ovlp_real_v,eh%evec_real_v)
            call elsi_init_elpa(eh%ph,eh%bh)
            call elsi_solve_elpa(eh%ph,eh%bh,eh%ham_real_v,eh%ovlp_real_v,&
                 eval(eh%ph%n_basis_c+1:eh%ph%n_basis_c+eh%ph%n_basis_v),&
                 eh%evec_real_v)
            call elsi_undo_fc_elpa(eh%ph,eh%bh,ham,ovlp,evec,eh%perm_fc,&
                 eval(1:eh%ph%n_basis_c),eh%evec_real_v)
         else
            call elsi_init_elpa(eh%ph,eh%bh)
            call elsi_solve_elpa(eh%ph,eh%bh,ham,ovlp,eval,evec)
         end if
      end if
   case(EIGENEXA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)
      call elsi_init_eigenexa(eh%ph,eh%bh)
      call elsi_solve_eigenexa(eh%ph,eh%bh,ham,ovlp,eval,evec)
   case(SIPS_SOLVER)
      call elsi_init_sips(eh%ph,eh%bh)

      if(eh%ph%n_calls > 1 .and. eh%ph%n_calls == eh%ph%sips_n_elpa+1) then
         ! Restore overlap
         ovlp(:,:) = eh%ovlp_real_copy
      end if

      if(.not. allocated(eh%row_ind_sp1)) then
         call elsi_allocate(eh%bh,eh%mask,eh%bh%n_lrow,eh%bh%n_lcol,"mask",&
              caller)

         call elsi_blacs_to_mask(eh%ph,eh%bh,ham,ovlp,eh%mask)
         call elsi_blacs_to_sips_hs_dim(eh%ph,eh%bh,eh%mask)

         if(eh%ph%unit_ovlp) then
            call elsi_allocate(eh%bh,eh%ovlp_real_sp,1,"ovlp_real_sp",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_real_sp,eh%bh%nnz_l_sp1,&
                 "ovlp_real_sp",caller)
         end if

         call elsi_allocate(eh%bh,eh%ham_real_sp,eh%bh%nnz_l_sp1,"ham_real_sp",&
              caller)
         call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lcol_sp1,eh%ph%n_states,&
              "evec_real",caller)
         call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,"row_ind_sp1",&
              caller)
         call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
              "col_ptr_sp1",caller)
      end if

      call elsi_blacs_to_sips_hs(eh%ph,eh%bh,ham,ovlp,eh%mask,eh%ham_real_sp,&
           eh%ovlp_real_sp,eh%row_ind_sp1,eh%col_ptr_sp1)
      call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
           eh%ham_real_sp,eh%ovlp_real_sp,eval,eh%evec_real)
      call elsi_sips_to_blacs_ev(eh%ph,eh%bh,eh%evec_real,evec)
   case(MAGMA_SOLVER)
      call elsi_init_magma(eh%ph)
      call elsi_solve_magma(eh%ph,eh%bh,ham,ovlp,eval,evec)
   case default
      write(msg,"(A)") "Unsupported eigensolver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! Compute the eigenvalues and eigenvectors of a Kohn-Sham problem. Note the
!! intent(inout), everything may be reused in the next call.
!!
subroutine elsi_ev_complex(eh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(inout) :: ham(eh%bh%n_lrow,eh%bh%n_lcol) !< Hamiltonian
   complex(kind=r8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   real(kind=r8), intent(inout) :: eval(eh%ph%n_basis) !< Eigenvalues
   complex(kind=r8), intent(out) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   real(kind=r8) :: t0
   character(len=29) :: dt0
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_ev_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)
   call elsi_decide_ev(eh%ph,eh%bh)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      if(eh%ph%parallel_mode == SINGLE_PROC) then
         call elsi_solve_lapack(eh%ph,eh%bh,ham,ovlp,eval,evec)
      else
         if(eh%ph%n_basis_c > 0) then
            if(.not. allocated(eh%ham_cmplx_v)) then
               call elsi_allocate(eh%bh,eh%ham_cmplx_v,eh%ph%n_lrow_v,&
                    eh%ph%n_lcol_v,"ham_cmplx_v",caller)
               call elsi_allocate(eh%bh,eh%ovlp_cmplx_v,eh%ph%n_lrow_v,&
                    eh%ph%n_lcol_v,"ovlp_cmplx_v",caller)
               call elsi_allocate(eh%bh,eh%evec_cmplx_v,eh%ph%n_lrow_v,&
                    eh%ph%n_lcol_v,"evec_cmplx_v",caller)
            end if

            call elsi_do_fc_elpa(eh%ph,eh%bh,ham,ovlp,evec,eh%perm_fc,&
                 eh%ham_cmplx_v,eh%ovlp_cmplx_v,eh%evec_cmplx_v)
            call elsi_init_elpa(eh%ph,eh%bh)
            call elsi_solve_elpa(eh%ph,eh%bh,eh%ham_cmplx_v,eh%ovlp_cmplx_v,&
                 eval(eh%ph%n_basis_c+1:eh%ph%n_basis_c+eh%ph%n_basis_v),&
                 eh%evec_cmplx_v)
            call elsi_undo_fc_elpa(eh%ph,eh%bh,ham,ovlp,evec,eh%perm_fc,&
                 eval(1:eh%ph%n_basis_c),eh%evec_cmplx_v)
         else
            call elsi_init_elpa(eh%ph,eh%bh)
            call elsi_solve_elpa(eh%ph,eh%bh,ham,ovlp,eval,evec)
         end if
      end if
   case(MAGMA_SOLVER)
      call elsi_init_magma(eh%ph)
      call elsi_solve_magma(eh%ph,eh%bh,ham,ovlp,eval,evec)
   case default
      write(msg,"(A)") "Unsupported eigensolver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! Compute the eigenvalues and eigenvectors of a Kohn-Sham problem. Note the
!! intent(inout), everything may be reused in the next call.
!!
subroutine elsi_ev_real_sparse(eh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(inout) :: ham(eh%bh%nnz_l_sp) !< Hamiltonian
   real(kind=r8), intent(inout) :: ovlp(eh%bh%nnz_l_sp) !< Overlap
   real(kind=r8), intent(inout) :: eval(eh%ph%n_basis) !< Eigenvalues
   real(kind=r8), intent(out) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   real(kind=r8) :: t0
   integer(kind=i4) :: solver
   character(len=29) :: dt0
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_ev_real_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)
   call elsi_decide_ev(eh%ph,eh%bh)

   eh%ph%n_calls = eh%ph%n_calls+1
   solver = eh%ph%solver

   if(eh%ph%solver == SIPS_SOLVER .and. eh%ph%n_calls <= eh%ph%sips_n_elpa) then
      solver = ELPA_SOLVER
   end if

   select case(solver)
   case(ELPA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)

      if(.not. allocated(eh%ham_real_den)) then
         call elsi_allocate(eh%bh,eh%ham_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "ham_real_den",caller)
      end if

      if(.not. allocated(eh%ovlp_real_den)) then
         if(.not. eh%ph%unit_ovlp) then
            call elsi_allocate(eh%bh,eh%ovlp_real_den,eh%bh%n_lrow,&
                 eh%bh%n_lcol,"ovlp_real_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_real_den,1,1,"ovlp_real_den",&
                 caller)
         end if
      end if

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%ham_real_den,eh%ovlp_real_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ham_real_den,eh%ovlp_real_den)
      case(GENERIC_COO)
         if(.not. allocated(eh%map_den)) then
            call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "map_den",caller)
         end if

         call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ham_real_den,eh%ovlp_real_den,eh%map_den)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_solve_elpa(eh%ph,eh%bh,eh%ham_real_den,eh%ovlp_real_den,eval,&
           evec)
   case(EIGENEXA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)
      call elsi_init_eigenexa(eh%ph,eh%bh)

      if(.not. allocated(eh%ham_real_den)) then
         call elsi_allocate(eh%bh,eh%ham_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "ham_real_den",caller)
      end if

      if(.not. allocated(eh%ovlp_real_den)) then
         if(.not. eh%ph%unit_ovlp) then
            call elsi_allocate(eh%bh,eh%ovlp_real_den,eh%bh%n_lrow,&
                 eh%bh%n_lcol,"ovlp_real_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_real_den,1,1,"ovlp_real_den",&
                 caller)
         end if
      end if

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%ham_real_den,eh%ovlp_real_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ham_real_den,eh%ovlp_real_den)
      case(GENERIC_COO)
         if(.not. allocated(eh%map_den)) then
            call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "map_den",caller)
         end if

         call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ham_real_den,eh%ovlp_real_den,eh%map_den)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_solve_eigenexa(eh%ph,eh%bh,eh%ham_real_den,eh%ovlp_real_den,&
           eval,evec)
   case(SIPS_SOLVER)
      call elsi_init_sips(eh%ph,eh%bh)

      if(allocated(eh%ham_real_den)) then
         call elsi_deallocate(eh%bh,eh%ham_real_den,"ham_real_den")
      end if

      if(allocated(eh%ovlp_real_den)) then
         call elsi_deallocate(eh%bh,eh%ovlp_real_den,"ovlp_real_den")
      end if

      if(.not. allocated(eh%evec_real)) then
         call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lcol_sp1,eh%ph%n_states,&
              "evec_real",caller)
      end if

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,ham,&
              ovlp,eval,eh%evec_real)
      case(SIESTA_CSC)
         if(.not. allocated(eh%row_ind_sp1)) then
            call elsi_siesta_to_sips_hs_dim(eh%ph,eh%bh,eh%col_ptr_sp2)

            if(eh%ph%unit_ovlp) then
               call elsi_allocate(eh%bh,eh%ovlp_real_sp,1,"ovlp_real_sp",caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_real_sp,eh%bh%nnz_l_sp1,&
                    "ovlp_real_sp",caller)
            end if

            call elsi_allocate(eh%bh,eh%ham_real_sp,eh%bh%nnz_l_sp1,&
                 "ham_real_sp",caller)
            call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                 "row_ind_sp1",caller)
            call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                 "col_ptr_sp1",caller)
         end if

         call elsi_siesta_to_sips_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ham_real_sp,eh%ovlp_real_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1)
         call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%ham_real_sp,eh%ovlp_real_sp,eval,eh%evec_real)
      case(GENERIC_COO)
         if(.not. allocated(eh%row_ind_sp1)) then
            call elsi_generic_to_sips_hs_dim(eh%ph,eh%bh,eh%col_ind_sp3)

            if(eh%ph%unit_ovlp) then
               call elsi_allocate(eh%bh,eh%ovlp_real_sp,1,"ovlp_real_sp",caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_real_sp,eh%bh%nnz_l_sp1,&
                    "ovlp_real_sp",caller)
            end if

            call elsi_allocate(eh%bh,eh%ham_real_sp,eh%bh%nnz_l_sp1,&
                 "ham_real_sp",caller)
            call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                 "row_ind_sp1",caller)
            call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                 "col_ptr_sp1",caller)
            call elsi_allocate(eh%bh,eh%map_sp1,eh%bh%nnz_l_sp1,"map_sp1",&
                 caller)
         end if

         call elsi_generic_to_sips_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ham_real_sp,eh%ovlp_real_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%map_sp1)
         call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%ham_real_sp,eh%ovlp_real_sp,eval,eh%evec_real)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_sips_to_blacs_ev(eh%ph,eh%bh,eh%evec_real,evec)
   case default
      write(msg,"(A)") "Unsupported eigensolver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! Compute the eigenvalues and eigenvectors of a Kohn-Sham problem. Note the
!! intent(inout), everything may be reused in the next call.
!!
subroutine elsi_ev_complex_sparse(eh,ham,ovlp,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(inout) :: ham(eh%bh%nnz_l_sp) !< Hamiltonian
   complex(kind=r8), intent(inout) :: ovlp(eh%bh%nnz_l_sp) !< Overlap
   real(kind=r8), intent(inout) :: eval(eh%ph%n_basis) !< Eigenvalues
   complex(kind=r8), intent(out) :: evec(eh%bh%n_lrow,eh%bh%n_lcol) !< Eigenvectors

   real(kind=r8) :: t0
   character(len=29) :: dt0
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_ev_complex_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)
   call elsi_decide_ev(eh%ph,eh%bh)

   eh%ph%n_calls = eh%ph%n_calls+1

   select case(eh%ph%solver)
   case(ELPA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)

      if(.not. allocated(eh%ham_cmplx_den)) then
         call elsi_allocate(eh%bh,eh%ham_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "ham_cmplx_den",caller)
      end if

      if(.not. allocated(eh%ovlp_cmplx_den)) then
         if(.not. eh%ph%unit_ovlp) then
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,eh%bh%n_lrow,&
                 eh%bh%n_lcol,"ovlp_cmplx_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,1,1,"ovlp_cmplx_den",&
                 caller)
         end if
      end if

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%ham_cmplx_den,eh%ovlp_cmplx_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ham_cmplx_den,eh%ovlp_cmplx_den)
      case(GENERIC_COO)
         if(.not. allocated(eh%map_den)) then
            call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "map_den",caller)
         end if

         call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ham_cmplx_den,eh%ovlp_cmplx_den,eh%map_den)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_solve_elpa(eh%ph,eh%bh,eh%ham_cmplx_den,eh%ovlp_cmplx_den,eval,&
           evec)
   case default
      write(msg,"(A)") "Unsupported eigensolver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! Compute the density matrix of a Kohn-Sham problem. Note the intent(inout),
!! everything may be reused in the next call.
!!
subroutine elsi_dm_real(eh,ham,ovlp,dm,ebs)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(inout) :: ham(eh%bh%n_lrow,eh%bh%n_lcol) !< Hamiltonian
   real(kind=r8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   real(kind=r8), intent(out) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix
   real(kind=r8), intent(out) :: ebs !< Band structure energy

   real(kind=r8) :: t0
   integer(kind=i4) :: solver
   character(len=29) :: dt0
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_dm_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)
   call elsi_decide_dm(eh%ph,eh%bh,ham)

   eh%ph%n_calls = eh%ph%n_calls+1
   eh%ph%ill_check = .false.
   solver = eh%ph%solver

   if(eh%ph%solver == SIPS_SOLVER .and. eh%ph%n_calls <= eh%ph%sips_n_elpa) then
      solver = ELPA_SOLVER
   end if

   if(eh%ph%solver == OMM_SOLVER .and. eh%ph%n_calls <= eh%ph%omm_n_elpa) then
      solver = ELPA_SOLVER
   end if

   if(eh%ph%solver /= solver .or. eh%ph%save_ovlp) then
      if(.not. allocated(eh%ovlp_real_copy)) then
         call elsi_allocate(eh%bh,eh%ovlp_real_copy,eh%bh%n_lrow,eh%bh%n_lcol,&
              "ovlp_real_copy",caller)
      end if

      if(eh%ph%n_calls == 1) then
         eh%ovlp_real_copy(:,:) = ovlp
      end if
   end if

   select case(solver)
   case(ELPA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)

      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      end if

      if(.not. allocated(eh%evec_real)) then
         call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lrow,eh%bh%n_lcol,&
              "evec_real",caller)
      end if

      if(.not. allocated(eh%occ)) then
         call elsi_allocate(eh%bh,eh%occ,eh%ph%n_states,eh%ph%n_spins,&
              eh%ph%n_kpts,"occ",caller)
      end if

      call elsi_solve_elpa(eh%ph,eh%bh,ham,ovlp,eh%eval,eh%evec_real)
      call elsi_get_occ_for_dm(eh%ph,eh%bh,eh%eval,eh%occ)
      call elsi_build_dm_edm(eh%ph,eh%bh,eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),&
           eh%evec_real,dm,GET_DM)
      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)

      eh%ph%eval_ready = .true.
      eh%ph%evec_ready = .true.
      eh%ph%occ_ready = .true.
   case(OMM_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)
      call elsi_init_omm(eh%ph,eh%bh)

      if(eh%ph%n_calls > 1 .and. eh%ph%n_calls == eh%ph%omm_n_elpa+1) then
         if(eh%ph%omm_flavor == 0) then
            ! Restore overlap
            ovlp(:,:) = eh%ovlp_real_copy
         end if
      end if

      if(.not. allocated(eh%omm_c_real)) then
         call elsi_allocate(eh%bh,eh%omm_c_real,eh%ph%omm_n_lrow,eh%bh%n_lcol,&
              "omm_c_real",caller)

         ! Initialize coefficient matrix with ELPA eigenvectors
         call pdtran(eh%ph%n_basis,eh%ph%n_basis,1.0_r8,eh%evec_real,1,1,&
              eh%bh%desc,0.0_r8,dm,1,1,eh%bh%desc)

         eh%omm_c_real(1:eh%ph%omm_n_lrow,:) = dm(1:eh%ph%omm_n_lrow,:)

         if(allocated(eh%evec_real)) then
            call elsi_deallocate(eh%bh,eh%evec_real,"evec_real")
         end if

         if(allocated(eh%eval)) then
            call elsi_deallocate(eh%bh,eh%eval,"eval")
         end if

         if(allocated(eh%occ)) then
            call elsi_deallocate(eh%bh,eh%occ,"occ")
         end if
      end if

      call elsi_solve_omm(eh%ph,eh%bh,ham,ovlp,eh%omm_c_real,dm)
      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(eh%ph,eh%bh)

      if(.not. allocated(eh%row_ind_sp1)) then
         call elsi_allocate(eh%bh,eh%mask,eh%bh%n_lrow,eh%bh%n_lcol,"mask",&
              caller)

         call elsi_blacs_to_mask(eh%ph,eh%bh,ham,ovlp,eh%mask)
         call elsi_blacs_to_pexsi_hs_dim(eh%ph,eh%bh,eh%mask)

         if(eh%ph%unit_ovlp) then
            call elsi_allocate(eh%bh,eh%ovlp_real_sp,1,"ovlp_real_sp",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_real_sp,eh%bh%nnz_l_sp1,&
                 "ovlp_real_sp",caller)
         end if

         call elsi_allocate(eh%bh,eh%ham_real_sp,eh%bh%nnz_l_sp1,"ham_real_sp",&
              caller)
         call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,"row_ind_sp1",&
              caller)
         call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
              "col_ptr_sp1",caller)
      end if

      call elsi_blacs_to_pexsi_hs(eh%ph,eh%bh,ham,ovlp,eh%mask,eh%ham_real_sp,&
           eh%ovlp_real_sp,eh%row_ind_sp1,eh%col_ptr_sp1)

      if(.not. allocated(eh%pexsi_ne_vec)) then
         call elsi_allocate(eh%bh,eh%pexsi_ne_vec,eh%ph%pexsi_options%nPoints,&
              "pexsi_ne_vec",caller)
      end if

      if(.not. allocated(eh%dm_real_sp)) then
         call elsi_allocate(eh%bh,eh%dm_real_sp,eh%bh%nnz_l_sp1,"dm_real_sp",&
              caller)
      end if

      eh%dm_real_sp(:) = 0.0_r8

      call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
           eh%pexsi_ne_vec,eh%ham_real_sp,eh%ovlp_real_sp,eh%dm_real_sp)
      call elsi_pexsi_to_blacs_dm(eh%ph,eh%bh,eh%dm_real_sp,eh%row_ind_sp1,&
           eh%col_ptr_sp1,dm)
      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
   case(EIGENEXA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)
      call elsi_init_eigenexa(eh%ph,eh%bh)

      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      end if

      if(.not. allocated(eh%evec_real)) then
         call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lrow,eh%bh%n_lcol,&
              "evec_real",caller)
      end if

      if(.not. allocated(eh%occ)) then
         call elsi_allocate(eh%bh,eh%occ,eh%ph%n_states,eh%ph%n_spins,&
              eh%ph%n_kpts,"occ",caller)
      end if

      call elsi_solve_eigenexa(eh%ph,eh%bh,ham,ovlp,eh%eval,eh%evec_real)
      call elsi_get_occ_for_dm(eh%ph,eh%bh,eh%eval,eh%occ)
      call elsi_build_dm_edm(eh%ph,eh%bh,eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),&
           eh%evec_real,dm,GET_DM)
      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)

      eh%ph%eval_ready = .true.
      eh%ph%evec_ready = .true.
      eh%ph%occ_ready = .true.
   case(SIPS_SOLVER)
      if(eh%ph%n_calls > 1 .and. eh%ph%n_calls == eh%ph%sips_n_elpa+1) then
         ! Restore overlap
         ovlp(:,:) = eh%ovlp_real_copy

         call elsi_deallocate(eh%bh,eh%evec_real,"evec_real")
      end if

      call elsi_init_sips(eh%ph,eh%bh)

      if(.not. allocated(eh%row_ind_sp1)) then
         call elsi_allocate(eh%bh,eh%mask,eh%bh%n_lrow,eh%bh%n_lcol,"mask",&
              caller)

         call elsi_blacs_to_mask(eh%ph,eh%bh,ham,ovlp,eh%mask)
         call elsi_blacs_to_sips_hs_dim(eh%ph,eh%bh,eh%mask)

         if(eh%ph%unit_ovlp) then
            call elsi_allocate(eh%bh,eh%ovlp_real_sp,1,"ovlp_real_sp",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_real_sp,eh%bh%nnz_l_sp1,&
                 "ovlp_real_sp",caller)
         end if

         call elsi_allocate(eh%bh,eh%ham_real_sp,eh%bh%nnz_l_sp1,"ham_real_sp",&
              caller)
         call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,"row_ind_sp1",&
              caller)
         call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
              "col_ptr_sp1",caller)
      end if

      call elsi_blacs_to_sips_hs(eh%ph,eh%bh,ham,ovlp,eh%mask,eh%ham_real_sp,&
           eh%ovlp_real_sp,eh%row_ind_sp1,eh%col_ptr_sp1)

      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      end if

      if(.not. allocated(eh%evec_real)) then
         call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lcol_sp1,eh%ph%n_states,&
              "evec_real",caller)
      end if

      if(.not. allocated(eh%occ)) then
         call elsi_allocate(eh%bh,eh%occ,eh%ph%n_states,eh%ph%n_spins,&
              eh%ph%n_kpts,"occ",caller)
      end if

      if(.not. allocated(eh%dm_real_sp)) then
         call elsi_allocate(eh%bh,eh%dm_real_sp,eh%bh%nnz_l_sp1,"dm_real_sp",&
              caller)
      end if

      eh%dm_real_sp(:) = 0.0_r8

      call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
           eh%ham_real_sp,eh%ovlp_real_sp,eh%eval,eh%evec_real)
      call elsi_get_occ_for_dm(eh%ph,eh%bh,eh%eval,eh%occ)
      call elsi_build_dm_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
           eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),eh%dm_real_sp,GET_DM)
      call elsi_sips_to_blacs_dm(eh%ph,eh%bh,eh%dm_real_sp,eh%row_ind_sp1,&
           eh%col_ptr_sp1,dm)
      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)

      eh%ph%eval_ready = .true.
      eh%ph%occ_ready = .true.
   case(NTPOLY_SOLVER)
      call elsi_init_ntpoly(eh%ph,eh%bh)

      if(.not. allocated(eh%mask)) then
         call elsi_allocate(eh%bh,eh%mask,eh%bh%n_lrow,eh%bh%n_lcol,"mask",&
              caller)

         call elsi_blacs_to_mask(eh%ph,eh%bh,ham,ovlp,eh%mask)
      end if

      call elsi_blacs_to_ntpoly_hs(eh%ph,eh%bh,ham,ovlp,eh%mask,eh%nt_ham,&
           eh%nt_ovlp)
      call elsi_solve_ntpoly(eh%ph,eh%bh,eh%nt_ham,eh%nt_ovlp,eh%nt_dm)
      call elsi_ntpoly_to_blacs_dm(eh%bh,eh%nt_dm,dm)
      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
   case default
      write(msg,"(A)") "Unsupported density matrix solver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%edm_ready = .true.

   if(eh%ph%save_ovlp) then
      if(.not. allocated(eh%dm_real_copy)) then
         call elsi_allocate(eh%bh,eh%dm_real_copy,eh%bh%n_lrow,eh%bh%n_lcol,&
              "dm_real_copy",caller)
      end if

      eh%dm_real_copy(:,:) = dm
   end if

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! Compute the density matrix of a Kohn-Sham problem. Note the intent(inout),
!! everything may be reused in the next call.
!!
subroutine elsi_dm_complex(eh,ham,ovlp,dm,ebs)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(inout) :: ham(eh%bh%n_lrow,eh%bh%n_lcol) !< Hamiltonian
   complex(kind=r8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol) !< Overlap
   complex(kind=r8), intent(out) :: dm(eh%bh%n_lrow,eh%bh%n_lcol) !< Density matrix
   real(kind=r8), intent(out) :: ebs !< Band structure energy

   real(kind=r8) :: t0
   integer(kind=i4) :: solver
   character(len=29) :: dt0
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_dm_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)
   call elsi_decide_dm(eh%ph,eh%bh,ham)

   eh%ph%n_calls = eh%ph%n_calls+1
   eh%ph%ill_check = .false.
   solver = eh%ph%solver

   if(eh%ph%solver == OMM_SOLVER .and. eh%ph%n_calls <= eh%ph%omm_n_elpa) then
      solver = ELPA_SOLVER
   end if

   if(eh%ph%solver /= solver .or. eh%ph%save_ovlp) then
      if(.not. allocated(eh%ovlp_cmplx_copy)) then
         call elsi_allocate(eh%bh,eh%ovlp_cmplx_copy,eh%bh%n_lrow,eh%bh%n_lcol,&
              "ovlp_cmplx_copy",caller)
      end if

      if(eh%ph%n_calls == 1) then
         eh%ovlp_cmplx_copy(:,:) = ovlp
      end if
   end if

   select case(solver)
   case(ELPA_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)

      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      end if

      if(.not. allocated(eh%evec_cmplx)) then
         call elsi_allocate(eh%bh,eh%evec_cmplx,eh%bh%n_lrow,eh%bh%n_lcol,&
              "evec_cmplx",caller)
      end if

      if(.not. allocated(eh%occ)) then
         call elsi_allocate(eh%bh,eh%occ,eh%ph%n_states,eh%ph%n_spins,&
              eh%ph%n_kpts,"occ",caller)
      end if

      call elsi_solve_elpa(eh%ph,eh%bh,ham,ovlp,eh%eval,eh%evec_cmplx)
      call elsi_get_occ_for_dm(eh%ph,eh%bh,eh%eval,eh%occ)
      call elsi_build_dm_edm(eh%ph,eh%bh,eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),&
           eh%evec_cmplx,dm,GET_DM)
      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)

      eh%ph%eval_ready = .true.
      eh%ph%evec_ready = .true.
      eh%ph%occ_ready = .true.
   case(OMM_SOLVER)
      call elsi_init_elpa(eh%ph,eh%bh)
      call elsi_init_omm(eh%ph,eh%bh)

      if(eh%ph%n_calls > 1 .and. eh%ph%n_calls == eh%ph%omm_n_elpa+1) then
         if(eh%ph%omm_flavor == 0) then
            ! Restore overlap
            ovlp(:,:) = eh%ovlp_cmplx_copy
         end if
      end if

      if(.not. allocated(eh%omm_c_cmplx)) then
         call elsi_allocate(eh%bh,eh%omm_c_cmplx,eh%ph%omm_n_lrow,eh%bh%n_lcol,&
              "omm_c_cmplx",caller)

         ! Initialize coefficient matrix with ELPA eigenvectors
         call pztranc(eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),&
              eh%evec_cmplx,1,1,eh%bh%desc,(0.0_r8,0.0_r8),dm,1,1,eh%bh%desc)

         eh%omm_c_cmplx(1:eh%ph%omm_n_lrow,:) = dm(1:eh%ph%omm_n_lrow,:)

         if(allocated(eh%evec_cmplx)) then
            call elsi_deallocate(eh%bh,eh%evec_cmplx,"evec_cmplx")
         end if

         if(allocated(eh%eval)) then
            call elsi_deallocate(eh%bh,eh%eval,"eval")
         end if

         if(allocated(eh%occ)) then
            call elsi_deallocate(eh%bh,eh%occ,"occ")
         end if
      end if

      call elsi_solve_omm(eh%ph,eh%bh,ham,ovlp,eh%omm_c_cmplx,dm)
      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(eh%ph,eh%bh)

      if(.not. allocated(eh%row_ind_sp1)) then
         call elsi_allocate(eh%bh,eh%mask,eh%bh%n_lrow,eh%bh%n_lcol,"mask",&
              caller)

         call elsi_blacs_to_mask(eh%ph,eh%bh,ham,ovlp,eh%mask)
         call elsi_blacs_to_pexsi_hs_dim(eh%ph,eh%bh,eh%mask)

         if(eh%ph%unit_ovlp) then
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_sp,1,"ovlp_cmplx_sp",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_sp,eh%bh%nnz_l_sp1,&
                 "ovlp_cmplx_sp",caller)
         end if

         call elsi_allocate(eh%bh,eh%ham_cmplx_sp,eh%bh%nnz_l_sp1,&
              "ham_cmplx_sp",caller)
         call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,"row_ind_sp1",&
              caller)
         call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
              "col_ptr_sp1",caller)
      end if

      call elsi_blacs_to_pexsi_hs(eh%ph,eh%bh,ham,ovlp,eh%mask,eh%ham_cmplx_sp,&
           eh%ovlp_cmplx_sp,eh%row_ind_sp1,eh%col_ptr_sp1)

      if(.not. allocated(eh%pexsi_ne_vec)) then
         call elsi_allocate(eh%bh,eh%pexsi_ne_vec,eh%ph%pexsi_options%nPoints,&
              "pexsi_ne_vec",caller)
      end if

      if(.not. allocated(eh%dm_cmplx_sp)) then
         call elsi_allocate(eh%bh,eh%dm_cmplx_sp,eh%bh%nnz_l_sp1,"dm_cmplx_sp",&
              caller)
      end if

      eh%dm_cmplx_sp(:) = (0.0_r8,0.0_r8)

      call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
           eh%pexsi_ne_vec,eh%ham_cmplx_sp,eh%ovlp_cmplx_sp,eh%dm_cmplx_sp)
      call elsi_pexsi_to_blacs_dm(eh%ph,eh%bh,eh%dm_cmplx_sp,eh%row_ind_sp1,&
           eh%col_ptr_sp1,dm)
      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
   case(NTPOLY_SOLVER)
      call elsi_init_ntpoly(eh%ph,eh%bh)

      if(.not. allocated(eh%mask)) then
         call elsi_allocate(eh%bh,eh%mask,eh%bh%n_lrow,eh%bh%n_lcol,"mask",&
              caller)

         call elsi_blacs_to_mask(eh%ph,eh%bh,ham,ovlp,eh%mask)
      end if

      call elsi_blacs_to_ntpoly_hs(eh%ph,eh%bh,ham,ovlp,eh%mask,eh%nt_ham,&
           eh%nt_ovlp)
      call elsi_solve_ntpoly(eh%ph,eh%bh,eh%nt_ham,eh%nt_ovlp,eh%nt_dm)
      call elsi_ntpoly_to_blacs_dm(eh%bh,eh%nt_dm,dm)
      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
   case default
      write(msg,"(A)") "Unsupported density matrix solver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%edm_ready = .true.

   if(eh%ph%save_ovlp) then
      if(.not. allocated(eh%dm_real_copy)) then
         call elsi_allocate(eh%bh,eh%dm_real_copy,eh%bh%n_lrow,eh%bh%n_lcol,&
              "dm_real_copy",caller)
      end if

      eh%dm_real_copy(:,:) = dm
   end if

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! Compute the density matrix of a Kohn-Sham problem. Note the intent(inout),
!! everything may be reused in the next call.
!!
subroutine elsi_dm_real_sparse(eh,ham,ovlp,dm,ebs)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(inout) :: ham(eh%bh%nnz_l_sp) !< Hamiltonian
   real(kind=r8), intent(inout) :: ovlp(eh%bh%nnz_l_sp) !< Overlap
   real(kind=r8), intent(out) :: dm(eh%bh%nnz_l_sp) !< Density matrix
   real(kind=r8), intent(out) :: ebs !< Band structure energy

   real(kind=r8) :: t0
   integer(kind=i4) :: solver
   character(len=29) :: dt0
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_dm_real_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)
   call elsi_decide_dm(eh%ph,eh%bh)

   eh%ph%n_calls = eh%ph%n_calls+1
   eh%ph%ill_check = .false.
   solver = eh%ph%solver

   if(eh%ph%solver == SIPS_SOLVER .and. eh%ph%n_calls <= eh%ph%sips_n_elpa) then
      solver = ELPA_SOLVER
   end if

   if(eh%ph%solver == OMM_SOLVER .and. eh%ph%n_calls <= eh%ph%omm_n_elpa) then
      solver = ELPA_SOLVER
   end if

   select case(solver)
   case(ELPA_SOLVER)
      call elsi_init_blacs(eh)
      call elsi_init_elpa(eh%ph,eh%bh)

      if(.not. allocated(eh%ham_real_den)) then
         call elsi_allocate(eh%bh,eh%ham_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "ham_real_den",caller)
      end if

      if(.not. allocated(eh%ovlp_real_den)) then
         if(.not. eh%ph%unit_ovlp) then
            call elsi_allocate(eh%bh,eh%ovlp_real_den,eh%bh%n_lrow,&
                 eh%bh%n_lcol,"ovlp_real_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_real_den,1,1,"ovlp_real_den",&
                 caller)
         end if
      end if

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%ham_real_den,eh%ovlp_real_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ham_real_den,eh%ovlp_real_den)
      case(GENERIC_COO)
         if(.not. allocated(eh%map_den)) then
            call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "map_den",caller)
         end if

         call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ham_real_den,eh%ovlp_real_den,eh%map_den)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      end if

      if(.not. allocated(eh%evec_real)) then
         call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lrow,eh%bh%n_lcol,&
              "evec_real",caller)
      end if

      if(.not. allocated(eh%dm_real_den)) then
         call elsi_allocate(eh%bh,eh%dm_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "dm_real_den",caller)
      end if

      if(.not. allocated(eh%occ)) then
         call elsi_allocate(eh%bh,eh%occ,eh%ph%n_states,eh%ph%n_spins,&
              eh%ph%n_kpts,"occ",caller)
      end if

      if(eh%ph%solver == OMM_SOLVER .and. eh%ph%omm_flavor == 0) then
         if(.not. allocated(eh%ovlp_real_copy)) then
            call elsi_allocate(eh%bh,eh%ovlp_real_copy,eh%bh%n_lrow,&
                 eh%bh%n_lcol,"ovlp_real_copy",caller)
         end if

         if(eh%ph%n_calls == 1) then
            eh%ovlp_real_copy(:,:) = eh%ovlp_real_den
         end if
      end if

      call elsi_solve_elpa(eh%ph,eh%bh,eh%ham_real_den,eh%ovlp_real_den,&
           eh%eval,eh%evec_real)
      call elsi_get_occ_for_dm(eh%ph,eh%bh,eh%eval,eh%occ)
      call elsi_build_dm_edm(eh%ph,eh%bh,eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),&
           eh%evec_real,eh%dm_real_den,GET_DM)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%dm_real_den,dm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%dm_real_den,dm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_blacs_to_generic_dm(eh%ph,eh%bh,eh%dm_real_den,eh%map_den,&
              dm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)

      eh%ph%eval_ready = .true.
      eh%ph%evec_ready = .true.
      eh%ph%occ_ready = .true.
   case(OMM_SOLVER)
      call elsi_init_blacs(eh)
      call elsi_init_elpa(eh%ph,eh%bh)
      call elsi_init_omm(eh%ph,eh%bh)

      if(.not. allocated(eh%ham_real_den)) then
         call elsi_allocate(eh%bh,eh%ham_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "ham_real_den",caller)
      end if

      if(.not. allocated(eh%ovlp_real_den)) then
         if(.not. eh%ph%unit_ovlp) then
            call elsi_allocate(eh%bh,eh%ovlp_real_den,eh%bh%n_lrow,&
                 eh%bh%n_lcol,"ovlp_real_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_real_den,1,1,"ovlp_real_den",&
                 caller)
         end if
      end if

      if(.not. allocated(eh%dm_real_den)) then
         call elsi_allocate(eh%bh,eh%dm_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "dm_real_den",caller)
      end if

      if(.not. allocated(eh%omm_c_real)) then
         call elsi_allocate(eh%bh,eh%omm_c_real,eh%ph%omm_n_lrow,eh%bh%n_lcol,&
              "omm_c_real",caller)

         ! Initialize coefficient matrix with ELPA eigenvectors
         call pdtran(eh%ph%n_basis,eh%ph%n_basis,1.0_r8,eh%evec_real,1,1,&
              eh%bh%desc,0.0_r8,eh%dm_real_den,1,1,eh%bh%desc)

         eh%omm_c_real(1:eh%ph%omm_n_lrow,:)&
            = eh%dm_real_den(1:eh%ph%omm_n_lrow,:)

         if(allocated(eh%evec_real)) then
            call elsi_deallocate(eh%bh,eh%evec_real,"evec_real")
         end if

         if(allocated(eh%eval)) then
            call elsi_deallocate(eh%bh,eh%eval,"eval")
         end if

         if(allocated(eh%occ)) then
            call elsi_deallocate(eh%bh,eh%occ,"occ")
         end if
      end if

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%ham_real_den,eh%ovlp_real_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ham_real_den,eh%ovlp_real_den)
      case(GENERIC_COO)
         if(.not. allocated(eh%map_den)) then
            call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "map_den",caller)
         end if

         call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ham_real_den,eh%ovlp_real_den,eh%map_den)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      if(eh%ph%n_calls > 1 .and. eh%ph%n_calls == eh%ph%omm_n_elpa+1) then
         if(eh%ph%omm_flavor == 0) then
            ! Restore overlap
            eh%ovlp_real_den(:,:) = eh%ovlp_real_copy
         end if
      end if

      call elsi_solve_omm(eh%ph,eh%bh,eh%ham_real_den,eh%ovlp_real_den,&
           eh%omm_c_real,eh%dm_real_den)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%dm_real_den,dm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%dm_real_den,dm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_blacs_to_generic_dm(eh%ph,eh%bh,eh%dm_real_den,eh%map_den,&
              dm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(eh%ph,eh%bh)

      if(.not. allocated(eh%pexsi_ne_vec)) then
         call elsi_allocate(eh%bh,eh%pexsi_ne_vec,eh%ph%pexsi_options%nPoints,&
              "pexsi_ne_vec",caller)
      end if

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%pexsi_ne_vec,ham,ovlp,dm)
      case(SIESTA_CSC)
         if(.not. allocated(eh%row_ind_sp1)) then
            call elsi_siesta_to_pexsi_hs_dim(eh%ph,eh%bh,eh%col_ptr_sp2)

            if(eh%ph%unit_ovlp) then
               call elsi_allocate(eh%bh,eh%ovlp_real_sp,1,"ovlp_real_sp",caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_real_sp,eh%bh%nnz_l_sp1,&
                    "ovlp_real_sp",caller)
            end if

            call elsi_allocate(eh%bh,eh%ham_real_sp,eh%bh%nnz_l_sp1,&
                 "ham_real_sp",caller)
            call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                 "row_ind_sp1",caller)
            call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                 "col_ptr_sp1",caller)
         end if

         call elsi_siesta_to_pexsi_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ham_real_sp,eh%ovlp_real_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1)

         if(.not. allocated(eh%dm_real_sp)) then
            call elsi_allocate(eh%bh,eh%dm_real_sp,eh%bh%nnz_l_sp1,&
                 "dm_real_sp",caller)
         end if

         eh%dm_real_sp(:) = 0.0_r8

         call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%pexsi_ne_vec,eh%ham_real_sp,eh%ovlp_real_sp,eh%dm_real_sp)
         call elsi_pexsi_to_siesta_dm(eh%ph,eh%bh,eh%dm_real_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,dm,eh%row_ind_sp2,eh%col_ptr_sp2)
      case(GENERIC_COO)
         if(.not. allocated(eh%row_ind_sp1)) then
            call elsi_generic_to_pexsi_hs_dim(eh%ph,eh%bh,eh%col_ind_sp3)

            if(eh%ph%unit_ovlp) then
               call elsi_allocate(eh%bh,eh%ovlp_real_sp,1,"ovlp_real_sp",caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_real_sp,eh%bh%nnz_l_sp1,&
                    "ovlp_real_sp",caller)
            end if

            call elsi_allocate(eh%bh,eh%ham_real_sp,eh%bh%nnz_l_sp1,&
                 "ham_real_sp",caller)
            call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                 "row_ind_sp1",caller)
            call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                 "col_ptr_sp1",caller)
            call elsi_allocate(eh%bh,eh%map_sp1,eh%bh%nnz_l_sp1,"map_sp1",&
                 caller)
         end if

         call elsi_generic_to_pexsi_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ham_real_sp,eh%ovlp_real_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%map_sp1)

         if(.not. allocated(eh%dm_real_sp)) then
            call elsi_allocate(eh%bh,eh%dm_real_sp,eh%bh%nnz_l_sp1,&
                 "dm_real_sp",caller)
         end if

         call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%pexsi_ne_vec,eh%ham_real_sp,eh%ovlp_real_sp,eh%dm_real_sp)
         call elsi_pexsi_to_generic_dm(eh%ph,eh%bh,eh%dm_real_sp,&
              eh%row_ind_sp1,eh%col_ptr_sp1,eh%map_sp1,dm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
   case(EIGENEXA_SOLVER)
      call elsi_init_blacs(eh)
      call elsi_init_elpa(eh%ph,eh%bh)
      call elsi_init_eigenexa(eh%ph,eh%bh)

      if(.not. allocated(eh%ham_real_den)) then
         call elsi_allocate(eh%bh,eh%ham_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "ham_real_den",caller)
      end if

      if(.not. allocated(eh%ovlp_real_den)) then
         if(.not. eh%ph%unit_ovlp) then
            call elsi_allocate(eh%bh,eh%ovlp_real_den,eh%bh%n_lrow,&
                 eh%bh%n_lcol,"ovlp_real_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_real_den,1,1,"ovlp_real_den",&
                 caller)
         end if
      end if

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%ham_real_den,eh%ovlp_real_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ham_real_den,eh%ovlp_real_den)
      case(GENERIC_COO)
         if(.not. allocated(eh%map_den)) then
            call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "map_den",caller)
         end if

         call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ham_real_den,eh%ovlp_real_den,eh%map_den)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      end if

      if(.not. allocated(eh%evec_real)) then
         call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lrow,eh%bh%n_lcol,&
              "evec_real",caller)
      end if

      if(.not. allocated(eh%dm_real_den)) then
         call elsi_allocate(eh%bh,eh%dm_real_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "dm_real_den",caller)
      end if

      if(.not. allocated(eh%occ)) then
         call elsi_allocate(eh%bh,eh%occ,eh%ph%n_states,eh%ph%n_spins,&
              eh%ph%n_kpts,"occ",caller)
      end if

      call elsi_solve_eigenexa(eh%ph,eh%bh,eh%ham_real_den,eh%ovlp_real_den,&
           eh%eval,eh%evec_real)
      call elsi_get_occ_for_dm(eh%ph,eh%bh,eh%eval,eh%occ)
      call elsi_build_dm_edm(eh%ph,eh%bh,eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),&
           eh%evec_real,eh%dm_real_den,GET_DM)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%dm_real_den,dm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%dm_real_den,dm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_blacs_to_generic_dm(eh%ph,eh%bh,eh%dm_real_den,eh%map_den,&
              dm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)

      eh%ph%eval_ready = .true.
      eh%ph%evec_ready = .true.
      eh%ph%occ_ready = .true.
   case(SIPS_SOLVER)
      call elsi_init_sips(eh%ph,eh%bh)

      if(allocated(eh%ham_real_den)) then
         call elsi_deallocate(eh%bh,eh%ham_real_den,"ham_real_den")
      end if

      if(allocated(eh%ovlp_real_den)) then
         call elsi_deallocate(eh%bh,eh%ovlp_real_den,"ovlp_real_den")
      end if

      if(allocated(eh%dm_real_den)) then
         call elsi_deallocate(eh%bh,eh%dm_real_den,"dm_real_den")
      end if

      if(.not. allocated(eh%occ)) then
         call elsi_allocate(eh%bh,eh%occ,eh%ph%n_states,eh%ph%n_spins,&
              eh%ph%n_kpts,"occ",caller)
      end if

      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      end if

      if(.not. allocated(eh%evec_real)) then
         call elsi_allocate(eh%bh,eh%evec_real,eh%bh%n_lcol_sp1,&
              eh%ph%n_states,"evec_real",caller)
      end if

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,ham,&
              ovlp,eh%eval,eh%evec_real)
         call elsi_get_occ_for_dm(eh%ph,eh%bh,eh%eval,eh%occ)
         call elsi_build_dm_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),dm,GET_DM)
         call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
      case(SIESTA_CSC)
         if(.not. allocated(eh%row_ind_sp1)) then
            call elsi_siesta_to_sips_hs_dim(eh%ph,eh%bh,eh%col_ptr_sp2)

            if(eh%ph%unit_ovlp) then
               call elsi_allocate(eh%bh,eh%ovlp_real_sp,1,"ovlp_real_sp",caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_real_sp,eh%bh%nnz_l_sp1,&
                    "ovlp_real_sp",caller)
            end if

            call elsi_allocate(eh%bh,eh%ham_real_sp,eh%bh%nnz_l_sp1,&
                 "ham_real_sp",caller)
            call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                 "row_ind_sp1",caller)
            call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                 "col_ptr_sp1",caller)
         end if

         call elsi_siesta_to_sips_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ham_real_sp,eh%ovlp_real_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1)

         if(.not. allocated(eh%dm_real_sp)) then
            call elsi_allocate(eh%bh,eh%dm_real_sp,eh%bh%nnz_l_sp1,&
                 "dm_real_sp",caller)
         end if

         eh%dm_real_sp(:) = 0.0_r8

         call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%ham_real_sp,eh%ovlp_real_sp,eh%eval,eh%evec_real)
         call elsi_get_occ_for_dm(eh%ph,eh%bh,eh%eval,eh%occ)
         call elsi_build_dm_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),eh%dm_real_sp,GET_DM)
         call elsi_sips_to_siesta_dm(eh%ph,eh%bh,eh%dm_real_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,dm,eh%row_ind_sp2,eh%col_ptr_sp2)
         call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
      case(GENERIC_COO)
         if(.not. allocated(eh%row_ind_sp1)) then
            call elsi_generic_to_sips_hs_dim(eh%ph,eh%bh,eh%col_ind_sp3)

            if(eh%ph%unit_ovlp) then
               call elsi_allocate(eh%bh,eh%ovlp_real_sp,1,"ovlp_real_sp",caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_real_sp,eh%bh%nnz_l_sp1,&
                    "ovlp_real_sp",caller)
            end if

            call elsi_allocate(eh%bh,eh%ham_real_sp,eh%bh%nnz_l_sp1,&
                 "ham_real_sp",caller)
            call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                 "row_ind_sp1",caller)
            call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                 "col_ptr_sp1",caller)
            call elsi_allocate(eh%bh,eh%map_sp1,eh%bh%nnz_l_sp1,"map_sp1",&
                 caller)
         end if

         call elsi_generic_to_sips_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ham_real_sp,eh%ovlp_real_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%map_sp1)

         if(.not. allocated(eh%dm_real_sp)) then
            call elsi_allocate(eh%bh,eh%dm_real_sp,eh%bh%nnz_l_sp1,&
                 "dm_real_sp",caller)
         end if

         call elsi_solve_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%ham_real_sp,eh%ovlp_real_sp,eh%eval,eh%evec_real)
         call elsi_get_occ_for_dm(eh%ph,eh%bh,eh%eval,eh%occ)
         call elsi_build_dm_edm_sips(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),eh%dm_real_sp,GET_DM)
         call elsi_sips_to_generic_dm(eh%ph,eh%bh,eh%dm_real_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%map_sp1,dm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      eh%ph%eval_ready = .true.
      eh%ph%occ_ready = .true.
   case(NTPOLY_SOLVER)
      call elsi_init_ntpoly(eh%ph,eh%bh)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_ntpoly_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%nt_ham,eh%nt_ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_ntpoly_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%nt_ham,eh%nt_ovlp)
      case(GENERIC_COO)
         call elsi_generic_to_ntpoly_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%nt_ham,eh%nt_ovlp,eh%nt_map)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_solve_ntpoly(eh%ph,eh%bh,eh%nt_ham,eh%nt_ovlp,eh%nt_dm)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_ntpoly_to_sips_dm(eh%ph,eh%bh,eh%nt_dm,dm,eh%row_ind_sp1,&
              eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_ntpoly_to_siesta_dm(eh%bh,eh%nt_dm,dm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_ntpoly_to_generic_dm(eh%ph,eh%bh,eh%nt_dm,eh%nt_map,dm,&
              eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
   case default
      write(msg,"(A)") "Unsupported density matrix solver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%edm_ready = .true.

   if(eh%ph%save_ovlp) then
      call elsi_init_ntpoly(eh%ph,eh%bh)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         if(eh%ph%n_calls == 1) then
            eh%ph%first_sips_to_ntpoly = .true.
         end if

         call elsi_sips_to_ntpoly_hs(eh%ph,eh%bh,dm,ovlp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%nt_dm_copy,eh%nt_ovlp_copy)
      case(SIESTA_CSC)
         if(eh%ph%n_calls == 1) then
            eh%ph%first_siesta_to_ntpoly = .true.
         end if

         call elsi_siesta_to_ntpoly_hs(eh%ph,eh%bh,dm,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%nt_dm_copy,eh%nt_ovlp_copy)
      case(GENERIC_COO)
         if(eh%ph%n_calls == 1) then
            eh%ph%first_generic_to_ntpoly = .true.
         end if

         call elsi_generic_to_ntpoly_hs(eh%ph,eh%bh,dm,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%nt_dm_copy,eh%nt_ovlp_copy,eh%nt_map)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select
   end if

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! Compute the density matrix of a Kohn-Sham problem. Note the intent(inout),
!! everything may be reused in the next call.
!!
subroutine elsi_dm_complex_sparse(eh,ham,ovlp,dm,ebs)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(inout) :: ham(eh%bh%nnz_l_sp) !< Hamiltonian
   complex(kind=r8), intent(inout) :: ovlp(eh%bh%nnz_l_sp) !< Overlap
   complex(kind=r8), intent(out) :: dm(eh%bh%nnz_l_sp) !< Density matrix
   real(kind=r8), intent(out) :: ebs !< Band structure energy

   real(kind=r8) :: t0
   integer(kind=i4) :: solver
   character(len=29) :: dt0
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_dm_complex_sparse"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)
   call elsi_decide_dm(eh%ph,eh%bh)

   eh%ph%n_calls = eh%ph%n_calls+1
   eh%ph%ill_check = .false.
   solver = eh%ph%solver

   if(eh%ph%solver == OMM_SOLVER .and. eh%ph%n_calls <= eh%ph%omm_n_elpa) then
      solver = ELPA_SOLVER
   end if

   select case(solver)
   case(ELPA_SOLVER)
      call elsi_init_blacs(eh)
      call elsi_init_elpa(eh%ph,eh%bh)

      if(.not. allocated(eh%ham_cmplx_den)) then
         call elsi_allocate(eh%bh,eh%ham_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "ham_cmplx_den",caller)
      end if

      if(.not. allocated(eh%ovlp_cmplx_den)) then
         if(.not. eh%ph%unit_ovlp) then
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,eh%bh%n_lrow,&
                 eh%bh%n_lcol,"ovlp_cmplx_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,1,1,"ovlp_cmplx_den",&
                 caller)
         end if
      end if

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%ham_cmplx_den,eh%ovlp_cmplx_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ham_cmplx_den,eh%ovlp_cmplx_den)
      case(GENERIC_COO)
         if(.not. allocated(eh%map_den)) then
            call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "map_den",caller)
         end if

         call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ham_cmplx_den,eh%ovlp_cmplx_den,eh%map_den)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      if(.not. allocated(eh%eval)) then
         call elsi_allocate(eh%bh,eh%eval,eh%ph%n_basis,"eval",caller)
      end if

      if(.not. allocated(eh%evec_cmplx)) then
         call elsi_allocate(eh%bh,eh%evec_cmplx,eh%bh%n_lrow,eh%bh%n_lcol,&
              "evec_cmplx",caller)
      end if

      if(.not. allocated(eh%dm_cmplx_den)) then
         call elsi_allocate(eh%bh,eh%dm_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "dm_cmplx_den",caller)
      end if

      if(.not. allocated(eh%occ)) then
         call elsi_allocate(eh%bh,eh%occ,eh%ph%n_states,eh%ph%n_spins,&
              eh%ph%n_kpts,"occ",caller)
      end if

      if(eh%ph%solver == OMM_SOLVER .and. eh%ph%omm_flavor == 0) then
         if(.not. allocated(eh%ovlp_cmplx_copy)) then
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_copy,eh%bh%n_lrow,&
                 eh%bh%n_lcol,"ovlp_cmplx_copy",caller)
         end if

         if(eh%ph%n_calls == 1) then
            eh%ovlp_cmplx_copy(:,:) = eh%ovlp_cmplx_den
         end if
      end if

      call elsi_solve_elpa(eh%ph,eh%bh,eh%ham_cmplx_den,eh%ovlp_cmplx_den,&
           eh%eval,eh%evec_cmplx)
      call elsi_get_occ_for_dm(eh%ph,eh%bh,eh%eval,eh%occ)
      call elsi_build_dm_edm(eh%ph,eh%bh,eh%occ(:,eh%ph%i_spin,eh%ph%i_kpt),&
           eh%evec_cmplx,eh%dm_cmplx_den,GET_DM)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%dm_cmplx_den,dm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%dm_cmplx_den,dm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_blacs_to_generic_dm(eh%ph,eh%bh,eh%dm_cmplx_den,eh%map_den,&
              dm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)

      eh%ph%eval_ready = .true.
      eh%ph%evec_ready = .true.
      eh%ph%occ_ready = .true.
   case(OMM_SOLVER)
      call elsi_init_blacs(eh)
      call elsi_init_elpa(eh%ph,eh%bh)
      call elsi_init_omm(eh%ph,eh%bh)

      if(.not. allocated(eh%ham_cmplx_den)) then
         call elsi_allocate(eh%bh,eh%ham_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "ham_cmplx_den",caller)
      end if

      if(.not. allocated(eh%ovlp_cmplx_den)) then
         if(.not. eh%ph%unit_ovlp) then
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,eh%bh%n_lrow,&
                 eh%bh%n_lcol,"ovlp_cmplx_den",caller)
         else
            call elsi_allocate(eh%bh,eh%ovlp_cmplx_den,1,1,"ovlp_cmplx_den",&
                 caller)
         end if
      end if

      if(.not. allocated(eh%dm_cmplx_den)) then
         call elsi_allocate(eh%bh,eh%dm_cmplx_den,eh%bh%n_lrow,eh%bh%n_lcol,&
              "dm_cmplx_den",caller)
      end if

      if(.not. allocated(eh%omm_c_cmplx)) then
         call elsi_allocate(eh%bh,eh%omm_c_cmplx,eh%ph%omm_n_lrow,eh%bh%n_lcol,&
              "omm_c_cmplx",caller)

         ! Initialize coefficient matrix with ELPA eigenvectors
         call pztranc(eh%ph%n_basis,eh%ph%n_basis,(1.0_r8,0.0_r8),&
              eh%evec_cmplx,1,1,eh%bh%desc,(0.0_r8,0.0_r8),eh%dm_cmplx_den,&
              1,1,eh%bh%desc)

         eh%omm_c_cmplx(1:eh%ph%omm_n_lrow,:)&
            = eh%dm_cmplx_den(1:eh%ph%omm_n_lrow,:)

         if(allocated(eh%evec_cmplx)) then
            call elsi_deallocate(eh%bh,eh%evec_cmplx,"evec_cmplx")
         end if

         if(allocated(eh%eval)) then
            call elsi_deallocate(eh%bh,eh%eval,"eval")
         end if

         if(allocated(eh%occ)) then
            call elsi_deallocate(eh%bh,eh%occ,"occ")
         end if
      end if

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%ham_cmplx_den,eh%ovlp_cmplx_den)
      case(SIESTA_CSC)
         call elsi_siesta_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ham_cmplx_den,eh%ovlp_cmplx_den)
      case(GENERIC_COO)
         if(.not. allocated(eh%map_den)) then
            call elsi_allocate(eh%bh,eh%map_den,eh%bh%n_lrow,eh%bh%n_lcol,&
                 "map_den",caller)
         end if

         call elsi_generic_to_blacs_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ham_cmplx_den,eh%ovlp_cmplx_den,eh%map_den)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      if(eh%ph%n_calls > 1 .and. eh%ph%n_calls == eh%ph%omm_n_elpa+1) then
         if(eh%ph%omm_flavor == 0) then
            ! Restore overlap
            eh%ovlp_cmplx_den(:,:) = eh%ovlp_cmplx_copy
         end if
      end if

      call elsi_solve_omm(eh%ph,eh%bh,eh%ham_cmplx_den,eh%ovlp_cmplx_den,&
           eh%omm_c_cmplx,eh%dm_cmplx_den)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_blacs_to_sips_dm(eh%ph,eh%bh,eh%dm_cmplx_den,dm,&
              eh%row_ind_sp1,eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_blacs_to_siesta_dm(eh%bh,eh%dm_cmplx_den,dm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_blacs_to_generic_dm(eh%ph,eh%bh,eh%dm_cmplx_den,eh%map_den,&
              dm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
   case(PEXSI_SOLVER)
      call elsi_init_pexsi(eh%ph,eh%bh)

      if(.not. allocated(eh%pexsi_ne_vec)) then
         call elsi_allocate(eh%bh,eh%pexsi_ne_vec,eh%ph%pexsi_options%nPoints,&
              "pexsi_ne_vec",caller)
      end if

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%pexsi_ne_vec,ham,ovlp,dm)
      case(SIESTA_CSC)
         if(.not. allocated(eh%row_ind_sp1)) then
            call elsi_siesta_to_pexsi_hs_dim(eh%ph,eh%bh,eh%col_ptr_sp2)

            if(eh%ph%unit_ovlp) then
               call elsi_allocate(eh%bh,eh%ovlp_cmplx_sp,1,"ovlp_cmplx_sp",&
                    caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_cmplx_sp,eh%bh%nnz_l_sp1,&
                    "ovlp_cmplx_sp",caller)
            end if

            call elsi_allocate(eh%bh,eh%ham_cmplx_sp,eh%bh%nnz_l_sp1,&
                 "ham_cmplx_sp",caller)
            call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                 "row_ind_sp1",caller)
            call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                 "col_ptr_sp1",caller)
         end if

         call elsi_siesta_to_pexsi_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%ham_cmplx_sp,eh%ovlp_cmplx_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1)

         if(.not. allocated(eh%dm_cmplx_sp)) then
            call elsi_allocate(eh%bh,eh%dm_cmplx_sp,eh%bh%nnz_l_sp1,&
                 "dm_cmplx_sp",caller)
         end if

         eh%dm_cmplx_sp(:) = (0.0_r8,0.0_r8)

         call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%pexsi_ne_vec,eh%ham_cmplx_sp,eh%ovlp_cmplx_sp,eh%dm_cmplx_sp)
         call elsi_pexsi_to_siesta_dm(eh%ph,eh%bh,eh%dm_cmplx_sp,&
              eh%row_ind_sp1,eh%col_ptr_sp1,dm,eh%row_ind_sp2,eh%col_ptr_sp2)
      case(GENERIC_COO)
         if(.not. allocated(eh%row_ind_sp1)) then
            call elsi_generic_to_pexsi_hs_dim(eh%ph,eh%bh,eh%col_ind_sp3)

            if(eh%ph%unit_ovlp) then
               call elsi_allocate(eh%bh,eh%ovlp_cmplx_sp,1,"ovlp_cmplx_sp",&
                    caller)
            else
               call elsi_allocate(eh%bh,eh%ovlp_cmplx_sp,eh%bh%nnz_l_sp1,&
                    "ovlp_cmplx_sp",caller)
            end if

            call elsi_allocate(eh%bh,eh%ham_cmplx_sp,eh%bh%nnz_l_sp1,&
                 "ham_cmplx_sp",caller)
            call elsi_allocate(eh%bh,eh%row_ind_sp1,eh%bh%nnz_l_sp1,&
                 "row_ind_sp1",caller)
            call elsi_allocate(eh%bh,eh%col_ptr_sp1,eh%bh%n_lcol_sp1+1,&
                 "col_ptr_sp1",caller)
            call elsi_allocate(eh%bh,eh%map_sp1,eh%bh%nnz_l_sp1,"map_sp1",&
                 caller)
         end if

         call elsi_generic_to_pexsi_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%ham_cmplx_sp,eh%ovlp_cmplx_sp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%map_sp1)

         if(.not. allocated(eh%dm_cmplx_sp)) then
            call elsi_allocate(eh%bh,eh%dm_cmplx_sp,eh%bh%nnz_l_sp1,&
                 "dm_cmplx_sp",caller)
         end if

         call elsi_solve_pexsi(eh%ph,eh%bh,eh%row_ind_sp1,eh%col_ptr_sp1,&
              eh%pexsi_ne_vec,eh%ham_cmplx_sp,eh%ovlp_cmplx_sp,eh%dm_cmplx_sp)
         call elsi_pexsi_to_generic_dm(eh%ph,eh%bh,eh%dm_cmplx_sp,&
              eh%row_ind_sp1,eh%col_ptr_sp1,eh%map_sp1,dm,eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
   case(NTPOLY_SOLVER)
      call elsi_init_ntpoly(eh%ph,eh%bh)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_sips_to_ntpoly_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%nt_ham,eh%nt_ovlp)
      case(SIESTA_CSC)
         call elsi_siesta_to_ntpoly_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%nt_ham,eh%nt_ovlp)
      case(GENERIC_COO)
         call elsi_generic_to_ntpoly_hs(eh%ph,eh%bh,ham,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%nt_ham,eh%nt_ovlp,eh%nt_map)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_solve_ntpoly(eh%ph,eh%bh,eh%nt_ham,eh%nt_ovlp,eh%nt_dm)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         call elsi_ntpoly_to_sips_dm(eh%ph,eh%bh,eh%nt_dm,dm,eh%row_ind_sp1,&
              eh%col_ptr_sp1)
      case(SIESTA_CSC)
         call elsi_ntpoly_to_siesta_dm(eh%bh,eh%nt_dm,dm,eh%row_ind_sp2,&
              eh%col_ptr_sp2)
      case(GENERIC_COO)
         call elsi_ntpoly_to_generic_dm(eh%ph,eh%bh,eh%nt_dm,eh%nt_map,dm,&
              eh%perm_sp3)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select

      call elsi_get_band_energy(eh%ph,eh%bh,ebs,solver)
   case default
      write(msg,"(A)") "Unsupported density matrix solver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   eh%ph%edm_ready = .true.

   if(eh%ph%save_ovlp) then
      call elsi_init_ntpoly(eh%ph,eh%bh)

      select case(eh%ph%matrix_format)
      case(PEXSI_CSC)
         if(eh%ph%n_calls == 1) then
            eh%ph%first_sips_to_ntpoly = .true.
         end if

         call elsi_sips_to_ntpoly_hs(eh%ph,eh%bh,dm,ovlp,eh%row_ind_sp1,&
              eh%col_ptr_sp1,eh%nt_dm_copy,eh%nt_ovlp_copy)
      case(SIESTA_CSC)
         if(eh%ph%n_calls == 1) then
            eh%ph%first_siesta_to_ntpoly = .true.
         end if

         call elsi_siesta_to_ntpoly_hs(eh%ph,eh%bh,dm,ovlp,eh%row_ind_sp2,&
              eh%col_ptr_sp2,eh%nt_dm_copy,eh%nt_ovlp_copy)
      case(GENERIC_COO)
         if(eh%ph%n_calls == 1) then
            eh%ph%first_generic_to_ntpoly = .true.
         end if

         call elsi_generic_to_ntpoly_hs(eh%ph,eh%bh,dm,ovlp,eh%row_ind_sp3,&
              eh%col_ind_sp3,eh%nt_dm_copy,eh%nt_ovlp_copy,eh%nt_map)
      case default
         write(msg,"(A)") "Unsupported matrix format"
         call elsi_stop(eh%bh,msg,caller)
      end select
   end if

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! Compute the eigenvalues and eigenvectors of a Bethe-Salpeter problem.
!!
subroutine elsi_bse_real(eh,mat_a,mat_b,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   real(kind=r8), intent(inout) :: mat_a(eh%bh%n_lrow,eh%bh%n_lcol) !< Matrix A
   real(kind=r8), intent(in) :: mat_b(eh%bh%n_lrow,eh%bh%n_lcol) !< Matrix B
   real(kind=r8), intent(out) :: eval(eh%ph%n_basis) !< Eigenvalues
   real(kind=r8), intent(out) :: evec(eh%ph%bse_n_lrow,eh%ph%bse_n_lcol) !< Eigenvectors

   real(kind=r8) :: t0
   character(len=29) :: dt0
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_bse_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   select case(eh%ph%solver)
   case(BSEPACK_SOLVER)
      call elsi_solve_bsepack(eh%ph,eh%bh,mat_a,mat_b,eval,evec)
   case default
      write(msg,"(A)") "Unsupported BSE eigensolver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! Compute the eigenvalues and eigenvectors of a Bethe-Salpeter problem.
!!
subroutine elsi_bse_complex(eh,mat_a,mat_b,eval,evec)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   complex(kind=r8), intent(inout) :: mat_a(eh%bh%n_lrow,eh%bh%n_lcol) !< Matrix A
   complex(kind=r8), intent(in) :: mat_b(eh%bh%n_lrow,eh%bh%n_lcol) !< Matrix B
   real(kind=r8), intent(out) :: eval(eh%ph%n_basis) !< Eigenvalues
   complex(kind=r8), intent(out) :: evec(eh%ph%bse_n_lrow,eh%ph%bse_n_lcol) !< Eigenvectors

   real(kind=r8) :: t0
   character(len=29) :: dt0
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_bse_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)
   call elsi_check(eh%ph,eh%bh,caller)
   call elsi_get_time(t0)
   call fjson_get_datetime_rfc3339(dt0)

   select case(eh%ph%solver)
   case(BSEPACK_SOLVER)
      call elsi_solve_bsepack(eh%ph,eh%bh,mat_a,mat_b,eval,evec)
   case default
      write(msg,"(A)") "Unsupported BSE eigensolver"
      call elsi_stop(eh%bh,msg,caller)
   end select

   call elsi_add_log(eh%ph,eh%bh,eh%jh,dt0,t0,caller)

end subroutine

!>
!! Construct density matrix from eigenvectors.
!!
subroutine elsi_compute_dm_real(eh,occ,evec,dm)

   implicit none

   type(elsi_handle), intent(in) :: eh
   real(kind=r8), intent(in) :: occ(eh%ph%n_states)
   real(kind=r8), intent(in) :: evec(eh%bh%n_lrow,eh%bh%n_lcol)
   real(kind=r8), intent(out) :: dm(eh%bh%n_lrow,eh%bh%n_lcol)

   character(len=*), parameter :: caller = "elsi_compute_dm_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_build_dm_edm(eh%ph,eh%bh,occ,evec,dm,GET_DM)

end subroutine

!>
!! Construct density matrix from eigenvectors.
!!
subroutine elsi_compute_dm_complex(eh,occ,evec,dm)

   implicit none

   type(elsi_handle), intent(in) :: eh
   real(kind=r8), intent(in) :: occ(eh%ph%n_states)
   complex(kind=r8), intent(in) :: evec(eh%bh%n_lrow,eh%bh%n_lcol)
   complex(kind=r8), intent(out) :: dm(eh%bh%n_lrow,eh%bh%n_lcol)

   character(len=*), parameter :: caller = "elsi_compute_dm_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_build_dm_edm(eh%ph,eh%bh,occ,evec,dm,GET_DM)

end subroutine

!>
!! Construct energy-weighted density matrix from eigenvectors.
!!
subroutine elsi_compute_edm_real(eh,eval,occ,evec,edm)

   implicit none

   type(elsi_handle), intent(in) :: eh
   real(kind=r8), intent(in) :: eval(eh%ph%n_states)
   real(kind=r8), intent(in) :: occ(eh%ph%n_states)
   real(kind=r8), intent(in) :: evec(eh%bh%n_lrow,eh%bh%n_lcol)
   real(kind=r8), intent(out) :: edm(eh%bh%n_lrow,eh%bh%n_lcol)

   real(kind=r8), allocatable :: factor(:)

   character(len=*), parameter :: caller = "elsi_compute_edm_real"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_allocate(eh%bh,factor,eh%ph%n_states,"factor",caller)

   factor(1:eh%ph%n_states_solve) = -occ(1:eh%ph%n_states_solve)&
      *eval(1:eh%ph%n_states_solve)

   call elsi_build_dm_edm(eh%ph,eh%bh,factor,evec,edm,GET_EDM)

   call elsi_deallocate(eh%bh,factor,"factor")

end subroutine

!>
!! Construct energy-weighted density matrix from eigenvectors.
!!
subroutine elsi_compute_edm_complex(eh,eval,occ,evec,edm)

   implicit none

   type(elsi_handle), intent(in) :: eh
   real(kind=r8), intent(in) :: eval(eh%ph%n_states)
   real(kind=r8), intent(in) :: occ(eh%ph%n_states)
   complex(kind=r8), intent(in) :: evec(eh%bh%n_lrow,eh%bh%n_lcol)
   complex(kind=r8), intent(out) :: edm(eh%bh%n_lrow,eh%bh%n_lcol)

   real(kind=r8), allocatable :: factor(:)

   character(len=*), parameter :: caller = "elsi_compute_edm_complex"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   call elsi_allocate(eh%bh,factor,eh%ph%n_states,"factor",caller)

   factor(1:eh%ph%n_states_solve) = -occ(1:eh%ph%n_states_solve)&
      *eval(1:eh%ph%n_states_solve)

   call elsi_build_dm_edm(eh%ph,eh%bh,factor,evec,edm,GET_EDM)

   call elsi_deallocate(eh%bh,factor,"factor")

end subroutine

!>
!! Compute the chemical potential and occupation numbers.
!!
subroutine elsi_compute_mu_and_occ(eh,n_electron,n_state,n_spin,n_kpt,k_wt,&
   eval,occ,mu)

   implicit none

   type(elsi_handle), intent(in) :: eh !< Handle
   real(kind=r8), intent(in) :: n_electron !< Number of electrons
   integer(kind=i4), intent(in) :: n_state !< Number of states
   integer(kind=i4), intent(in) :: n_spin !< Number of spins
   integer(kind=i4), intent(in) :: n_kpt !< Number of k-points
   real(kind=r8), intent(in) :: k_wt(n_kpt) !< K-points weights
   real(kind=r8), intent(in) :: eval(n_state,n_spin,n_kpt) !< Eigenvalues
   real(kind=r8), intent(out) :: occ(n_state,n_spin,n_kpt) !< Occupation members
   real(kind=r8), intent(out) :: mu !< Chemical potential

   character(len=*), parameter :: caller = "elsi_compute_mu_and_occ"

   call elsi_mu_and_occ(eh%ph,eh%bh,n_electron,n_state,n_spin,n_kpt,k_wt,eval,&
        occ,mu)

end subroutine

!>
!! Compute the electronic entropy.
!!
subroutine elsi_compute_entropy(eh,n_state,n_spin,n_kpt,k_wt,eval,occ,mu,ts)

   implicit none

   type(elsi_handle), intent(in) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_state !< Number of states
   integer(kind=i4), intent(in) :: n_spin !< Number of spins
   integer(kind=i4), intent(in) :: n_kpt !< Number of k-points
   real(kind=r8), intent(in) :: k_wt(n_kpt) !< K-points weights
   real(kind=r8), intent(in) :: eval(n_state,n_spin,n_kpt) !< Eigenvalues
   real(kind=r8), intent(in) :: occ(n_state,n_spin,n_kpt) !< Occupation numbers
   real(kind=r8), intent(in) :: mu !< Input chemical potential
   real(kind=r8), intent(out) :: ts !< Entropy

   character(len=*), parameter :: caller = "elsi_compute_entropy"

   call elsi_entropy(eh%ph,n_state,n_spin,n_kpt,k_wt,eval,occ,mu,ts)

end subroutine

end module ELSI_SOLVER
