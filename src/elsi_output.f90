! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Manage output to stdout and/or a JSON log file.
!!
module ELSI_OUTPUT

   use ELSI_CONSTANTS, only: MULTI_PROC,SINGLE_PROC,BLACS_DENSE,PEXSI_CSC,&
       SIESTA_CSC,GENERIC_COO,ELPA_SOLVER,PEXSI_SOLVER,SIPS_SOLVER,OMM_SOLVER,&
       NTPOLY_SOLVER
   use ELSI_DATATYPE, only: elsi_param_t,elsi_basic_t
   use ELSI_PRECISION, only: r8,i4,i8
   use FORTJSON, only: fjson_write_name_value,fjson_reset_fj_handle,&
       fjson_start_name_object,fjson_start_array,fjson_finish_object,&
       fjson_open_file,fjson_handle,fjson_start_object,&
       fjson_get_datetime_rfc3339,fjson_close_file,fjson_finish_array

   implicit none

   private

   public :: elsi_say
   public :: elsi_add_log
   public :: elsi_get_time
   public :: elsi_final_print
   public :: fjson_close_file
   public :: fjson_finish_array
   public :: fjson_get_datetime_rfc3339
   public :: fjson_reset_fj_handle

contains

!>
!! Print a message.
!!
subroutine elsi_say(bh,msg)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   character(len=*), intent(in) :: msg

   character(len=*), parameter :: caller = "elsi_say"

   if(bh%print_info > 0) then
      write(bh%print_unit,"(2X,A)") trim(msg)
   end if

end subroutine

!>
!! Add an entry to the log file.
!!
subroutine elsi_add_log(ph,bh,jh,dt0,t0,caller)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(inout) :: bh
   type(fjson_handle), intent(inout) :: jh
   character(len=*), intent(in) :: dt0
   real(kind=r8), intent(in) :: t0
   character(len=*), intent(in) :: caller

   integer(kind=i4) :: solver_use
   real(kind=r8) :: t1
   character(len=20) :: solver_tag
   character(len=29) :: dt_record

   if(bh%print_json > 0) then
      if(.not. bh%json_init) then
         call fjson_open_file(jh,66,"elsi_log.json")
         call fjson_start_array(jh)

         if(.not. bh%uuid_ready) then
            call elsi_gen_uuid(bh%uuid)
            bh%uuid_ready = .true.
         end if

         bh%json_init = .true.
      end if

      call elsi_get_time(t1)
      call fjson_get_datetime_rfc3339(dt_record)

      solver_use = 1

      select case(ph%solver)
      case(ELPA_SOLVER)
         if(ph%parallel_mode == SINGLE_PROC) then
            solver_tag = "LAPACK"
         else
            solver_tag = "ELPA"
         end if
      case(OMM_SOLVER)
         if(ph%n_calls <= ph%omm_n_elpa) then
            solver_tag = "ELPA"
         else
            solver_tag = "LIBOMM"
            solver_use = 2
         end if
      case(PEXSI_SOLVER)
         solver_tag = "PEXSI"
         solver_use = 3
      case(SIPS_SOLVER)
         if(ph%n_calls <= ph%sips_n_elpa) then
            solver_tag = "ELPA"
         else
            solver_tag = "SLEPC_SIPS"
            solver_use = 5
         end if
      case(NTPOLY_SOLVER)
         solver_tag = "NTPOLY"
         solver_use = 6
      end select

      call fjson_start_object(jh)
      call elsi_print_version_summary(bh%uuid,jh)
      call fjson_write_name_value(jh,"step",ph%n_calls)
      call fjson_write_name_value(jh,"total_step",ph%n_calls+ph%n_calls_all)

      if(caller(6:6) == "e") then
         call fjson_write_name_value(jh,"output_type","EIGENSOLUTION")
      else
         call fjson_write_name_value(jh,"output_type","DENSITY MATRIX")
      end if

      if(caller(9:9) == "r") then
         call fjson_write_name_value(jh,"data_type","REAL")
      else
         call fjson_write_name_value(jh,"data_type","COMPLEX")
      end if

      call fjson_write_name_value(jh,"elsi_tag",trim(solver_tag))
      call fjson_write_name_value(jh,"user_tag",trim(bh%user_tag))
      call fjson_write_name_value(jh,"start_datetime",dt0)
      call fjson_write_name_value(jh,"record_datetime",dt_record)
      call fjson_write_name_value(jh,"total_time",t1-t0)
      call elsi_print_handle_summary(ph,bh,jh)
      call elsi_print_ovlp_summary(ph,jh)
      call fjson_write_name_value(jh,"solver_used",trim(solver_tag))

      select case(solver_use)
      case(NTPOLY_SOLVER)
         call elsi_print_ntpoly_settings(ph,jh)
      case(ELPA_SOLVER)
         call elsi_print_elpa_settings(ph,jh)
      case(OMM_SOLVER)
         call elsi_print_omm_settings(ph,jh)
      case(PEXSI_SOLVER)
         call elsi_print_pexsi_settings(ph,jh)
      case(SIPS_SOLVER)
         call elsi_print_sips_settings(ph,jh)
      end select

      if(ph%matrix_format == BLACS_DENSE) then
         call elsi_print_dense_settings(bh,jh)
      else
         call elsi_print_sparse_settings(bh,jh)
      end if

      call fjson_finish_object(jh)
   end if

end subroutine

!>
!! Print the state of the handle.
!!
subroutine elsi_print_handle_summary(ph,bh,jh)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   type(fjson_handle), intent(inout) :: jh

   real(kind=r8) :: sparsity

   character(len=*), parameter :: caller = "elsi_print_handle_summary"

   call fjson_write_name_value(jh,"n_electrons",ph%n_electrons)

   if(ph%parallel_mode == MULTI_PROC) then
      call fjson_write_name_value(jh,"n_spin",ph%n_spins)
      call fjson_write_name_value(jh,"n_kpts",ph%n_kpts)
   end if

   if(ph%solver == ELPA_SOLVER .or. ph%solver == SIPS_SOLVER) then
      call fjson_write_name_value(jh,"n_states",ph%n_states)
   end if

   select case(ph%matrix_format)
   case(BLACS_DENSE)
      call fjson_write_name_value(jh,"matrix_format","BLACS_DENSE")
   case(PEXSI_CSC)
      call fjson_write_name_value(jh,"matrix_format","PEXSI_CSC")
   case(SIESTA_CSC)
      call fjson_write_name_value(jh,"matrix_format","SIESTA_CSC")
   case(GENERIC_COO)
      call fjson_write_name_value(jh,"matrix_format","GENERIC_COO")
   end select

   call fjson_write_name_value(jh,"n_basis",ph%n_basis)

   if(ph%parallel_mode == MULTI_PROC) then
      sparsity = 1.0_r8-(1.0_r8*bh%nnz_g/ph%n_basis/ph%n_basis)
      call fjson_write_name_value(jh,"sparsity",sparsity)
      call fjson_write_name_value(jh,"nnz_g",bh%nnz_g)
   end if

   if(ph%parallel_mode == MULTI_PROC) then
      call fjson_write_name_value(jh,"parallel_mode","MULTI_PROC")
      call fjson_write_name_value(jh,"n_procs",bh%n_procs)
      call fjson_write_name_value(jh,"n_procs_all",bh%n_procs_all)
   else if(ph%parallel_mode == SINGLE_PROC) then
      call fjson_write_name_value(jh,"parallel_mode","SINGLE_PROC")
   end if

   select case(ph%solver)
   case(ELPA_SOLVER)
      call fjson_write_name_value(jh,"solver_chosen","ELPA")
   case(NTPOLY_SOLVER)
      call fjson_write_name_value(jh,"solver_chosen","NTPOLY")
   case(OMM_SOLVER)
      call fjson_write_name_value(jh,"solver_chosen","libOMM")
   case(PEXSI_SOLVER)
      call fjson_write_name_value(jh,"solver_chosen","PEXSI")
   case(SIPS_SOLVER)
      call fjson_write_name_value(jh,"solver_chosen","SLEPc_SIPs")
   end select

end subroutine

!>
!! Print versioning information.
!!
subroutine elsi_print_version_summary(uuid,jh)

   implicit none

   character(len=*), intent(in) :: uuid
   type(fjson_handle), intent(inout) :: jh

   character(len=8) :: VERSION
   character(len=8) :: DATESTAMP
   character(len=8) :: COMMIT
   character(len=40) :: HOSTNAME
   character(len=20) :: DATETIME

   character(len=*), parameter :: caller = "elsi_print_version_summary"

   call elsi_version_info(VERSION,DATESTAMP,COMMIT,HOSTNAME,DATETIME)

   call fjson_write_name_value(jh,"data_source","ELSI")
   call fjson_write_name_value(jh,"code_version",trim(VERSION))
   call fjson_write_name_value(jh,"code_date_stamp",trim(DATESTAMP))
   call fjson_write_name_value(jh,"git_commit_abbrev",trim(COMMIT))
   call fjson_write_name_value(jh,"compiled_on_hostname",trim(HOSTNAME))
   call fjson_write_name_value(jh,"compiled_at_datetime",trim(DATETIME))
   call fjson_write_name_value(jh,"uuid",trim(uuid))

end subroutine

!>
!! Print information about the overlap matrix.
!!
subroutine elsi_print_ovlp_summary(ph,jh)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(fjson_handle), intent(inout) :: jh

   character(len=*), parameter :: caller = "elsi_print_ovlp_summary"

   call fjson_write_name_value(jh,"save_ovlp",ph%save_ovlp)
   call fjson_write_name_value(jh,"unit_ovlp",ph%unit_ovlp)
   call fjson_write_name_value(jh,"ill_check",ph%ill_check)
   call fjson_write_name_value(jh,"ill_ovlp",ph%ill_ovlp)
   call fjson_write_name_value(jh,"ill_tol",ph%ill_tol)
   call fjson_write_name_value(jh,"n_illcond",ph%n_basis-ph%n_good)
   call fjson_write_name_value(jh,"n_states_solve",ph%n_states_solve)
   call fjson_write_name_value(jh,"ovlp_ev_min",ph%ovlp_ev_min)
   call fjson_write_name_value(jh,"ovlp_ev_max",ph%ovlp_ev_max)

end subroutine

!>
!! Print settings for ELPA.
!!
subroutine elsi_print_elpa_settings(ph,jh)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(fjson_handle), intent(inout) :: jh

   character(len=*), parameter :: caller = "elsi_print_elpa_settings"

   call fjson_start_name_object(jh,"solver_settings")
   call fjson_write_name_value(jh,"elpa_solver",ph%elpa_solver)
   call fjson_write_name_value(jh,"elpa_n_states",ph%n_states)
   call fjson_write_name_value(jh,"elpa_n_single",ph%elpa_n_single)
   call fjson_write_name_value(jh,"elpa_gpu",ph%elpa_gpu)
   call fjson_write_name_value(jh,"elpa_gpu_kernels",ph%elpa_gpu_kernels)
   call fjson_write_name_value(jh,"elpa_autotune",ph%elpa_autotune)
   call fjson_finish_object(jh)

end subroutine

!>
!! Print settings for NTPoly.
!!
subroutine elsi_print_ntpoly_settings(ph,jh)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(fjson_handle), intent(inout) :: jh

   character(len=*), parameter :: caller = "elsi_print_ntpoly_settings"

   call fjson_start_name_object(jh,"solver_settings")
   call fjson_write_name_value(jh,"nt_n_layers",ph%nt_n_layers)
   call fjson_write_name_value(jh,"nt_n_prow",ph%nt_n_prow)
   call fjson_write_name_value(jh,"nt_n_pcol",ph%nt_n_pcol)
   call fjson_write_name_value(jh,"nt_method",ph%nt_method)
   call fjson_write_name_value(jh,"nt_tol",ph%nt_tol)
   call fjson_write_name_value(jh,"nt_filter",ph%nt_filter)
   call fjson_write_name_value(jh,"nt_max_iter",ph%nt_max_iter)
   call fjson_finish_object(jh)

end subroutine

!>
!! Print settings for libOMM.
!!
subroutine elsi_print_omm_settings(ph,jh)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(fjson_handle), intent(inout) :: jh

   character(len=*), parameter :: caller = "elsi_print_omm_settings"

   call fjson_start_name_object(jh,"solver_settings")
   call fjson_write_name_value(jh,"omm_n_states",ph%omm_n_states)
   call fjson_write_name_value(jh,"omm_n_elpa",ph%omm_n_elpa)
   call fjson_write_name_value(jh,"omm_flavor",ph%omm_flavor)
   call fjson_write_name_value(jh,"omm_tol",ph%omm_tol)
   call fjson_finish_object(jh)

end subroutine

!>
!! Print settings for PEXSI.
!!
subroutine elsi_print_pexsi_settings(ph,jh)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(fjson_handle), intent(inout) :: jh

   character(len=*), parameter :: caller = "elsi_print_pexsi_settings"

   call fjson_start_name_object(jh,"solver_settings")
   call fjson_write_name_value(jh,"pexsi_n_pole",ph%pexsi_options%numPole)
   call fjson_write_name_value(jh,"pexsi_n_point",ph%pexsi_options%nPoints)
   call fjson_write_name_value(jh,"pexsi_np_per_pole",ph%pexsi_np_per_pole)
   call fjson_write_name_value(jh,"pexsi_np_per_point",ph%pexsi_np_per_point)
   call fjson_write_name_value(jh,"pexsi_n_prow",ph%pexsi_n_prow)
   call fjson_write_name_value(jh,"pexsi_n_pcol",ph%pexsi_n_pcol)
   call fjson_write_name_value(jh,"pexsi_temperature",&
        ph%pexsi_options%temperature)
   call fjson_write_name_value(jh,"pexsi_delta_e",ph%pexsi_options%deltaE)
   call fjson_write_name_value(jh,"pexsi_gap",ph%pexsi_options%gap)
   call fjson_write_name_value(jh,"pexsi_mu_min",ph%pexsi_options%muMin0)
   call fjson_write_name_value(jh,"pexsi_mu_max",ph%pexsi_options%muMax0)
   call fjson_write_name_value(jh,"pexsi_do_inertia",&
        ph%pexsi_options%isInertiaCount)
   call fjson_write_name_value(jh,"pexsi_intertia_tol",&
        ph%pexsi_options%muInertiaTolerance)
   call fjson_write_name_value(jh,"pexsi_tol",&
        ph%pexsi_options%numElectronPEXSITolerance)
   call fjson_write_name_value(jh,"pexsi_np_symbfact",&
        ph%pexsi_options%npSymbFact)
   call fjson_write_name_value(jh,"pexsi_reordering",ph%pexsi_options%ordering)
   call fjson_finish_object(jh)

end subroutine

!>
!! Print settings for SLEPc-SIPs.
!!
subroutine elsi_print_sips_settings(ph,jh)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(fjson_handle), intent(inout) :: jh

   character(len=*), parameter :: caller = "elsi_print_sips_settings"

   call fjson_start_name_object(jh,"solver_settings")
   call fjson_write_name_value(jh,"sips_n_states",ph%n_states)
   call fjson_write_name_value(jh,"sips_n_elpa",ph%sips_n_elpa)
   call fjson_write_name_value(jh,"sips_slice_type",ph%sips_slice_type)
   call fjson_write_name_value(jh,"sips_n_slices",ph%sips_n_slices)
   call fjson_write_name_value(jh,"sips_buffer",ph%sips_buffer)
   call fjson_finish_object(jh)

end subroutine

!>
!! Print settings for the dense matrix format.
!!
subroutine elsi_print_dense_settings(bh,jh)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   type(fjson_handle), intent(inout) :: jh

   character(len=*), parameter :: caller = "elsi_print_dense_settings"

   call fjson_start_name_object(jh,"matrix_format_settings")
   call fjson_write_name_value(jh,"blk",bh%blk)
   call fjson_write_name_value(jh,"n_prow",bh%n_prow)
   call fjson_write_name_value(jh,"n_pcol",bh%n_pcol)
   call fjson_write_name_value(jh,"blacs_ready",bh%blacs_ready)
   call fjson_finish_object(jh)

end subroutine

!>
!! Print settings for the sparse matrix formats.
!!
subroutine elsi_print_sparse_settings(bh,jh)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   type(fjson_handle), intent(inout) :: jh

   character(len=*), parameter :: caller = "elsi_print_sparse_settings"

   call fjson_start_name_object(jh,"matrix_format_settings")
   call fjson_write_name_value(jh,"def0",bh%def0)
   call fjson_write_name_value(jh,"blk_sp2",bh%blk_sp2)
   call fjson_write_name_value(jh,"pexsi_csc_ready",bh%pexsi_csc_ready)
   call fjson_write_name_value(jh,"siesta_csc_ready",bh%siesta_csc_ready)
   call fjson_write_name_value(jh,"generic_coo_ready",bh%generic_coo_ready)
   call fjson_finish_object(jh)

end subroutine

!>
!! Print a final summary.
!!
subroutine elsi_final_print(ph,bh)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh

   real(kind=r8) :: sparsity
   character(len=200) :: ll
   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_final_print"

   write(ll,"(A)") "|------------------------------------------------------"
   call elsi_say(bh,ll)

   write(msg,"(A)") "| Final ELSI Output"
   call elsi_say(bh,msg)

   call elsi_say(bh,ll)

   write(msg,"(A)") "|"
   call elsi_say(bh,msg)

   write(msg,"(A)") "| Physical Properties"
   call elsi_say(bh,msg)

   write(msg,"(A,E22.8)") "|   Number of electrons       :",ph%n_electrons
   call elsi_say(bh,msg)

   if(ph%parallel_mode == MULTI_PROC) then
      write(msg,"(A,I22)") "|   Number of spins           :",ph%n_spins
      call elsi_say(bh,msg)

      write(msg,"(A,I22)") "|   Number of k-points        :",ph%n_kpts
      call elsi_say(bh,msg)
   end if

   if(ph%solver == ELPA_SOLVER .or. ph%solver == SIPS_SOLVER) then
      write(msg,"(A,I22)") "|   Number of states          :",ph%n_states
      call elsi_say(bh,msg)
   end if

   write(msg,"(A)") "|"
   call elsi_say(bh,msg)

   write(msg,"(A)") "| Matrix Properties"
   call elsi_say(bh,msg)

   select case(ph%matrix_format)
   case(BLACS_DENSE)
      write(msg,"(A,A22)") "|   Matrix format             :","BLACS_DENSE"
      call elsi_say(bh,msg)
   case(PEXSI_CSC)
      write(msg,"(A,A22)") "|   Matrix format             :","PEXSI_CSC"
      call elsi_say(bh,msg)
   case(SIESTA_CSC)
      write(msg,"(A,A22)") "|   Matrix format             :","SIESTA_CSC"
      call elsi_say(bh,msg)
   case(GENERIC_COO)
      write(msg,"(A,A22)") "|   Matrix format             :","GENERIC_COO"
      call elsi_say(bh,msg)
   end select

   write(msg,"(A,I22)") "|   Number of basis functions :",ph%n_basis
   call elsi_say(bh,msg)

   if(ph%parallel_mode == MULTI_PROC) then
      sparsity = 1.0_r8-(1.0_r8*bh%nnz_g/ph%n_basis/ph%n_basis)

      write(msg,"(A,E22.8)") "|   Matrix sparsity           :",sparsity
      call elsi_say(bh,msg)
   end if

   write(msg,"(A)") "|"
   call elsi_say(bh,msg)

   write(msg,"(A)") "| Computational Details"
   call elsi_say(bh,msg)

   if(ph%parallel_mode == MULTI_PROC) then
      write(msg,"(A,A22)") "|   Parallel mode             :","MULTI_PROC"
      call elsi_say(bh,msg)

      write(msg,"(A,I22)") "|   Number of MPI tasks       :",bh%n_procs_all
      call elsi_say(bh,msg)
   else if(ph%parallel_mode == SINGLE_PROC) then
      write(msg,"(A,A22)") "|   Parallel mode             :","SINGLE_PROC"
      call elsi_say(bh,msg)
   end if

   select case(ph%solver)
   case(ELPA_SOLVER)
      write(msg,"(A,A22)") "|   Solver requested          :","ELPA"
      call elsi_say(bh,msg)
   case(OMM_SOLVER)
      write(msg,"(A,A22)") "|   Solver requested          :","libOMM"
      call elsi_say(bh,msg)
   case(PEXSI_SOLVER)
      write(msg,"(A,A22)") "|   Solver requested          :","PEXSI"
      call elsi_say(bh,msg)
   case(SIPS_SOLVER)
      write(msg,"(A,A22)") "|   Solver requested          :","SLEPc-SIPs"
      call elsi_say(bh,msg)
   case(NTPOLY_SOLVER)
      write(msg,"(A,A22)") "|   Solver requested          :","NTPoly"
      call elsi_say(bh,msg)
   end select

   write(msg,"(A,I22)") "|   Number of ELSI calls      :",ph%n_calls_all
   call elsi_say(bh,msg)

   write(msg,"(A)") "|"
   call elsi_say(bh,msg)

   call elsi_say(bh,ll)

   write(msg,"(A)") "| ELSI Project (c)  elsi-interchange.org"
   call elsi_say(bh,msg)

   call elsi_say(bh,ll)

end subroutine

!>
!! Get the current wallclock time.
!! (Taken from FHI-aims with permission of copyright holders)
!!
subroutine elsi_get_time(wtime)

   implicit none

   real(kind=r8), intent(out) :: wtime

   character(len=8) :: cdate
   character(len=10) :: ctime
   integer(kind=i4) :: val
   integer(kind=i4) :: int_year
   real(kind=r8) :: year
   real(kind=r8) :: day
   real(kind=r8) :: hour
   real(kind=r8) :: minute
   real(kind=r8) :: second
   real(kind=r8) :: millisecond

   character(len=*), parameter :: caller = "elsi_get_time"

   call date_and_time(cdate,ctime)

   read(cdate(1:4),"(I4)") val

   int_year = val
   year = real(val,kind=r8)-2009.0_r8 ! 2009 is an arbitrary zero
   day = year*365+floor(year/4.0_r8)

   read(cdate(5:6),"(I2)") val

   val = val-1

   do while(val > 0)
      select case(val)
      case(1)
         day = day+31
      case(2)
         if(mod(int_year,4) == 0) then
            day = day+29
         else
            day = day+28
         end if
      case(3)
         day = day+31
      case(4)
         day = day+30
      case(5)
         day = day+31
      case(6)
         day = day+30
      case(7)
         day = day+31
      case(8)
         day = day+31
      case(9)
         day = day+30
      case(10)
         day = day+31
      case(11)
         day = day+30
      end select

      val = val-1
   end do

   read(cdate(7:8),"(I2)") val
   day = day+real(val,kind=r8)-1

   read(ctime(1:2),"(I2)") val
   hour = real(val,kind=r8)

   read(ctime(3:4),"(I2)") val
   minute = real(val,kind=r8)

   read(ctime(5:6),"(I2)") val
   second = real(val,kind=r8)

   read(ctime(8:10),"(I3)") val
   millisecond = real(val,kind=r8)

   wtime = day*86400.0_r8+hour*3600.0_r8+minute*60.0_r8+second&
      +millisecond*0.001_r8

end subroutine

!>
!! Generate a UUID (unique identifier) in RFC 4122 format.
!! (Taken from FHI-aims with permission of copyright holders)
!!
subroutine elsi_gen_uuid(uuid)

   implicit none

   character(len=36), intent(out) :: uuid

   integer(kind=i4) :: ii3
   integer(kind=i4) :: ii4
   integer(kind=i4) :: i_entry
   real(kind=r8) :: rr(8)
   character(len=3) :: ss3
   character(len=4) :: ss4

   character(len=*), parameter :: caller = "elsi_gen_uuid"

   call elsi_init_random_seed()

   ii3 = 4095
   ii4 = 65535

   call random_number(rr)

   do i_entry=1,8
      write(ss3,"(Z3.3)") transfer(int(rr(i_entry)*ii3),16)
      write(ss4,"(Z4.4)") transfer(int(rr(i_entry)*ii4),16)

      select case(i_entry)
      case(1)
         write(uuid,"(A)") ss4
      case(2)
         write(uuid,"(2A)") trim(uuid),ss4
      case(3)
         write(uuid,"(3A)") trim(uuid),"-",ss4
      case(4)
         write(uuid,"(3A)") trim(uuid),"-4",ss3
      case(5)
         write(uuid,"(4A)") trim(uuid),"-A",ss3,"-"
      case default
         write(uuid,"(2A)") trim(uuid),ss4
      end select
   end do

end subroutine

!>
!! Linear congruential generator.
!! (Taken from FHI-aims with permission of copyright holders)
!!
integer(kind=i4) function lcg(s)

   implicit none

   integer(kind=i8) :: s

   if(s == 0) then
      s = 104729
   else
      s = mod(s,int(huge(0_2)*2,kind=i8))
   end if

   s = mod(s*int(huge(0_2),kind=i8),int(huge(0_2)*2,kind=i8))
   lcg = int(mod(s,int(huge(0),kind=i8)),kind(0))

end function

!>
!! Set the seed in the built-in random_seed subroutine using the system clock
!! modified by lcg().
!! (Taken from FHI-aims with permission of copyright holders)
!!
subroutine elsi_init_random_seed()

   implicit none

   real(kind=r8) :: wtime
   integer(kind=i4) :: i
   integer(kind=i4) :: n
   integer(kind=i8) :: t

   integer(kind=i4), allocatable :: seed(:)

   character(len=*), parameter :: caller = "elsi_init_random_seed"

   call elsi_get_time(wtime)

   t = nint(wtime*1000.0_r8,kind=i8)

   call random_seed(size=n)

   allocate(seed(n))

   do i = 1,n
      seed(i) = lcg(t)
   end do

   call random_seed(put=seed)

   deallocate(seed)

end subroutine

end module ELSI_OUTPUT
