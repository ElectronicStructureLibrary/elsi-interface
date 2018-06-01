! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module performs IO to stdout and files.
!!
module ELSI_IO

   use ELSI_CONSTANTS, only: MULTI_PROC,SINGLE_PROC,BLACS_DENSE,PEXSI_CSC,&
                             SIESTA_CSC,UNSET,STR_LEN,UUID_LEN,ELPA_SOLVER,&
                             PEXSI_SOLVER,SIPS_SOLVER,OMM_SOLVER,DMP_SOLVER,&
                             DATETIME_LEN
   use ELSI_DATATYPE,  only: elsi_param_t,elsi_basic_t
   use ELSI_PRECISION, only: r8,i4,i8
   use FORTJSON,       only: fjson_write_name_value,fjson_reset_fj_handle,&
                             fjson_start_name_object,fjson_start_array,&
                             fjson_finish_object,fjson_open_file,fjson_handle,&
                             fjson_start_object,fjson_get_datetime_rfc3339,&
                             fjson_close_file,fjson_finish_array

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
!! This routine prints a message.
!!
subroutine elsi_say(bh,info_str)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   character(len=*),   intent(in) :: info_str

   character(len=40), parameter :: caller = "elsi_say"

   if(bh%print_info > 0) then
      write(bh%print_unit,"(A)") trim(info_str)
   endif

end subroutine

!>
!! This routine adds an entry to the log file.
!!
subroutine elsi_add_log(ph,bh,jh,dt0,t0,caller)

   implicit none

   type(elsi_param_t),          intent(in)    :: ph
   type(elsi_basic_t),          intent(inout) :: bh
   type(fjson_handle),          intent(inout) :: jh
   character(len=DATETIME_LEN), intent(in)    :: dt0
   real(kind=r8),               intent(in)    :: t0
   character(len=*),            intent(in)    :: caller

   real(kind=r8)               :: t1
   real(kind=r8)               :: t_total
   character(len=STR_LEN)      :: solver_tag
   character(len=STR_LEN)      :: elsi_tag
   character(len=STR_LEN)      :: user_tag
   character(len=DATETIME_LEN) :: dt_record

   if(bh%print_json > 0) then
      if(.not. bh%json_init) then
         call fjson_open_file(jh,66,"elsi_log.json")
         call fjson_start_array(jh)

         if(.not. bh%uuid_ready) then
            call elsi_gen_uuid(bh%uuid)
            bh%uuid_ready = .true.
         endif

         bh%json_init = .true.
      endif

      call elsi_get_time(t1)
      call fjson_get_datetime_rfc3339(dt_record)

      select case(ph%solver)
      case(ELPA_SOLVER)
         if(ph%parallel_mode == SINGLE_PROC) then
            solver_tag = "LAPACK"
         else
            solver_tag = "ELPA"
         endif
      case(OMM_SOLVER)
         if(ph%n_calls <= ph%omm_n_elpa) then
            solver_tag = "ELPA"
         else
            solver_tag = "LIBOMM"
         endif
      case(PEXSI_SOLVER)
         solver_tag = "PEXSI"
      case(SIPS_SOLVER)
         if(ph%n_calls <= ph%sips_n_elpa) then
            solver_tag = "ELPA"
         else
            solver_tag = "SIPS"
         endif
      case(DMP_SOLVER)
         solver_tag = "DMP"
      end select

      t_total  = t1-t0
      elsi_tag = adjustr(trim(solver_tag))
      user_tag = adjustr(trim(bh%user_tag))

      call fjson_start_object(jh)

      call elsi_print_versioning(bh%uuid,jh)
      call fjson_write_name_value(jh,"iteration",ph%n_calls)
      if(caller(6:6) == "e") then
         call fjson_write_name_value(jh,"output_type","EIGENSOLUTION")
      else
         call fjson_write_name_value(jh,"output_type","DENSITY MATRIX")
      endif
      if(caller(9:9) == "r") then
         call fjson_write_name_value(jh,"data_type","REAL")
      else
         call fjson_write_name_value(jh,"data_type","COMPLEX")
      endif
      call fjson_write_name_value(jh,"elsi_tag",elsi_tag)
      call fjson_write_name_value(jh,"user_tag",user_tag)
      call fjson_write_name_value(jh,"start_datetime",dt0)
      call fjson_write_name_value(jh,"record_datetime",dt_record)
      call fjson_write_name_value(jh,"total_time",t_total)

      call elsi_print_handle_summary(ph,bh,jh)

      if(ph%matrix_format == BLACS_DENSE) then
         call elsi_print_den_settings(bh,jh)
      else
         call elsi_print_csc_settings(bh,jh)
      endif

      select case(ph%solver)
      case(DMP_SOLVER)
         call elsi_print_dmp_settings(ph,jh)
      case(ELPA_SOLVER)
         call elsi_print_elpa_settings(ph,jh)
      case(OMM_SOLVER)
         call elsi_print_omm_settings(ph,jh)
      case(PEXSI_SOLVER)
         call elsi_print_pexsi_settings(ph,jh)
      case(SIPS_SOLVER)
         call elsi_print_sips_settings(ph,jh)
      end select

      call fjson_finish_object(jh)
   endif

end subroutine

!>
!! This routine prints the state of the handle.
!!
subroutine elsi_print_handle_summary(ph,bh,jh)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(elsi_basic_t), intent(in)    :: bh
   type(fjson_handle), intent(inout) :: jh

   real(kind=r8) :: sparsity

   character(len=40), parameter :: caller = "elsi_print_handle_summary"

   call fjson_write_name_value(jh,"n_electrons",ph%n_electrons)
   if(ph%parallel_mode == MULTI_PROC) then
      call fjson_write_name_value(jh,"n_spin",ph%n_spins)
      call fjson_write_name_value(jh,"n_kpts",ph%n_kpts)
   endif
   if(ph%solver == ELPA_SOLVER .or. ph%solver == SIPS_SOLVER) then
      call fjson_write_name_value(jh,"n_states",ph%n_states)
   endif
   if(ph%matrix_format == BLACS_DENSE) then
      call fjson_write_name_value(jh,"matrix_format","BLACS_DENSE")
   elseif(ph%matrix_format == PEXSI_CSC) then
      call fjson_write_name_value(jh,"matrix_format","PEXSI_CSC")
   elseif(ph%matrix_format == SIESTA_CSC) then
      call fjson_write_name_value(jh,"matrix_format","SIESTA_CSC")
   endif
   call fjson_write_name_value(jh,"n_basis",ph%n_basis)
   if(ph%parallel_mode == MULTI_PROC) then
      sparsity = 1.0_r8-(1.0_r8*bh%nnz_g/ph%n_basis/ph%n_basis)
      call fjson_write_name_value(jh,"sparsity",sparsity)
   endif
   if(ph%parallel_mode == MULTI_PROC) then
      call fjson_write_name_value(jh,"parallel_mode","MULTI_PROC")
   elseif(ph%parallel_mode == SINGLE_PROC) then
      call fjson_write_name_value(jh,"parallel_mode","SINGLE_PROC")
   endif
   call fjson_write_name_value(jh,"n_procs",bh%n_procs)
   if(ph%solver == ELPA_SOLVER) then
      call fjson_write_name_value(jh,"solver","ELPA")
   elseif(ph%solver == OMM_SOLVER) then
      call fjson_write_name_value(jh,"solver","libOMM")
   elseif(ph%solver == PEXSI_SOLVER) then
      call fjson_write_name_value(jh,"solver","PEXSI")
   elseif(ph%solver == SIPS_SOLVER) then
      call fjson_write_name_value(jh,"solver","SLEPc_SIPs")
   elseif(ph%solver == DMP_SOLVER) then
      call fjson_write_name_value(jh,"solver","DMP")
   endif

end subroutine

!>
!! This routine prints versioning information.
!!
subroutine elsi_print_versioning(uuid,jh)

   implicit none

   character(len=UUID_LEN), intent(in)    :: uuid
   type(fjson_handle),      intent(inout) :: jh

   logical            :: COMMIT_MODIFIED
   character(len=10)  :: DATE_STAMP
   character(len=40)  :: COMMIT
   character(len=8)   :: COMMIT_ABBREV
   character(len=40)  :: COMMIT_MSG_ABBREV
   character(len=40)  :: HOSTNAME
   character(len=20)  :: DATETIME

   character(len=40), parameter :: caller = "elsi_print_versioning"

   call elsi_version_info(DATE_STAMP,COMMIT,COMMIT_ABBREV,COMMIT_MODIFIED,&
           COMMIT_MSG_ABBREV,HOSTNAME,DATETIME)

   call fjson_write_name_value(jh,"data_source","ELSI")
   call fjson_write_name_value(jh,"date_stamp",trim(adjustl(DATE_STAMP)))
   call fjson_write_name_value(jh,"git_commit",trim(adjustl(COMMIT)))
   call fjson_write_name_value(jh,"git_commit_modified",COMMIT_MODIFIED)
   call fjson_write_name_value(jh,"git_message_abbrev",&
           trim(adjustl(COMMIT_MSG_ABBREV)))
   call fjson_write_name_value(jh,"source_created_on_hostname",&
           trim(adjustl(HOSTNAME)))
   call fjson_write_name_value(jh,"source_created_at_datetime",&
           trim(adjustl(DATETIME)))
   call fjson_write_name_value(jh,"uuid",trim(adjustl(uuid)))

end subroutine

!>
!! This routine prints out settings for DMP.
!!
subroutine elsi_print_dmp_settings(ph,jh)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(fjson_handle), intent(inout) :: jh

   character(len=40), parameter :: caller = "elsi_print_dmp_settings"

   call fjson_start_name_object(jh,"solver_settings")
   call fjson_write_name_value(jh,"dmp_n_states",ph%dmp_n_states)
   call fjson_write_name_value(jh,"dmp_method",ph%dmp_method)
   call fjson_write_name_value(jh,"dmp_max_power",ph%dmp_max_power)
   call fjson_write_name_value(jh,"dmp_max_iter",ph%dmp_max_iter)
   call fjson_write_name_value(jh,"dmp_tol",ph%dmp_tol)
   call fjson_finish_object(jh)

end subroutine

!>
!! This routine prints out settings for ELPA.
!!
subroutine elsi_print_elpa_settings(ph,jh)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(fjson_handle), intent(inout) :: jh

   character(len=40), parameter :: caller = "elsi_print_elpa_settings"

   call fjson_start_name_object(jh,"solver_settings")
   call fjson_write_name_value(jh,"elpa_solver",ph%elpa_solver)
   call fjson_write_name_value(jh,"elpa_n_states",ph%n_states)
   call fjson_write_name_value(jh,"elpa_gpu",ph%elpa_gpu)
   call fjson_write_name_value(jh,"elpa_gpu_kernels",ph%elpa_gpu_kernels)
   call fjson_finish_object(jh)

end subroutine

!>
!! This routine prints out settings for libOMM.
!!
subroutine elsi_print_omm_settings(ph,jh)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(fjson_handle), intent(inout) :: jh

   character(len=40), parameter :: caller = "elsi_print_omm_settings"

   call fjson_start_name_object(jh,"solver_settings")
   call fjson_write_name_value(jh,"omm_n_states",ph%omm_n_states)
   call fjson_write_name_value(jh,"omm_n_elpa",ph%omm_n_elpa)
   call fjson_write_name_value(jh,"omm_flavor",ph%omm_flavor)
   call fjson_write_name_value(jh,"omm_tol",ph%omm_tol)
   call fjson_finish_object(jh)

end subroutine

!>
!! This routine prints out settings for PEXSI.
!!
subroutine elsi_print_pexsi_settings(ph,jh)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(fjson_handle), intent(inout) :: jh

   character(len=40), parameter :: caller = "elsi_print_pexsi_settings"

   call fjson_start_name_object(jh,"solver_settings")
   call fjson_write_name_value(jh,"pexsi_n_pole",ph%pexsi_options%numPole)
   call fjson_write_name_value(jh,"pexsi_n_point",ph%pexsi_options%nPoints)
   call fjson_write_name_value(jh,"pexsi_np_per_pole",ph%pexsi_np_per_pole)
   call fjson_write_name_value(jh,"pexsi_np_per_point",ph%pexsi_np_per_point)
   call fjson_write_name_value(jh,"pexsi_n_prow_pexsi",ph%pexsi_n_prow)
   call fjson_write_name_value(jh,"pexsi_n_pcol_pexsi",ph%pexsi_n_pcol)
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
!! This routine prints out settings for SIPS.
!!
subroutine elsi_print_sips_settings(ph,jh)

   implicit none

   type(elsi_param_t), intent(in)    :: ph
   type(fjson_handle), intent(inout) :: jh

   character(len=40), parameter :: caller = "elsi_print_sips_settings"

   call fjson_start_name_object(jh,"solver_settings")
   call fjson_write_name_value(jh,"sips_n_states",ph%n_states)
   call fjson_write_name_value(jh,"sips_n_elpa",ph%sips_n_elpa)
   call fjson_write_name_value(jh,"sips_slice_type",ph%sips_slice_type)
   call fjson_write_name_value(jh,"sips_n_slices",ph%sips_n_slices)
   call fjson_write_name_value(jh,"sips_np_per_slice",ph%sips_np_per_slice)
   call fjson_write_name_value(jh,"sips_buffer",ph%sips_buffer)
   call fjson_finish_object(jh)

end subroutine

!>
!! This routine prints out settings for the dense matrix format.
!!
subroutine elsi_print_den_settings(bh,jh)

   implicit none

   type(elsi_basic_t), intent(in)    :: bh
   type(fjson_handle), intent(inout) :: jh

   character(len=40), parameter :: caller = "elsi_print_den_settings"

   call fjson_start_name_object(jh,"matrix_format_settings")
   call fjson_write_name_value(jh,"blk",bh%blk)
   call fjson_write_name_value(jh,"n_prow",bh%n_prow)
   call fjson_write_name_value(jh,"n_pcol",bh%n_pcol)
   call fjson_write_name_value(jh,"blacs_ready",bh%blacs_ready)
   call fjson_finish_object(jh)

end subroutine

!>
!! This routine prints out settings for the sparse matrix format.
!!
subroutine elsi_print_csc_settings(bh,jh)

   implicit none

   type(elsi_basic_t), intent(in)    :: bh
   type(fjson_handle), intent(inout) :: jh

   character(len=40), parameter :: caller = "elsi_print_csc_settings"

   call fjson_start_name_object(jh,"matrix_format_settings")
   call fjson_write_name_value(jh,"zero_def",bh%zero_def)
   call fjson_write_name_value(jh,"blk_sp2",bh%blk_sp2)
   call fjson_write_name_value(jh,"pexsi_csc_ready",bh%pexsi_csc_ready)
   call fjson_write_name_value(jh,"siesta_csc_ready",bh%siesta_csc_ready)
   call fjson_finish_object(jh)

end subroutine

!>
!! This routine prints a final output.
!!
subroutine elsi_final_print(ph,bh)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh

   real(kind=r8)      :: sparsity
   character(len=200) :: ll
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_final_print"

   write(ll,"(2X,A)") "|------------------------------------------------------"
   call elsi_say(bh,ll)

   write(info_str,"(2X,A)") "| Final ELSI Output"
   call elsi_say(bh,info_str)

   call elsi_say(bh,ll)

   write(info_str,"(2X,A)") "|"
   call elsi_say(bh,info_str)

   write(info_str,"(2X,A)") "| Physical Properties"
   call elsi_say(bh,info_str)

   write(info_str,"(2X,A,E22.8)") "|   Number of electrons       :",&
      ph%n_electrons
   call elsi_say(bh,info_str)

   if(ph%parallel_mode == MULTI_PROC) then
      write(info_str,"(2X,A,I22)") "|   Number of spins           :",ph%n_spins
      call elsi_say(bh,info_str)

      write(info_str,"(2X,A,I22)") "|   Number of k-points        :",ph%n_kpts
      call elsi_say(bh,info_str)
   endif

   if(ph%solver == ELPA_SOLVER .or. ph%solver == SIPS_SOLVER) then
      write(info_str,"(2X,A,I22)") "|   Number of states          :",ph%n_states
      call elsi_say(bh,info_str)
   endif

   write(info_str,"(2X,A)") "|"
   call elsi_say(bh,info_str)

   write(info_str,"(2X,A)") "| Matrix Properties"
   call elsi_say(bh,info_str)

   if(ph%matrix_format == BLACS_DENSE) then
      write(info_str,"(2X,A,A22)") "|   Matrix format             :",&
         "BLACS_DENSE"
   elseif(ph%matrix_format == PEXSI_CSC) then
      write(info_str,"(2X,A,A22)") "|   Matrix format             :",&
         "PEXSI_CSC"
   elseif(ph%matrix_format == SIESTA_CSC) then
      write(info_str,"(2X,A,A22)") "|   Matrix format             :",&
         "SIESTA_CSC"
   endif
   call elsi_say(bh,info_str)

   write(info_str,"(2X,A,I22)") "|   Number of basis functions :",ph%n_basis
   call elsi_say(bh,info_str)

   if(ph%parallel_mode == MULTI_PROC) then
      sparsity = 1.0_r8-(1.0_r8*bh%nnz_g/ph%n_basis/ph%n_basis)

      write(info_str,"(2X,A,E22.8)") "|   Matrix sparsity           :",sparsity
      call elsi_say(bh,info_str)
   endif

   write(info_str,"(2X,A)") "|"
   call elsi_say(bh,info_str)

   write(info_str,"(2X,A)") "| Computational Details"
   call elsi_say(bh,info_str)

   if(ph%parallel_mode == MULTI_PROC) then
      write(info_str,"(2X,A,A22)") "|   Parallel mode             :",&
         "MULTI_PROC"
   elseif(ph%parallel_mode == SINGLE_PROC) then
      write(info_str,"(2X,A,A22)") "|   Parallel mode             :",&
         "SINGLE_PROC"
   endif
   call elsi_say(bh,info_str)

   write(info_str,"(2X,A,I22)") "|   Number of MPI tasks       :",bh%n_procs
   call elsi_say(bh,info_str)

   if(ph%solver == ELPA_SOLVER) then
      write(info_str,"(2X,A,A22)") "|   Solver requested          :","ELPA"
   elseif(ph%solver == OMM_SOLVER) then
      write(info_str,"(2X,A,A22)") "|   Solver requested          :","libOMM"
   elseif(ph%solver == PEXSI_SOLVER) then
      write(info_str,"(2X,A,A22)") "|   Solver requested          :","PEXSI"
   elseif(ph%solver == SIPS_SOLVER) then
      write(info_str,"(2X,A,A22)") "|   Solver requested          :","SIPs"
   elseif(ph%solver == DMP_SOLVER) then
      write(info_str,"(2X,A,A22)") "|   Solver requested          :","DMP"
   endif
   call elsi_say(bh,info_str)

   write(info_str,"(2X,A,I22)") "|   Number of ELSI calls      :",ph%n_calls
   call elsi_say(bh,info_str)

   write(info_str,"(2X,A)") "|"
   call elsi_say(bh,info_str)

   call elsi_say(bh,ll)

   write(info_str,"(2X,A)") "| ELSI Project (c)  elsi-interchange.org"
   call elsi_say(bh,info_str)

   call elsi_say(bh,ll)

end subroutine

!>
!! This routine gets the current wallclock time.
!! (Taken from FHI-aims with permission of copyright holders)
!!
subroutine elsi_get_time(wtime)

   implicit none

   real(kind=r8), intent(out) :: wtime

   character(len=8)  :: cdate
   character(len=10) :: ctime
   integer(kind=i4)  :: val
   integer(kind=i4)  :: int_year
   real(kind=r8)     :: year
   real(kind=r8)     :: day
   real(kind=r8)     :: hour
   real(kind=r8)     :: minute
   real(kind=r8)     :: second
   real(kind=r8)     :: millisecond

   character(len=40), parameter :: caller = "elsi_get_time"

   call date_and_time(cdate,ctime)

   read(cdate(1:4),"(I4)") val

   int_year = val
   year     = real(val,kind=r8)-2009.0_r8 ! 2009 is an arbitrary zero
   day      = year*365+floor(year/4.0_r8)

   read(cdate(5:6),"(I2)") val

   val = val-1

   do while(val > 0)
      if(val == 1) then
         day = day+31
      elseif(val == 2) then
         if(mod(int_year,4) == 0) then
            day = day+29
         else
            day = day+28
         endif
      elseif(val == 3) then
         day = day+31
      elseif(val == 4) then
         day = day+30
      elseif(val == 5) then
         day = day+31
      elseif(val == 6) then
         day = day+30
      elseif(val == 7) then
         day = day+31
      elseif(val == 8) then
         day = day+31
      elseif(val == 9) then
         day = day+30
      elseif(val == 10) then
         day = day+31
      elseif(val == 11) then
         day = day+30
      endif

      val = val-1
   enddo

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

   wtime = day*24.0_r8*3600.0_r8+hour*3600.0_r8+minute*60.0_r8+second+&
              millisecond*0.001_r8

end subroutine

!>
!! Generate a UUID (unique identifier) in RFC 4122 format.
!! (Taken from FHI-aims with permission of copyright holders)
!!
subroutine elsi_gen_uuid(uuid)

   implicit none

   character(len=UUID_LEN), intent(out) :: uuid

   integer(kind=i4) :: ii3
   integer(kind=i4) :: ii4
   integer(kind=i4) :: i_entry
   real(kind=r8)    :: rr(8)
   character(len=3) :: ss3
   character(len=4) :: ss4

   character(len=40), parameter :: caller = "elsi_gen_uuid"

   call elsi_init_random_seed()

   ii3 = 4095
   ii4 = 65535

   call random_number(rr)

   do i_entry=1,8
      write(ss3,"(Z3.3)") transfer(int(rr(i_entry)*ii3),16)
      write(ss4,"(Z4.4)") transfer(int(rr(i_entry)*ii4),16)

      if(i_entry == 1) then
         write(uuid,"(A)") ss4
      elseif(i_entry == 2) then
         write(uuid,"(2A)") trim(uuid),ss4
      elseif(i_entry == 3) then
         write(uuid,"(3A)") trim(uuid),"-",ss4
      elseif(i_entry == 4) then
         write(uuid,"(3A)") trim(uuid),"-4",ss3
      elseif(i_entry == 5) then
         write(uuid,"(4A)") trim(uuid),"-A",ss3,"-"
      else
         write(uuid,"(2A)") trim(uuid),ss4
      endif
   enddo

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
   endif

   s   = mod(s*int(huge(0_2),kind=i8),int(huge(0_2)*2,kind=i8))
   lcg = int(mod(s,int(huge(0),kind=i8)),kind(0))

end function

!>
!! Set the seed in the built-in random_seed subroutine using the system clock
!! modified by lcg().
!! (Taken from FHI-aims with permission of copyright holders)
!!
subroutine elsi_init_random_seed()

   implicit none

   integer(kind=i4) :: i
   integer(kind=i4) :: n
   integer(kind=i4) :: dt(8)
   integer(kind=i8) :: t

   integer(kind=i4), allocatable :: seed(:)

   character(len=40), parameter :: caller = "elsi_init_random_seed"

   call random_seed(size=n)
   call date_and_time(values=dt)

   t = (dt(1)-1970)*365*24*60*60*1000+dt(2)*31*24*60*60*1000+&
          dt(3)*24*60*60*1000+dt(5)*60*60*1000+dt(6)*60*1000+dt(7)*1000+dt(8)

   allocate(seed(n))

   ! Writting the array with seeds
   do i = 1,n
      seed(i) = lcg(t)
   enddo

   call random_seed(put=seed)

   deallocate(seed)

end subroutine

end module ELSI_IO
