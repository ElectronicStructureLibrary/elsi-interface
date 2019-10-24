! Copyright (c) 2015-2019, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! Read input parameters from a file.
!!
module ELSI_INPUT

   use ELSI_DATATYPE, only: elsi_handle,elsi_basic_t
   use ELSI_MPI, only: elsi_stop
   use ELSI_PRECISION, only: r8,i4
   use ELSI_SET, only: elsi_set_output,elsi_set_output_log,elsi_set_save_ovlp,&
       elsi_set_unit_ovlp,elsi_set_zero_def,elsi_set_illcond_check,&
       elsi_set_illcond_tol,elsi_set_energy_gap,elsi_set_spectrum_width,&
       elsi_set_dimensionality,elsi_set_extrapolation,elsi_set_elpa_solver,&
       elsi_set_elpa_n_single,elsi_set_elpa_gpu,elsi_set_elpa_autotune,&
       elsi_set_omm_flavor,elsi_set_omm_n_elpa,elsi_set_omm_tol,&
       elsi_set_pexsi_n_mu,elsi_set_pexsi_n_pole,elsi_set_pexsi_np_per_pole,&
       elsi_set_pexsi_np_symbo,elsi_set_pexsi_inertia_tol,&
       elsi_set_eigenexa_method,elsi_set_sips_n_elpa,elsi_set_sips_n_slice,&
       elsi_set_sips_ev_min,elsi_set_sips_ev_max,elsi_set_ntpoly_method,&
       elsi_set_ntpoly_tol,elsi_set_ntpoly_filter,elsi_set_magma_solver,&
       elsi_set_mu_broaden_scheme,elsi_set_mu_broaden_width,elsi_set_mu_tol,&
       elsi_set_mu_mp_order
   use ELSI_UTIL, only: elsi_check_init

   implicit none

   private

   public :: elsi_set_input_file

contains

!>
!! Set runtime parameters from a file.
!!
subroutine elsi_set_input_file(eh,f_name)

   implicit none

   type(elsi_handle), intent(inout) :: eh !< Handle
   character(len=*), intent(in) :: f_name !< File name

   real(kind=r8) :: val_r8
   integer(kind=i4) :: i
   integer(kind=i4) :: ierr
   integer(kind=i4) :: val_i4
   logical :: kwd_found
   character(len=200) :: msg
   character(len=20) :: kwd
   character(len=20) :: val_str

   character(len=*), parameter :: caller = "elsi_set_input_file"

   call elsi_check_init(eh%bh,eh%handle_init,caller)

   open(313,file=f_name,status="old",iostat=ierr)

   if(ierr /= 0) then
      write(msg,"(2A)") "Failed to open input file ",trim(f_name)
      call elsi_stop(eh%bh,msg,caller)
   end if

   do
      ! Read in a line
      read(313,"(A)",iostat=ierr) msg

      if(ierr < 0) then
         ! EOF
         exit
      else if(ierr > 0) then
         write(msg,"(2A)") "Failed to parse input file ",trim(f_name)
         call elsi_stop(eh%bh,msg,caller)
      end if

      ! Remove leading blanks
      msg = adjustl(msg)

      ! Skip comment line
      if(msg(1:1) == "#") then
         cycle
      end if

      ! Parse
      kwd_found = .false.

      do i = 1,len(trim(msg))
         select case(msg(i:i))
         case("A":"Z")
            msg(i:i) = achar(iachar(msg(i:i))+32)
            kwd_found = .true.
         case("a":"z")
            kwd_found = .true.
         case("0":"9","+","-",".","_")
            if(.not. kwd_found) then
               msg(i:i) = " "
            end if
         case default
            msg(i:i) = " "
         end select
      end do

      read(msg,*,iostat=ierr) kwd

      call elsi_check_read(eh%bh,ierr,kwd)

      select case(kwd)
      case("output")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_output(eh,val_i4)
      case("output_log")
         read(msg,*,iostat=ierr) kwd,val_str

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_str_to_int(val_str,val_i4)
         call elsi_set_output_log(eh,val_i4)
      case("save_ovlp")
         read(msg,*,iostat=ierr) kwd,val_str

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_str_to_int(val_str,val_i4)
         call elsi_set_save_ovlp(eh,val_i4)
      case("unit_ovlp")
         read(msg,*,iostat=ierr) kwd,val_str

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_str_to_int(val_str,val_i4)
         call elsi_set_unit_ovlp(eh,val_i4)
      case("zero_def")
         read(msg,*,iostat=ierr) kwd,val_r8

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_zero_def(eh,val_r8)
      case("illcond_check")
         read(msg,*,iostat=ierr) kwd,val_str

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_str_to_int(val_str,val_i4)
         call elsi_set_illcond_check(eh,val_i4)
      case("illcond_tol")
         read(msg,*,iostat=ierr) kwd,val_r8

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_illcond_tol(eh,val_r8)
      case("energy_gap")
         read(msg,*,iostat=ierr) kwd,val_r8

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_energy_gap(eh,val_r8)
      case("spectrum_width")
         read(msg,*,iostat=ierr) kwd,val_r8

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_spectrum_width(eh,val_r8)
      case("dimensionality")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_dimensionality(eh,val_i4)
      case("extrapolation")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_extrapolation(eh,val_i4)
      case("elpa_solver")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_elpa_solver(eh,val_i4)
      case("elpa_n_single")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_elpa_n_single(eh,val_i4)
      case("elpa_gpu")
         read(msg,*,iostat=ierr) kwd,val_str

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_str_to_int(val_str,val_i4)
         call elsi_set_elpa_gpu(eh,val_i4)
      case("elpa_autotune")
         read(msg,*,iostat=ierr) kwd,val_str

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_str_to_int(val_str,val_i4)
         call elsi_set_elpa_autotune(eh,val_i4)
      case("omm_flavor")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_omm_flavor(eh,val_i4)
      case("omm_n_elpa")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_omm_n_elpa(eh,val_i4)
      case("omm_tol")
         read(msg,*,iostat=ierr) kwd,val_r8

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_omm_tol(eh,val_r8)
      case("pexsi_n_mu")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_pexsi_n_mu(eh,val_i4)
      case("pexsi_n_pole")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_pexsi_n_pole(eh,val_i4)
      case("pexsi_np_per_pole")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_pexsi_np_per_pole(eh,val_i4)
      case("pexsi_np_symbo")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_pexsi_np_symbo(eh,val_i4)
      case("pexsi_inertia_tol")
         read(msg,*,iostat=ierr) kwd,val_r8

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_pexsi_inertia_tol(eh,val_r8)
      case("eigenexa_method")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_eigenexa_method(eh,val_i4)
      case("sips_n_elpa")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_sips_n_elpa(eh,val_i4)
      case("sips_n_slice")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_sips_n_slice(eh,val_i4)
      case("sips_ev_min")
         read(msg,*,iostat=ierr) kwd,val_r8

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_sips_ev_min(eh,val_r8)
      case("sips_ev_max")
         read(msg,*,iostat=ierr) kwd,val_r8

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_sips_ev_max(eh,val_r8)
      case("ntpoly_method")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_ntpoly_method(eh,val_i4)
      case("ntpoly_tol")
         read(msg,*,iostat=ierr) kwd,val_r8

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_ntpoly_tol(eh,val_r8)
      case("ntpoly_filter")
         read(msg,*,iostat=ierr) kwd,val_r8

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_ntpoly_filter(eh,val_r8)
      case("magma_solver")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_magma_solver(eh,val_i4)
      case("mu_broaden_scheme")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_mu_broaden_scheme(eh,val_i4)
      case("mu_broaden_width")
         read(msg,*,iostat=ierr) kwd,val_r8

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_mu_broaden_width(eh,val_r8)
      case("mu_tol")
         read(msg,*,iostat=ierr) kwd,val_r8

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_mu_tol(eh,val_r8)
      case("mu_mp_order")
         read(msg,*,iostat=ierr) kwd,val_i4

         call elsi_check_read(eh%bh,ierr,kwd)
         call elsi_set_mu_mp_order(eh,val_i4)
      end select
   end do

   close(313)

end subroutine

!>
!! Check if a Fortran read statement is successful.
!!
subroutine elsi_check_read(bh,ierr,kwd)

   implicit none

   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: ierr
   character(len=*), intent(in) :: kwd

   character(len=200) :: msg

   character(len=*), parameter :: caller = "elsi_check_read"

   if(ierr /= 0) then
      write(msg,"(2A)") "Failed to parse input keyword ",trim(kwd)
      call elsi_stop(bh,msg,caller)
   end if

end subroutine

!>
!! Convert a string to an integer. Return 0 if the string itself is "0" or
!! contains "false" as a substring. Return 1 otherwise.
!!
subroutine elsi_str_to_int(val_str,val_i4)

   implicit none

   character(len=*), intent(in) :: val_str
   integer(kind=i4), intent(out) :: val_i4

   character(len=*), parameter :: caller = "elsi_str_to_int"

   select case(trim(adjustl(val_str)))
   case("false",".false.","no","0")
      val_i4 = 0
   case default
      val_i4 = 1
   end select

end subroutine

end module ELSI_INPUT
