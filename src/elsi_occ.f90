! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module provides routines for chemical potential determination.
!!
module ELSI_OCC

   use ELSI_CONSTANTS, only: GAUSSIAN,FERMI,METHFESSEL_PAXTON,COLD,CUBIC,&
       SQRT_PI,INVERT_SQRT_PI
   use ELSI_DATATYPE, only: elsi_handle,elsi_param_t,elsi_basic_t
   use ELSI_IO, only: elsi_say
   use ELSI_MALLOC, only: elsi_allocate,elsi_deallocate
   use ELSI_MPI, only: elsi_stop
   use ELSI_PRECISION, only: r8,i4

   implicit none

   private

   public :: elsi_mu_and_occ
   public :: elsi_entropy
   public :: elsi_compute_mu_and_occ
   public :: elsi_compute_entropy

contains

!>
!! This routine computes the chemical potential and occupation numbers.
!!
subroutine elsi_compute_mu_and_occ(eh,n_electron,n_state,n_spin,n_kpt,&
              k_weights,evals,occ_nums,mu)

   implicit none

   type(elsi_handle), intent(in) :: eh !< Handle
   real(kind=r8), intent(in) :: n_electron !< Number of electrons
   integer(kind=i4), intent(in) :: n_state !< Number of states
   integer(kind=i4), intent(in) :: n_spin !< Number of spins
   integer(kind=i4), intent(in) :: n_kpt !< Number of k-points
   real(kind=r8), intent(in) :: k_weights(n_kpt) !< K-points weights
   real(kind=r8), intent(in) :: evals(n_state,n_spin,n_kpt) !< Eigenvalues
   real(kind=r8), intent(out) :: occ_nums(n_state,n_spin,n_kpt) !< Occupation members
   real(kind=r8), intent(out) :: mu !< Chemical potential

   character(len=*), parameter :: caller = "elsi_compute_mu_and_occ"

   call elsi_mu_and_occ(eh%ph,eh%bh,n_electron,n_state,n_spin,n_kpt,k_weights,&
           evals,occ_nums,mu)

end subroutine

!>
!! This routine computes the chemical potential and occupation numbers.
!!
subroutine elsi_mu_and_occ(ph,bh,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: n_electron
   integer(kind=i4), intent(in) :: n_state
   integer(kind=i4), intent(in) :: n_spin
   integer(kind=i4), intent(in) :: n_kpt
   real(kind=r8), intent(in) :: k_weights(n_kpt)
   real(kind=r8), intent(in) :: evals(n_state,n_spin,n_kpt)
   real(kind=r8), intent(out) :: occ_nums(n_state,n_spin,n_kpt)
   real(kind=r8), intent(out) :: mu

   real(kind=r8) :: e_low
   real(kind=r8) :: e_high
   real(kind=r8) :: mu_lower
   real(kind=r8) :: mu_upper
   real(kind=r8) :: diff_ne_lower ! Electron count error on lower bound
   real(kind=r8) :: diff_ne_upper ! Electron count error on upper bound
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: n_steps

   character(len=*), parameter :: caller = "elsi_mu_and_occ"

   ! Determine smallest and largest eivenvalues
   e_low = evals(1,1,1)
   e_high = evals(n_state,1,1)

   do i_kpt = 1,n_kpt
      do i_spin = 1,n_spin
         do i_state = 1,n_state
            if(evals(i_state,i_spin,i_kpt) < e_low) then
               e_low = evals(i_state,i_spin,i_kpt)
            end if
            if(evals(i_state,i_spin,i_kpt) > e_high) then
               e_high = evals(i_state,i_spin,i_kpt)
            end if
         end do
      end do
   end do

   ! Determine upper and lower bounds of mu
   mu_lower = e_low
   mu_upper = e_high

   if(mu_upper - mu_lower < ph%mu_tol) then
      mu_lower = mu_lower-1.0_r8
      mu_upper = mu_upper+1.0_r8
   end if

   occ_nums = 0.0_r8

   ! Compute electron count error
   call elsi_check_electrons(ph,n_electron,n_state,n_spin,n_kpt,k_weights,&
           evals,occ_nums,mu_lower,diff_ne_lower)
   call elsi_check_electrons(ph,n_electron,n_state,n_spin,n_kpt,k_weights,&
           evals,occ_nums,mu_upper,diff_ne_upper)

   ! Enlarge the interval towards both sides if solution not found
   n_steps = 0

   do while(diff_ne_lower*diff_ne_upper > 0.0_r8)
      n_steps = n_steps+1

      if(n_steps > ph%mu_max_steps) then
         call elsi_stop(bh,"Chemical potential not found.",caller)
      end if

      mu_lower = mu_lower-0.5_r8*abs(e_high-e_low)
      mu_upper = mu_upper+0.5_r8*abs(e_high-e_low)

      call elsi_check_electrons(ph,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_lower,diff_ne_lower)
      call elsi_check_electrons(ph,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_upper,diff_ne_upper)
   end do

   ! Perform bisection
   call elsi_find_mu(ph,bh,n_electron,n_state,n_spin,n_kpt,k_weights,evals,&
           occ_nums,mu_lower,mu_upper,mu)

end subroutine

!>
!! This routine computes the number of electrons for a given chemical potential,
!! and returns the error in the number of electrons. The occupation numbers are
!! updated as well.
!!
subroutine elsi_check_electrons(ph,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_in,diff_ne_out)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   real(kind=r8), intent(in) :: n_electron
   integer(kind=i4), intent(in) :: n_state
   integer(kind=i4), intent(in) :: n_spin
   integer(kind=i4), intent(in) :: n_kpt
   real(kind=r8), intent(in) :: k_weights(n_kpt)
   real(kind=r8), intent(in) :: evals(n_state,n_spin,n_kpt)
   real(kind=r8), intent(out) :: occ_nums(n_state,n_spin,n_kpt)
   real(kind=r8), intent(in) :: mu_in
   real(kind=r8), intent(out) :: diff_ne_out

   real(kind=r8) :: spin_degen
   real(kind=r8) :: invert_width
   real(kind=r8) :: delta
   real(kind=r8) :: max_exp ! Maximum possible exponent
   real(kind=r8) :: arg
   real(kind=r8) :: weight
   real(kind=r8) :: A
   real(kind=r8) :: H_even
   real(kind=r8) :: H_odd
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: i_mp

   character(len=*), parameter :: caller = "elsi_check_electrons"

   invert_width = 1.0_r8/ph%mu_width
   diff_ne_out = 0.0_r8

   if(.not. ph%spin_is_set) then
      if(n_spin == 2) then
         spin_degen = 1.0_r8
      else
         spin_degen = 2.0_r8
      end if
   else
      spin_degen = ph%spin_degen
   end if

   select case(ph%mu_scheme)
   case(GAUSSIAN)
      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               occ_nums(i_state,i_spin,i_kpt) = spin_degen*0.5_r8*&
                  (1.0_r8-erf((evals(i_state,i_spin,i_kpt)-mu_in)*invert_width))

               diff_ne_out = diff_ne_out+occ_nums(i_state,i_spin,i_kpt)*&
                                k_weights(i_kpt)
            end do
         end do
      end do
   case(FERMI)
      max_exp = maxexponent(mu_in)*log(2.0_r8)

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg = (evals(i_state,i_spin,i_kpt)-mu_in)*invert_width

               if(arg < max_exp) then
                  occ_nums(i_state,i_spin,i_kpt) = spin_degen/(1.0_r8+exp(arg))

                  diff_ne_out = diff_ne_out+occ_nums(i_state,i_spin,i_kpt)*&
                                   k_weights(i_kpt)
               else
                  occ_nums(i_state,i_spin,i_kpt) = 0.0_r8
               end if
            end do
         end do
      end do
   case(METHFESSEL_PAXTON)
      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg = (evals(i_state,i_spin,i_kpt)-mu_in)*invert_width
               weight = exp(-arg**2)

               occ_nums(i_state,i_spin,i_kpt) = 0.5_r8*(1.0_r8-erf(arg))*&
                                                   spin_degen

               if(ph%mu_mp_order > 0) then
                  A = -0.25_r8*INVERT_SQRT_PI
                  H_even = 1.0_r8
                  H_odd = 2.0_r8*arg

                  occ_nums(i_state,i_spin,i_kpt) =&
                     occ_nums(i_state,i_spin,i_kpt)+A*H_odd*weight*spin_degen
               end if

               if(ph%mu_mp_order > 1) then
                  do i_mp = 2,ph%mu_mp_order
                     A = -0.25_r8/real(i_mp,kind=r8)*A
                     H_even = 2.0_r8*arg*H_odd-2.0_r8*real(i_mp,kind=r8)*H_even
                     H_odd = 2.0_r8*arg*H_even-2.0_r8*real(i_mp+1,kind=r8)*H_odd

                     occ_nums(i_state,i_spin,i_kpt) =&
                        occ_nums(i_state,i_spin,i_kpt)+A*H_odd*weight*spin_degen
                  end do
               end if

               diff_ne_out = diff_ne_out+occ_nums(i_state,i_spin,i_kpt)*&
                                k_weights(i_kpt)
            end do
         end do
      end do
   case(CUBIC)
      ! To have a consistent slope of the occupation function at the chemical
      ! potential, the parameters for GAUSSIAN and CUBIC should be related as:
      delta = 0.75_r8*SQRT_PI*ph%mu_width

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg = (evals(i_state,i_spin,i_kpt)-mu_in)/delta

               if(arg <= -1.0_r8) then
                  occ_nums(i_state,i_spin,i_kpt) = spin_degen
               else if(arg >= 1.0_r8) then
                  occ_nums(i_state,i_spin,i_kpt) = 0.0_r8
               else
                  occ_nums(i_state,i_spin,i_kpt) = spin_degen*0.25_r8*&
                     (arg+2.0_r8)*(arg-1.0_r8)**2
               end if

               diff_ne_out = diff_ne_out+occ_nums(i_state,i_spin,i_kpt)*&
                                k_weights(i_kpt)
            end do
         end do
      end do
   case(COLD)
      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg = (evals(i_state,i_spin,i_kpt)-mu_in)*invert_width
               arg = -arg-sqrt(0.5_r8)

               occ_nums(i_state,i_spin,i_kpt) = (0.5_r8-erf(arg)*0.5_r8-&
                  INVERT_SQRT_PI*sqrt(0.5_r8)*exp(-arg**2))*spin_degen

               diff_ne_out = diff_ne_out+occ_nums(i_state,i_spin,i_kpt)*&
                                k_weights(i_kpt)
            end do
         end do
      end do
   end select

   diff_ne_out = diff_ne_out-n_electron

end subroutine

!>
!! This routine computes the chemical potential using a bisection algorithm.
!!
subroutine elsi_find_mu(ph,bh,n_electron,n_state,n_spin,n_kpt,k_weights,evals,&
              occ_nums,mu_lower_in,mu_upper_in,mu_out)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   real(kind=r8), intent(in) :: n_electron
   integer(kind=i4), intent(in) :: n_state
   integer(kind=i4), intent(in) :: n_spin
   integer(kind=i4), intent(in) :: n_kpt
   real(kind=r8), intent(in) :: k_weights(n_kpt)
   real(kind=r8), intent(in) :: evals(n_state,n_spin,n_kpt)
   real(kind=r8), intent(out) :: occ_nums(n_state,n_spin,n_kpt)
   real(kind=r8), intent(in) :: mu_lower_in
   real(kind=r8), intent(in) :: mu_upper_in
   real(kind=r8), intent(out) :: mu_out

   real(kind=r8) :: mu_left
   real(kind=r8) :: mu_right
   real(kind=r8) :: mu_mid
   real(kind=r8) :: diff_left ! Electron count error on left bound
   real(kind=r8) :: diff_right ! Electron count error on right bound
   real(kind=r8) :: diff_mid ! Electron count error on middle point
   logical :: found_mu
   integer(kind=i4) :: n_steps
   character(len=200) :: info_str

   character(len=*), parameter :: caller = "elsi_find_mu"

   n_steps = 0
   found_mu = .false.
   mu_left = mu_lower_in
   mu_right = mu_upper_in

   call elsi_check_electrons(ph,n_electron,n_state,n_spin,n_kpt,k_weights,&
           evals,occ_nums,mu_left,diff_left)
   call elsi_check_electrons(ph,n_electron,n_state,n_spin,n_kpt,k_weights,&
           evals,occ_nums,mu_right,diff_right)

   if(abs(diff_left) < ph%mu_tol) then
      mu_out = mu_left
      found_mu = .true.
   else if(abs(diff_right) < ph%mu_tol) then
      mu_out = mu_right
      found_mu = .true.
   end if

   do while(.not. found_mu .and. n_steps < ph%mu_max_steps)
      n_steps = n_steps+1
      mu_mid = 0.5_r8*(mu_left+mu_right)

      call elsi_check_electrons(ph,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_mid,diff_mid)

      if(abs(diff_mid) < ph%mu_tol) then
         mu_out = mu_mid
         found_mu = .true.
      else if(diff_mid < 0.0_r8) then
         mu_left = mu_mid
      else if(diff_mid > 0.0_r8) then
         mu_right = mu_mid
      end if
   end do

   if(found_mu) then
      call elsi_check_electrons(ph,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_out,diff_right)
   else
      ! Use mu of the right bound...
      call elsi_check_electrons(ph,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_right,diff_right)

      mu_out = mu_right

      ! ...with adjusted occupation numbers
      write(info_str,"(2X,A)")&
         "Chemical potential cannot reach the required accuracy."
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A,E10.2,A)") "| Residual error :",diff_right
      call elsi_say(bh,info_str)
      write(info_str,"(2X,A)")&
         "The error will be removed from the highest occupied states."
      call elsi_say(bh,info_str)

      call elsi_adjust_occ(ph,bh,n_state,n_spin,n_kpt,k_weights,evals,occ_nums,&
              diff_right)
   end if

end subroutine

!>
!! This routine cancels the small error in number of electrons.
!!
subroutine elsi_adjust_occ(ph,bh,n_state,n_spin,n_kpt,k_weights,evals,occ_nums,&
              diff_ne)

   implicit none

   type(elsi_param_t), intent(in) :: ph
   type(elsi_basic_t), intent(in) :: bh
   integer(kind=i4), intent(in) :: n_state
   integer(kind=i4), intent(in) :: n_spin
   integer(kind=i4), intent(in) :: n_kpt
   real(kind=r8), intent(in) :: k_weights(n_kpt)
   real(kind=r8), intent(in) :: evals(n_state,n_spin,n_kpt)
   real(kind=r8), intent(inout) :: occ_nums(n_state,n_spin,n_kpt)
   real(kind=r8), intent(inout) :: diff_ne

   real(kind=r8) :: min_eval
   integer(kind=i4) :: max_id
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: i_val
   integer(kind=i4) :: n_total

   real(kind=r8), allocatable :: eval_aux(:)

   character(len=*), parameter :: caller = "elsi_adjust_occ"

   n_total = n_state*n_spin*n_kpt

   call elsi_allocate(bh,eval_aux,n_total,"eval_aux",caller)

   ! Put evals into a 1D array
   i_val = 0

   do i_kpt = 1,n_kpt
      do i_spin = 1,n_spin
         do i_state = 1,n_state
            i_val = i_val+1
            eval_aux(i_val) = evals(i_state,i_spin,i_kpt)
         end do
      end do
   end do

   min_eval = minval(eval_aux,1)

   ! Remove error
   do i_val = 1,n_total
      max_id = maxloc(eval_aux,1)
      eval_aux(max_id) = min_eval-1.0_r8

      i_kpt = (max_id-1)/(n_spin*n_state)+1
      i_spin = mod((max_id-1)/n_state,n_spin)+1
      i_state = mod(max_id-1,n_state)+1

      if(k_weights(i_kpt)*occ_nums(i_state,i_spin,i_kpt) > diff_ne) then
         occ_nums(i_state,i_spin,i_kpt) = occ_nums(i_state,i_spin,i_kpt)-&
                                             diff_ne/k_weights(i_kpt)
         diff_ne = 0.0_r8
      else
         diff_ne = diff_ne-k_weights(i_kpt)*occ_nums(i_state,i_spin,i_kpt)
         occ_nums(i_state,i_spin,i_kpt) = 0.0_r8
      end if

      if(diff_ne <= ph%mu_tol) then
         exit
      end if
   end do

   call elsi_deallocate(bh,eval_aux,"eval_aux")

end subroutine

!>
!! This routine computes the entropy.
!!
subroutine elsi_compute_entropy(eh,n_state,n_spin,n_kpt,k_weights,evals,&
              occ_nums,mu,ts)

   implicit none

   type(elsi_handle), intent(in) :: eh !< Handle
   integer(kind=i4), intent(in) :: n_state !< Number of states
   integer(kind=i4), intent(in) :: n_spin !< Number of spins
   integer(kind=i4), intent(in) :: n_kpt !< Number of k-points
   real(kind=r8), intent(in) :: k_weights(n_kpt) !< K-points weights
   real(kind=r8), intent(in) :: evals(n_state,n_spin,n_kpt) !< Eigenvalues
   real(kind=r8), intent(in) :: occ_nums(n_state,n_spin,n_kpt) !< Occupation numbers
   real(kind=r8), intent(in) :: mu !< Input chemical potential
   real(kind=r8), intent(out) :: ts !< Entropy

   character(len=*), parameter :: caller = "elsi_compute_entropy"

   call elsi_entropy(eh%ph,n_state,n_spin,n_kpt,k_weights,evals,occ_nums,mu,ts)

end subroutine

!>
!! This routine computes the entropy.
!!
subroutine elsi_entropy(ph,n_state,n_spin,n_kpt,k_weights,evals,occ_nums,mu,ts)

   implicit none

   type(elsi_param_t), intent(in) :: ph !< Parameters
   integer(kind=i4), intent(in) :: n_state !< Number of states
   integer(kind=i4), intent(in) :: n_spin !< Number of spins
   integer(kind=i4), intent(in) :: n_kpt !< Number of k-points
   real(kind=r8), intent(in) :: k_weights(n_kpt) !< K-points weights
   real(kind=r8), intent(in) :: evals(n_state,n_spin,n_kpt) !< Eigenvalues
   real(kind=r8), intent(in) :: occ_nums(n_state,n_spin,n_kpt) !< Occupation numbers
   real(kind=r8), intent(in) :: mu !< Input chemical potential
   real(kind=r8), intent(out) :: ts !< Entropy

   real(kind=r8) :: spin_degen
   real(kind=r8) :: invert_width
   real(kind=r8) :: delta
   real(kind=r8) :: pre
   real(kind=r8) :: arg
   real(kind=r8) :: weight
   real(kind=r8) :: A
   real(kind=r8) :: H_even
   real(kind=r8) :: H_odd
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: i_mp

   real(kind=r8), parameter :: ts_thr = 1.0e-15_r8
   character(len=*), parameter :: caller = "elsi_entropy"

   invert_width = 1.0_r8/ph%mu_width
   ts = 0.0_r8
   pre = 0.0_r8

   if(.not. ph%spin_is_set) then
      if(n_spin == 2) then
         spin_degen = 1.0_r8
      else
         spin_degen = 2.0_r8
      end if
   else
      spin_degen = ph%spin_degen
   end if

   select case(ph%mu_scheme)
   case(GAUSSIAN)
      pre = 0.5_r8*spin_degen*ph%mu_width*INVERT_SQRT_PI

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg = (evals(i_state,i_spin,i_kpt)-mu)*invert_width
               ts = ts+exp(-arg**2)*k_weights(i_kpt)
            end do
         end do
      end do
   case(FERMI)
      pre = spin_degen*ph%mu_width

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg = occ_nums(i_state,i_spin,i_kpt)/spin_degen

               if(1.0_r8-arg > ts_thr .and. arg > ts_thr) then
                  ts = ts-(arg*log(arg)+(1.0_r8-arg)*log(1.0_r8-arg))*&
                          k_weights(i_kpt)
               end if
            end do
         end do
      end do
   case(METHFESSEL_PAXTON)
      pre = 0.5_r8*spin_degen*ph%mu_width

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg = (evals(i_state,i_spin,i_kpt)-mu)*invert_width
               weight = exp(-arg**2)
               A = INVERT_SQRT_PI
               H_even = 1.0_r8
               H_odd = 2.0_r8*arg
               ts = ts+INVERT_SQRT_PI*weight*k_weights(i_kpt)

               do i_mp = 1,ph%mu_mp_order
                  A = -0.25_r8/real(i_mp,kind=r8)*A
                  H_even = 2.0_r8*arg*H_odd-2.0_r8*real(i_mp,kind=r8)*H_even
                  H_odd = 2.0_r8*arg*H_even-2.0_r8*real(i_mp+1,kind=r8)*H_odd
                  ts = ts+A*H_even*weight*k_weights(i_kpt)
               end do
            end do
         end do
      end do
   case(CUBIC)
      delta = 0.75_r8*ph%mu_width*SQRT_PI
      pre = 0.1875_r8*spin_degen*ph%mu_width

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg = (evals(i_state,i_spin,i_kpt)-mu)/delta

               if(arg > -1.0_r8 .and. arg < 1.0_r8) then
                  ts = ts+(((arg**2)-1.0_r8)**2)*k_weights(i_kpt)
               end if
            end do
         end do
      end do
   case(COLD)
      pre = 0.5_r8*spin_degen*ph%mu_width*INVERT_SQRT_PI

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg = (evals(i_state,i_spin,i_kpt)-mu)*invert_width
               ts = ts+exp(-(arg-sqrt(0.5_r8))**2)*(1.0_r8-sqrt(2.0_r8)*arg)
            end do
         end do
      end do
   end select

   ts = pre*ts

end subroutine

end module ELSI_OCC
