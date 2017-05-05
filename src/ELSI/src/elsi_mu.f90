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
!! This module provides routines for chemical potential determination.
!!

module ELSI_MU

   use ELSI_PRECISION, only : dp
   use ELSI_CONSTANTS, only : INVERT_SQRT_PI, GAUSSIAN, FERMI, METHFESSEL_PAXTON_0, METHFESSEL_PAXTON_1
   use ELSI_DIMENSIONS
   use ELSI_UTILS

   implicit none
   private

   !> Public routines
   public  :: elsi_compute_mu_and_occ

   real(kind=dp) :: n_electron
   integer       :: n_state
   integer       :: n_spin
   integer       :: n_kpoint

contains

!>
!! This routine computes the chemical potential and occupation numbers
!! from eigenvalues.
!!
subroutine elsi_compute_mu_and_occ(n_electron_in,n_state_in,n_spin_in,n_kpoint_in,&
                                   kpoint_weights,eigenvalues,occ_numbers,mu)

   implicit none

   real(kind=dp), intent(in)  :: n_electron_in
   integer,       intent(in)  :: n_state_in
   integer,       intent(in)  :: n_spin_in
   integer,       intent(in)  :: n_kpoint_in
   real(kind=dp), intent(in)  :: kpoint_weights(n_kpoint_in)
   real(kind=dp), intent(in)  :: eigenvalues(n_state_in,n_spin_in,n_kpoint_in)
   real(kind=dp), intent(out) :: occ_numbers(n_state_in,n_spin_in,n_kpoint_in)
   real(kind=dp), intent(out) :: mu

   real(kind=dp) :: e_low         !< Lowest eigenvalue
   real(kind=dp) :: e_high        !< Highest eigenvalue
   real(kind=dp) :: mu_lower      !< Lower bound of chemical potential
   real(kind=dp) :: mu_upper      !< Upper bound of chemical potential
   real(kind=dp) :: diff_ne_lower !< Difference in number of electrons on lower bound
   real(kind=dp) :: diff_ne_upper !< Difference in number of electrons on upper bound

   integer :: i_state  !< State index
   integer :: i_kpoint !< K-point index
   integer :: i_spin   !< Spin index
   integer :: n_steps  !< Number of steps to find chemical potential interval

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_compute_mu_and_occ"

   n_electron = n_electron_in
   n_state    = n_state_in
   n_spin     = n_spin_in
   n_kpoint   = n_kpoint_in

   ! Determine the smallest and largest eivenvalues
   e_low = eigenvalues(1,1,1)
   e_high = eigenvalues(n_state,1,1)

   do i_kpoint = 1,n_kpoint
      do i_spin = 1,n_spin
         do i_state = 1,n_state
            if(eigenvalues(i_state,i_spin,i_kpoint) < e_low) then
               e_low = eigenvalues(i_state,i_spin,i_kpoint)
            endif
            if(eigenvalues(i_state,i_spin,i_kpoint) > e_high) then
               e_high = eigenvalues(i_state,i_spin,i_kpoint)
            endif
         enddo
      enddo
   enddo

   ! Determine the upper and lower bounds for chemical potential
   mu_lower = e_low

   if(e_low == e_high) then
      mu_upper = 0.0_dp
   else
      mu_upper = e_high
   endif

   occ_numbers = 0.0_dp

   ! Compute the difference of number of electrons
   call elsi_check_electrons(kpoint_weights,eigenvalues,occ_numbers,&
                             mu_lower,diff_ne_lower)
   call elsi_check_electrons(kpoint_weights,eigenvalues,occ_numbers,&
                             mu_upper,diff_ne_upper)

   ! If diff_ne_lower*diff_ne_upper > 0, it means that the solution is
   ! not in this interval.
   ! Enlarge the interval towards both sides, then recheck the condition.
   n_steps = 0
   do while(diff_ne_lower*diff_ne_upper > 0)
      n_steps = n_steps+1
      if(n_steps > max_mu_steps) then
         write(info_str,"(A,I13,A)") " Chemical potential not found in ",&
            max_mu_steps," iterations! Exiting..."
         call elsi_stop(info_str,caller)
      endif

      mu_lower = mu_lower-0.5_dp*ABS(e_high-e_low)
      mu_upper = mu_upper+0.5_dp*ABS(e_high-e_low)

      call elsi_check_electrons(kpoint_weights,eigenvalues,occ_numbers,&
                                mu_lower,diff_ne_lower)
      call elsi_check_electrons(kpoint_weights,eigenvalues,occ_numbers,&
                                mu_upper,diff_ne_upper)
   enddo

   ! At this point we should have the correct interval for chemical potential.
   ! Use simple bisection algorithm to find the solution.
   call elsi_find_mu(kpoint_weights,eigenvalues,occ_numbers,&
                     mu_lower,mu_upper,mu)

end subroutine

!>
!! This routine computes the number of electrons using a given chemical potential,
!! and returns the difference in number of electrons. The occupation numbers will
!! be updated as well.
!!
subroutine elsi_check_electrons(kpoint_weights,eigenvalues,occ_numbers,&
                                mu_in,diff_ne_out)

   implicit none

   real(kind=dp), intent(in)  :: kpoint_weights(n_kpoint)
   real(kind=dp), intent(in)  :: eigenvalues(n_state,n_spin,n_kpoint)
   real(kind=dp), intent(out) :: occ_numbers(n_state,n_spin,n_kpoint)
   real(kind=dp), intent(in)  :: mu_in       !< Input chemical potential
   real(kind=dp), intent(out) :: diff_ne_out !< Difference in number of electrons

   real(kind=dp) :: invert_width !< 1/broaden_width
   integer       :: n_steps      !< Number of steps to find chemical potential interval
   real(kind=dp) :: max_exp      !< Maximum possible exponent
   real(kind=dp) :: this_exp     !< Exponent argument in this step
   real(kind=dp) :: this_hermite !< Hermite argument in this step
   integer :: i_state  !< State index
   integer :: i_kpoint !< K-point index
   integer :: i_spin   !< Spin index

   character*40, parameter :: caller = "elsi_check_electrons"

   if(broaden_width .le. 0.0_dp) then
      call elsi_stop(" Broadening width in chemical potential determination must"//&
                     " be a positive number. Exiting...",caller)
   endif

   invert_width = 1.0_dp/broaden_width
   diff_ne_out = 0.0_dp

   if((spin_degen /= 1.0_dp) .and. (spin_degen /= 2.0_dp)) then ! Not set by user
      if(n_spin == 2) then
         spin_degen = 1.0_dp
      else
         spin_degen = 2.0_dp
      endif
   endif

   select case (broaden_method)
      case(GAUSSIAN)
         do i_kpoint = 1,n_kpoint
            do i_spin = 1,n_spin
               do i_state = 1,n_state
                  occ_numbers(i_state,i_spin,i_kpoint) = spin_degen*0.5_dp*&
                     (1.0_dp-ERF((eigenvalues(i_state,i_spin,i_kpoint)-mu_in)*invert_width))

                  diff_ne_out = diff_ne_out+&
                     occ_numbers(i_state,i_spin,i_kpoint)*kpoint_weights(i_kpoint)
               enddo
            enddo
         enddo

      case(FERMI)
         max_exp = MAXEXPONENT(mu_in)*LOG(2.0_dp)

         do i_kpoint = 1,n_kpoint
            do i_spin = 1,n_spin
               do i_state = 1,n_state
                  this_exp = (eigenvalues(i_state,i_spin,i_kpoint)-mu_in)*invert_width

                  if(this_exp < max_exp) then
                     occ_numbers(i_state,i_spin,i_kpoint) = spin_degen/(1.0_dp+EXP(this_exp))

                     diff_ne_out = diff_ne_out+&
                        occ_numbers(i_state,i_spin,i_kpoint)*kpoint_weights(i_kpoint)
                  else ! Exponent in this step is larger than the largest possible exponent
                     occ_numbers(i_state,i_spin,i_kpoint) = 0.0_dp
                  endif
               enddo
            enddo
         enddo

      case(METHFESSEL_PAXTON_0)
         do i_kpoint = 1,n_kpoint
            do i_spin = 1,n_spin
               do i_state = 1,n_state
                  occ_numbers(i_state,i_spin,i_kpoint) = spin_degen*0.5_dp*&
                     (1.0_dp-ERF((eigenvalues(i_state,i_spin,i_kpoint)-mu_in)*invert_width))

                  diff_ne_out = diff_ne_out+&
                     occ_numbers(i_state,i_spin,i_kpoint)*kpoint_weights(i_kpoint)
               enddo
            enddo
         enddo

      case(METHFESSEL_PAXTON_1)
         do i_kpoint = 1,n_kpoint
            do i_spin = 1,n_spin
               do i_state = 1,n_state
                  this_hermite = (eigenvalues(i_state,i_spin,i_kpoint)-mu_in)*invert_width

                  occ_numbers(i_state,i_spin,i_kpoint) = spin_degen*0.5_dp*&
                     (1.0_dp-ERF(this_hermite))-0.5_dp*invert_sqrt_pi*this_hermite*&
                     EXP(-this_hermite*this_hermite)

                  diff_ne_out = diff_ne_out+&
                     occ_numbers(i_state,i_spin,i_kpoint)*kpoint_weights(i_kpoint)
               enddo
            enddo
         enddo

      case DEFAULT
         call elsi_stop(" No supperted broadening scheme has been chosen."//&
                        " Please choose GAUSSIAN, FERMI, METHFESSEL_PAXTON_0,"//&
                        " or METHFESSEL_PAXTON_1 broadening scheme."//&
                        " Exiting...",caller)
   end select

   diff_ne_out = diff_ne_out-n_electron

end subroutine

!>
!! This routine computes the chemical potential using bisection algorithm.
!!
subroutine elsi_find_mu(kpoint_weights,eigenvalues,occ_numbers,&
                        mu_lower_in,mu_upper_in,mu_out)

   implicit none

   real(kind=dp), intent(in)  :: kpoint_weights(n_kpoint)
   real(kind=dp), intent(in)  :: eigenvalues(n_state,n_spin,n_kpoint)
   real(kind=dp), intent(out) :: occ_numbers(n_state,n_spin,n_kpoint)
   real(kind=dp), intent(in)  :: mu_lower_in !< Lower bound of chemical potential
   real(kind=dp), intent(in)  :: mu_upper_in !< Upper bound of chemical potential
   real(kind=dp), intent(out) :: mu_out      !< Solution of chemical potential

   real(kind=dp) :: mu_left    !< Left bound of chemical potential interval
   real(kind=dp) :: mu_right   !< Right bound of chemical potential interval
   real(kind=dp) :: mu_mid     !< Middle point of chemical potential interval
   real(kind=dp) :: diff_left  !< Difference in number of electrons on left bound
   real(kind=dp) :: diff_right !< Difference in number of electrons on right bound
   real(kind=dp) :: diff_mid   !< Difference in number of electrons on middle point
   logical :: found_mu !< Is chemical potential found?
   integer :: n_steps  !< Number of steps to find chemical potential

   character*200 :: info_str
   character*40, parameter :: caller = "elsi_find_mu"

   n_steps = 0
   found_mu = .false.

   mu_left = mu_lower_in
   mu_right = mu_upper_in

   do while((.not.found_mu) .and. (n_steps < max_mu_steps))
      call elsi_check_electrons(kpoint_weights,eigenvalues,occ_numbers,&
                                mu_left,diff_left)
      call elsi_check_electrons(kpoint_weights,eigenvalues,occ_numbers,&
                                mu_right,diff_right)

      if(ABS(diff_left) < occ_tolerance) then
         mu_out = mu_left
         found_mu = .true.
      elseif(ABS(diff_right) < occ_tolerance) then
         mu_out = mu_right
         found_mu = .true.
      else
         n_steps = n_steps+1

         mu_mid = 0.5_dp*(mu_left+mu_right)

         call elsi_check_electrons(kpoint_weights,eigenvalues,occ_numbers,&
                                   mu_mid,diff_mid)

         if(ABS(diff_mid) < occ_tolerance) then
            mu_out = mu_mid
            found_mu = .true.
         elseif(diff_mid < 0) then
            mu_left = mu_mid
         elseif(diff_mid > 0) then
            mu_right = mu_mid
         endif
      endif
   enddo

   ! Special treatment if mu does not reach the required accuracy
   if(.not.found_mu) then
      ! Use the chemical potential of the right bound
      call elsi_check_electrons(kpoint_weights,eigenvalues,occ_numbers,&
                                mu_right,diff_right)

      mu_out = mu_right

      ! With adjusted occupation numbers
      call elsi_statement_print("  Chemical potential cannot reach the"//&
                                " required accuracy by bisection method."//&
                                " The error will be arbitrarily canceled"//&
                                " in the following state(s):")

      call elsi_adjust_occ(kpoint_weights,eigenvalues,occ_numbers,diff_right)
   endif

end subroutine

!>
!! This routine cancels the small error in number of electrons.
!!
subroutine elsi_adjust_occ(kpoint_weights,eigenvalues,occ_numbers,diff_ne)

   implicit none

   real(kind=dp), intent(in)    :: kpoint_weights(n_kpoint)
   real(kind=dp), intent(in)    :: eigenvalues(n_state,n_spin,n_kpoint)
   real(kind=dp), intent(inout) :: occ_numbers(n_state,n_spin,n_kpoint)
   real(kind=dp), intent(inout) :: diff_ne !< Error in number of electrons before adjusting

   real(kind=dp) :: tmp_real
   integer       :: tmp_int
   integer       :: min_id
   integer       :: i_state  !< State index
   integer       :: i_kpoint !< K-point index
   integer       :: i_spin   !< Spin index
   integer       :: i_val    !< Index for 1D array
   integer       :: n_total  !< Total number of states

   real(kind=dp), allocatable :: eval_aux(:) !< Aux 1D array
   real(kind=dp), allocatable :: occ_aux(:)  !< Aux 1D array
   integer,       allocatable :: idx_aux(:)  !< Aux 1D array

   character*200 :: info_str

   character*40, parameter :: caller = "elsi_adjust_occ"

   n_total = n_state*n_spin*n_kpoint

   call elsi_allocate(eval_aux,n_total,"eval_aux",caller)
   call elsi_allocate(occ_aux,n_total,"occ_aux",caller)
   call elsi_allocate(idx_aux,n_total,"idx_aux",caller)

   ! Put eigenvalues and occupation numbers into 1D arrays
   i_val = 0

   do i_kpoint = 1,n_kpoint
      do i_spin = 1,n_spin
         do i_state = 1,n_state
            i_val = i_val+1
            eval_aux(i_val) = eigenvalues(i_state,i_spin,i_kpoint)
            occ_aux(i_val)  = occ_numbers(i_state,i_spin,i_kpoint)
            idx_aux(i_val)  = i_val
         enddo
      enddo
   enddo

   ! Sort eval_aux; move occ_aux and idx_aux accordingly
   do i_val = 1,n_total
      min_id = MINLOC(eval_aux(i_val:n_total),1)+i_val-1

      tmp_int = idx_aux(i_val)
      idx_aux(i_val) = idx_aux(min_id)
      idx_aux(min_id) = tmp_int

      tmp_real = eval_aux(i_val)
      eval_aux(i_val) = eval_aux(min_id)
      eval_aux(min_id) = tmp_real

      tmp_real = occ_aux(i_val)
      occ_aux(i_val) = occ_aux(min_id)
      occ_aux(min_id) = tmp_real
   enddo

   deallocate(eval_aux)

   do i_val = n_total,1,-1
      if(occ_aux(i_val) > 0.0_dp) then
         i_kpoint = (idx_aux(i_val)-1)/(n_spin*n_state)+1
         i_spin   = MOD((idx_aux(i_val)-1)/6,2)+1
         i_state  = MOD(idx_aux(i_val)-1,n_state)+1

         write(info_str,"(A,I5,A,I5,A,I5)") "  | k-point ",i_kpoint,&
               ", spin channel ",i_spin,", state ",i_state
         call elsi_statement_print(info_str)

         if(kpoint_weights(i_kpoint)*occ_aux(i_val) > diff_ne) then
            occ_aux(i_val) = occ_aux(i_val)-diff_ne/kpoint_weights(i_kpoint)
            diff_ne = 0.0_dp
         else
            diff_ne = diff_ne-kpoint_weights(i_kpoint)*occ_aux(i_val)
            occ_aux(i_val) = 0.0_dp
         endif

         occ_numbers(i_state,i_spin,i_kpoint) = occ_aux(i_val)
      endif

      if((diff_ne < occ_tolerance) .or. (diff_ne == 0.0_dp)) exit
   enddo

   deallocate(occ_aux)
   deallocate(idx_aux)

end subroutine

end module ELSI_MU
