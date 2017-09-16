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

   use ELSI_CONSTANTS, only: GAUSSIAN,FERMI,METHFESSEL_PAXTON_0,&
                             METHFESSEL_PAXTON_1,INVERT_SQRT_PI
   use ELSI_DATATYPE
   use ELSI_MALLOC
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS

   implicit none

   private

   public :: elsi_compute_mu_and_occ

contains

!>
!! This routine computes the chemical potential and occupation numbers.
!!
subroutine elsi_compute_mu_and_occ(e_h,n_electron,n_state,n_spin,n_kpt,&
              k_weights,evals,occ_nums,mu)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                            !< Handle
   real(kind=r8),     intent(in)    :: n_electron                     !< Number of electrons
   integer(kind=i4),  intent(in)    :: n_state                        !< Number of states
   integer(kind=i4),  intent(in)    :: n_spin                         !< Number of spins
   integer(kind=i4),  intent(in)    :: n_kpt                          !< Number of k-points
   real(kind=r8),     intent(in)    :: k_weights(n_kpt)               !< K-points weights
   real(kind=r8),     intent(in)    :: evals(n_state,n_spin,n_kpt)    !< Eigenvalues
   real(kind=r8),     intent(out)   :: occ_nums(n_state,n_spin,n_kpt) !< Occupation members
   real(kind=r8),     intent(out)   :: mu                             !< Chemical potential

   real(kind=r8)    :: e_low         ! Lowest eigenvalue
   real(kind=r8)    :: e_high        ! Highest eigenvalue
   real(kind=r8)    :: mu_lower      ! Lower bound of chemical potential
   real(kind=r8)    :: mu_upper      ! Upper bound of chemical potential
   real(kind=r8)    :: diff_ne_lower ! Difference in number of electrons on lower bound
   real(kind=r8)    :: diff_ne_upper ! Difference in number of electrons on upper bound
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: n_steps

   character*40, parameter :: caller = "elsi_compute_mu_and_occ"

   ! Determine smallest and largest eivenvalues
   e_low = evals(1,1,1)
   e_high = evals(n_state,1,1)

   do i_kpt = 1,n_kpt
      do i_spin = 1,n_spin
         do i_state = 1,n_state
            if(evals(i_state,i_spin,i_kpt) < e_low) then
               e_low = evals(i_state,i_spin,i_kpt)
            endif
            if(evals(i_state,i_spin,i_kpt) > e_high) then
               e_high = evals(i_state,i_spin,i_kpt)
            endif
         enddo
      enddo
   enddo

   ! Determine upper and lower bounds of mu
   mu_lower = e_low

   if(e_low == e_high) then
      mu_upper = 0.0_r8
   else
      mu_upper = e_high
   endif

   occ_nums = 0.0_r8

   ! Compute electron count error
   call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
           evals,occ_nums,mu_lower,diff_ne_lower)
   call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
           evals,occ_nums,mu_upper,diff_ne_upper)

   ! Enlarge the interval towards both sides if solution not found
   n_steps = 0
   do while(diff_ne_lower*diff_ne_upper > 0)
      n_steps = n_steps+1

      if(n_steps > e_h%max_mu_steps) then
         call elsi_stop(" Chemical potential not found.",e_h,caller)
      endif

      mu_lower = mu_lower-0.5_r8*abs(e_high-e_low)
      mu_upper = mu_upper+0.5_r8*abs(e_high-e_low)

      call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_lower,diff_ne_lower)
      call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_upper,diff_ne_upper)
   enddo

   ! Perform bisection
   call elsi_find_mu(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,evals,&
           occ_nums,mu_lower,mu_upper,mu)

end subroutine

!>
!! This routine computes the number of electrons for a given chemical potential,
!! and returns the error in the number of electrons. The occupation numbers are
!! updated as well.
!!
subroutine elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_in,diff_ne_out)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                            !< Handle
   real(kind=r8),     intent(in)    :: n_electron                     !< Number of electrons
   integer(kind=i4),  intent(in)    :: n_state                        !< Number of states
   integer(kind=i4),  intent(in)    :: n_spin                         !< Number of spins
   integer(kind=i4),  intent(in)    :: n_kpt                          !< Number of k-points
   real(kind=r8),     intent(in)    :: k_weights(n_kpt)               !< K-points weights
   real(kind=r8),     intent(in)    :: evals(n_state,n_spin,n_kpt)    !< Eigenvalues
   real(kind=r8),     intent(out)   :: occ_nums(n_state,n_spin,n_kpt) !< Occupation numbers
   real(kind=r8),     intent(in)    :: mu_in                          !< Input chemical potential
   real(kind=r8),     intent(out)   :: diff_ne_out                    !< Electron count error

   real(kind=r8) :: invert_width ! 1/broaden_width
   real(kind=r8) :: max_exp ! Maximum possible exponent
   real(kind=r8) :: this_exp
   real(kind=r8) :: this_hermite

   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin

   character*40, parameter :: caller = "elsi_check_electrons"

   invert_width = 1.0_r8/e_h%broaden_width
   diff_ne_out = 0.0_r8

   if(.not. e_h%spin_is_set) then
      if(n_spin == 2) then
         e_h%spin_degen = 1.0_r8
      else
         e_h%spin_degen = 2.0_r8
      endif
   endif

   select case(e_h%broaden_scheme)
   case(GAUSSIAN)
      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               occ_nums(i_state,i_spin,i_kpt) = e_h%spin_degen*0.5_r8*&
                  (1.0_r8-erf((evals(i_state,i_spin,i_kpt)-mu_in)*invert_width))

               diff_ne_out = diff_ne_out+&
                  occ_nums(i_state,i_spin,i_kpt)*k_weights(i_kpt)
            enddo
         enddo
      enddo
   case(FERMI)
      max_exp = maxexponent(mu_in)*log(2.0_r8)

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               this_exp = (evals(i_state,i_spin,i_kpt)-mu_in)*invert_width

               if(this_exp < max_exp) then
                  occ_nums(i_state,i_spin,i_kpt) = e_h%spin_degen/&
                                                      (1.0_r8+exp(this_exp))

                  diff_ne_out = diff_ne_out+&
                     occ_nums(i_state,i_spin,i_kpt)*k_weights(i_kpt)
               else ! Exponent larger than maximum value
                  occ_nums(i_state,i_spin,i_kpt) = 0.0_r8
               endif
            enddo
         enddo
      enddo
   case(METHFESSEL_PAXTON_0)
      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               occ_nums(i_state,i_spin,i_kpt) = e_h%spin_degen*0.5_r8*&
                  (1.0_r8-erf((evals(i_state,i_spin,i_kpt)-mu_in)*invert_width))

               diff_ne_out = diff_ne_out+occ_nums(i_state,i_spin,i_kpt)*&
                                k_weights(i_kpt)
            enddo
         enddo
      enddo
   case(METHFESSEL_PAXTON_1)
      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               this_hermite = (evals(i_state,i_spin,i_kpt)-mu_in)*invert_width

               occ_nums(i_state,i_spin,i_kpt) = e_h%spin_degen*0.5_r8*&
                  (1.0_r8-erf(this_hermite))-0.5_r8*INVERT_SQRT_PI*&
                  this_hermite*exp(-this_hermite*this_hermite)

               diff_ne_out = diff_ne_out+occ_nums(i_state,i_spin,i_kpt)*&
                                k_weights(i_kpt)
            enddo
         enddo
      enddo
   end select

   diff_ne_out = diff_ne_out-n_electron

end subroutine

!>
!! This routine computes the chemical potential using a bisection algorithm.
!!
subroutine elsi_find_mu(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,evals,&
              occ_nums,mu_lower_in,mu_upper_in,mu_out)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                            !< Handle
   real(kind=r8),     intent(in)    :: n_electron                     !< Number of electrons
   integer(kind=i4),  intent(in)    :: n_state                        !< Number of states
   integer(kind=i4),  intent(in)    :: n_spin                         !< Number of spins
   integer(kind=i4),  intent(in)    :: n_kpt                          !< Number of k-points
   real(kind=r8),     intent(in)    :: k_weights(n_kpt)               !< K-points weights
   real(kind=r8),     intent(in)    :: evals(n_state,n_spin,n_kpt)    !< Eigenvalues
   real(kind=r8),     intent(out)   :: occ_nums(n_state,n_spin,n_kpt) !< Occupation numbers
   real(kind=r8),     intent(in)    :: mu_lower_in                    !< Lower bound of mu
   real(kind=r8),     intent(in)    :: mu_upper_in                    !< Upper bound of mu
   real(kind=r8),     intent(out)   :: mu_out                         !< Solution

   real(kind=r8)    :: mu_left    ! Left bound of chemical potential interval
   real(kind=r8)    :: mu_right   ! Right bound of chemical potential interval
   real(kind=r8)    :: mu_mid     ! Middle point of chemical potential interval
   real(kind=r8)    :: diff_left  ! Difference in number of electrons on left bound
   real(kind=r8)    :: diff_right ! Difference in number of electrons on right bound
   real(kind=r8)    :: diff_mid   ! Difference in number of electrons on middle point
   logical          :: found_mu
   integer(kind=i4) :: n_steps

   character*40, parameter :: caller = "elsi_find_mu"

   n_steps = 0
   found_mu = .false.

   mu_left = mu_lower_in
   mu_right = mu_upper_in

   do while(.not. found_mu .and. n_steps < e_h%max_mu_steps)
      call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_left,diff_left)
      call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_right,diff_right)

      if(abs(diff_left) < e_h%occ_tolerance) then
         mu_out = mu_left
         found_mu = .true.
      elseif(abs(diff_right) < e_h%occ_tolerance) then
         mu_out = mu_right
         found_mu = .true.
      else
         n_steps = n_steps+1

         mu_mid = 0.5_r8*(mu_left+mu_right)

         call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,&
                 k_weights,evals,occ_nums,mu_mid,diff_mid)

         if(abs(diff_mid) < e_h%occ_tolerance) then
            mu_out = mu_mid
            found_mu = .true.
         elseif(diff_mid < 0) then
            mu_left = mu_mid
         elseif(diff_mid > 0) then
            mu_right = mu_mid
         endif
      endif
   enddo

   ! Special treatment if mu cannot reach the required accuracy
   if(.not. found_mu) then
      ! Use the chemical potential of the right bound...
      call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_right,diff_right)

      mu_out = mu_right

      ! ...with adjusted occupation numbers
      call elsi_say("  Chemical potential cannot reach the"//&
              " required accuracy by bisection method.",e_h)
      call elsi_say("  The error will be arbitrarily removed"//&
              " from the highest occupied states.",e_h)

      call elsi_adjust_occ(e_h,n_state,n_spin,n_kpt,k_weights,evals,occ_nums,&
              diff_right)
   endif

end subroutine

!>
!! This routine cancels the small error in number of electrons.
!!
subroutine elsi_adjust_occ(e_h,n_state,n_spin,n_kpt,k_weights,evals,occ_nums,&
              diff_ne)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                            !< Handle
   integer(kind=i4),  intent(in)    :: n_state                        !< Number of states
   integer(kind=i4),  intent(in)    :: n_spin                         !< Number of spins
   integer(kind=i4),  intent(in)    :: n_kpt                          !< Number of k-points
   real(kind=r8),     intent(in)    :: k_weights(n_kpt)               !< K-points weights
   real(kind=r8),     intent(in)    :: evals(n_state,n_spin,n_kpt)    !< Eigenvalues
   real(kind=r8),     intent(inout) :: occ_nums(n_state,n_spin,n_kpt) !< Occupation numbers
   real(kind=r8),     intent(inout) :: diff_ne                        !< Electron count error

   real(kind=r8), allocatable :: eval_aux(:)

   real(kind=r8)    :: min_eval
   integer(kind=i4) :: max_id
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: i_val
   integer(kind=i4) :: n_total

   character*40, parameter :: caller = "elsi_adjust_occ"

   n_total = n_state*n_spin*n_kpt

   call elsi_allocate(e_h,eval_aux,n_total,"eval_aux",caller)

   ! Put evals into a 1D array
   i_val = 0

   do i_kpt = 1,n_kpt
      do i_spin = 1,n_spin
         do i_state = 1,n_state
            i_val = i_val+1
            eval_aux(i_val) = evals(i_state,i_spin,i_kpt)
         enddo
      enddo
   enddo

   min_eval = minval(eval_aux,1)

   ! Remove error
   do i_val = 1,n_total
      max_id           = maxloc(eval_aux,1)
      eval_aux(max_id) = min_eval-1.0_r8

      i_kpt = (i_val-1)/(n_spin*n_state)+1
      i_spin   = mod((i_val-1)/n_state,n_spin)+1
      i_state  = mod(i_val-1,n_state)+1

      if(k_weights(i_kpt)*occ_nums(i_state,i_spin,i_kpt) > diff_ne) then
         occ_nums(i_state,i_spin,i_kpt) = occ_nums(i_state,i_spin,i_kpt)-&
                                             diff_ne/k_weights(i_kpt)
         diff_ne = 0.0_r8
      else
         diff_ne = diff_ne-k_weights(i_kpt)*occ_nums(i_state,i_spin,i_kpt)
         occ_nums(i_state,i_spin,i_kpt) = 0.0_r8
      endif

      if(diff_ne <= e_h%occ_tolerance) exit
   enddo

   call elsi_deallocate(e_h,eval_aux,"eval_aux")

end subroutine

end module ELSI_MU
