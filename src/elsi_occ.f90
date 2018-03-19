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
   use ELSI_DATATYPE,  only: elsi_handle
   use ELSI_IO,        only: elsi_say
   use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
   use ELSI_MPI,       only: elsi_stop
   use ELSI_PRECISION, only: r8,i4

   implicit none

   private

   public :: elsi_compute_mu_and_occ
   public :: elsi_compute_entropy

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
   real(kind=r8)    :: diff_ne_lower ! Electron count error on lower bound
   real(kind=r8)    :: diff_ne_upper ! Electron count error on upper bound
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: n_steps

   character(len=40), parameter :: caller = "elsi_compute_mu_and_occ"

   ! Determine smallest and largest eivenvalues
   e_low  = evals(1,1,1)
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

   real(kind=r8)    :: invert_width ! 1/broaden_width
   real(kind=r8)    :: delta
   real(kind=r8)    :: max_exp ! Maximum possible exponent
   real(kind=r8)    :: arg
   real(kind=r8)    :: weight
   real(kind=r8)    :: A
   real(kind=r8)    :: H_even
   real(kind=r8)    :: H_odd
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: i_mp

   character(len=40), parameter :: caller = "elsi_check_electrons"

   invert_width = 1.0_r8/e_h%broaden_width
   diff_ne_out  = 0.0_r8

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
               arg = (evals(i_state,i_spin,i_kpt)-mu_in)*invert_width

               if(arg < max_exp) then
                  occ_nums(i_state,i_spin,i_kpt) = e_h%spin_degen/&
                                                      (1.0_r8+exp(arg))

                  diff_ne_out = diff_ne_out+&
                     occ_nums(i_state,i_spin,i_kpt)*k_weights(i_kpt)
               else ! Exponent larger than maximum value
                  occ_nums(i_state,i_spin,i_kpt) = 0.0_r8
               endif
            enddo
         enddo
      enddo
   case(METHFESSEL_PAXTON)
      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg    = (evals(i_state,i_spin,i_kpt)-mu_in)*invert_width
               weight = exp(-arg*arg)

               occ_nums(i_state,i_spin,i_kpt) = &
                  0.5_r8*(1.0_r8-erf(arg))*e_h%spin_degen

               if(e_h%mp_order > 0) then ! 1st order
                  A      = -0.25_r8*INVERT_SQRT_PI
                  H_even = 1.0_r8
                  H_odd  = 2.0_r8*arg

                  occ_nums(i_state,i_spin,i_kpt) = &
                     occ_nums(i_state,i_spin,i_kpt)+&
                     A*H_odd*weight*e_h%spin_degen
               endif

               if(e_h%mp_order > 1) then ! higher order
                  do i_mp = 2,e_h%mp_order
                     A      = -1.0_r8/real(4*i_mp,kind=r8)*A
                     H_even = 2.0_r8*arg*H_odd-2.0_r8*i_mp*H_even
                     H_odd  = 2.0_r8*arg*H_even-2.0_r8*(i_mp+1)*H_odd

                     occ_nums(i_state,i_spin,i_kpt) = &
                        occ_nums(i_state,i_spin,i_kpt)+&
                        A*H_odd*weight*e_h%spin_degen
                  enddo
               endif

               diff_ne_out = diff_ne_out+occ_nums(i_state,i_spin,i_kpt)*&
                                k_weights(i_kpt)
            enddo
         enddo
      enddo
   case(CUBIC)
      ! To have a consistent slope of the occupation function at the chemical
      ! potential, the parameters for GAUSSIAN and CUBIC should be related as:
      delta = 0.75_r8*SQRT_PI*e_h%broaden_width

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               if(evals(i_state,i_spin,i_kpt) <= mu_in-delta) then
                  occ_nums(i_state,i_spin,i_kpt) = e_h%spin_degen*1.0_r8
               elseif(evals(i_state,i_spin,i_kpt) >= mu_in+delta) then
                  occ_nums(i_state,i_spin,i_kpt) = 0.0_r8
               else
                  occ_nums(i_state,i_spin,i_kpt) = e_h%spin_degen*&
                     (0.25_r8/delta**3)*&
                     (evals(i_state,i_spin,i_kpt)-mu_in+2*delta)*&
                     (evals(i_state,i_spin,i_kpt)-mu_in-delta)**2
               endif

               diff_ne_out = diff_ne_out+occ_nums(i_state,i_spin,i_kpt)*&
                                k_weights(i_kpt)
            enddo
         enddo
      enddo
   case(COLD)
      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg = (evals(i_state,i_spin,i_kpt)-mu_in)*invert_width
               arg = -arg-sqrt(0.5_r8)

               occ_nums(i_state,i_spin,i_kpt) = (0.5_r8+erf(arg)*0.5_r8+&
                  INVERT_SQRT_PI*sqrt(0.5_r8)*exp(-arg**2))*e_h%spin_degen

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

   real(kind=r8)      :: mu_left
   real(kind=r8)      :: mu_right
   real(kind=r8)      :: mu_mid
   real(kind=r8)      :: diff_left  ! Electron count error on left bound
   real(kind=r8)      :: diff_right ! Electron count error on right bound
   real(kind=r8)      :: diff_mid   ! Electron count error on middle point
   logical            :: found_mu
   integer(kind=i4)   :: n_steps
   character(len=200) :: info_str

   character(len=40), parameter :: caller = "elsi_find_mu"

   n_steps  = 0
   found_mu = .false.
   mu_left  = mu_lower_in
   mu_right = mu_upper_in

   call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
           evals,occ_nums,mu_left,diff_left)
   call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
           evals,occ_nums,mu_right,diff_right)

   if(abs(diff_left) < e_h%occ_tolerance) then
      mu_out   = mu_left
      found_mu = .true.
   elseif(abs(diff_right) < e_h%occ_tolerance) then
      mu_out   = mu_right
      found_mu = .true.
   endif

   do while(.not. found_mu .and. n_steps < e_h%max_mu_steps)
      n_steps = n_steps+1
      mu_mid  = 0.5_r8*(mu_left+mu_right)

      call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_mid,diff_mid)

      if(abs(diff_mid) < e_h%occ_tolerance) then
         mu_out   = mu_mid
         found_mu = .true.
      elseif(diff_mid < 0) then
         mu_left = mu_mid
      elseif(diff_mid > 0) then
         mu_right = mu_mid
      endif
   enddo

   if(found_mu) then
      call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_out,diff_right)
   else ! mu cannot reach required accuracy
      ! Use mu of the right bound...
      call elsi_check_electrons(e_h,n_electron,n_state,n_spin,n_kpt,k_weights,&
              evals,occ_nums,mu_right,diff_right)

      mu_out = mu_right

      ! ...with adjusted occupation numbers
      call elsi_say(e_h,"  Chemical potential cannot reach the required"//&
              " accuracy by bisection method.")
      write(info_str,"('  | Residual error :',E10.2)") diff_right
      call elsi_say(e_h,info_str)
      call elsi_say(e_h,"  The error will be arbitrarily removed from the"//&
              " highest occupied states.")

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

   character(len=40), parameter :: caller = "elsi_adjust_occ"

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

      i_kpt   = (max_id-1)/(n_spin*n_state)+1
      i_spin  = mod((max_id-1)/n_state,n_spin)+1
      i_state = mod(max_id-1,n_state)+1

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

!>
!! This routine computes the entropy.
!!
subroutine elsi_compute_entropy(e_h,n_state,n_spin,n_kpt,k_weights,evals,&
              occ_nums,mu,entropy)

   implicit none

   type(elsi_handle), intent(inout) :: e_h                            !< Handle
   integer(kind=i4),  intent(in)    :: n_state                        !< Number of states
   integer(kind=i4),  intent(in)    :: n_spin                         !< Number of spins
   integer(kind=i4),  intent(in)    :: n_kpt                          !< Number of k-points
   real(kind=r8),     intent(in)    :: k_weights(n_kpt)               !< K-points weights
   real(kind=r8),     intent(in)    :: evals(n_state,n_spin,n_kpt)    !< Eigenvalues
   real(kind=r8),     intent(in)    :: occ_nums(n_state,n_spin,n_kpt) !< Occupation numbers
   real(kind=r8),     intent(in)    :: mu                             !< Input chemical potential
   real(kind=r8),     intent(out)   :: entropy                        !< Entropy

   real(kind=r8)    :: invert_width ! 1/broaden_width
   real(kind=r8)    :: delta
   real(kind=r8)    :: const
   real(kind=r8)    :: pre
   real(kind=r8)    :: arg
   real(kind=r8)    :: weight
   real(kind=r8)    :: A
   real(kind=r8)    :: H_even
   real(kind=r8)    :: H_odd
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: i_mp

   real(kind=r8),     parameter :: entropy_thr = 1.0e-15_r8
   character(len=40), parameter :: caller = "elsi_compute_entropy"

   invert_width = 1.0_r8/e_h%broaden_width
   entropy      = 0.0_r8

   if(.not. e_h%spin_is_set) then
      if(n_spin == 2) then
         e_h%spin_degen = 1.0_r8
      else
         e_h%spin_degen = 2.0_r8
      endif
   endif

   select case(e_h%broaden_scheme)
   case(GAUSSIAN)
      pre = e_h%spin_degen*0.25_r8*e_h%broaden_width*INVERT_SQRT_PI

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg     = (evals(i_state,i_spin,i_kpt)-mu)*invert_width
               arg     = -(arg*arg)
               entropy = entropy+exp(arg)*k_weights(i_kpt)
            enddo
         enddo
      enddo
   case(FERMI)
      pre = e_h%spin_degen*0.5_r8*e_h%broaden_width

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg = occ_nums(i_state,i_spin,i_kpt)/e_h%spin_degen

               if(1-arg > entropy_thr .and. arg > entropy_thr) then
                  entropy = entropy+(arg*log(arg)+(1-arg)*log(1-arg))*&
                               k_weights(i_kpt)
               endif
            enddo
         enddo
      enddo
   case(METHFESSEL_PAXTON)
      pre = e_h%spin_degen*0.25_r8*e_h%broaden_width

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg     = (evals(i_state,i_spin,i_kpt)-mu)*invert_width
               weight  = exp(-arg*arg)
               A       = INVERT_SQRT_PI
               H_even  = 1.0_r8
               H_odd   = 2.0_r8*arg
               entropy = entropy+INVERT_SQRT_PI*weight*k_weights(i_kpt)

               do i_mp = 1,e_h%mp_order
                  A      = -1.0_r8/real(4*i_mp,kind=i4)*A
                  H_even = 2.0_r8*arg*H_odd-2.0_r8*i_mp*H_even
                  H_odd  = 2.0_r8*arg*H_even-2.0_r8*(i_mp+1)*H_odd

                  entropy = entropy+A*H_even*weight*k_weights(i_kpt)
               enddo
            enddo
         enddo
      enddo
   case(CUBIC)
      delta = 0.75_r8*SQRT_PI*e_h%broaden_width
      pre   = e_h%spin_degen*0.5_r8*e_h%broaden_width*0.1875_r8/(delta**4)
      const = delta**2

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               if(evals(i_state,i_spin,i_kpt) > mu-delta .and. &
                  evals(i_state,i_spin,i_kpt) < mu+delta) then
                  arg     = evals(i_state,i_spin,i_kpt)-mu
                  entropy = entropy+k_weights(i_kpt)*(((arg**2)-const)**2)
               endif
            enddo
         enddo
      enddo
   case(COLD)
      pre = e_h%spin_degen*INVERT_SQRT_PI*sqrt(0.5_r8)*e_h%broaden_width

      do i_kpt = 1,n_kpt
         do i_spin = 1,n_spin
            do i_state = 1,n_state
               arg = (evals(i_state,i_spin,i_kpt)-mu)*invert_width
               arg = -arg-sqrt(0.5_r8)

               entropy = entropy-arg*exp(-arg**2)
            enddo
         enddo
      enddo
   end select

   entropy = -(pre*entropy)

end subroutine

end module ELSI_OCC
