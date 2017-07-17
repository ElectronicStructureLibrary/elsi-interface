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

   use ISO_C_BINDING
   use ELSI_CONSTANTS, only: GAUSSIAN,FERMI,METHFESSEL_PAXTON_0,&
                             METHFESSEL_PAXTON_1,INVERT_SQRT_PI
   use ELSI_DIMENSIONS, only: elsi_handle
   use ELSI_PRECISION, only: r8,i4
   use ELSI_UTILS

   implicit none
   private

   public :: elsi_compute_mu_and_occ

contains

!>
!! This routine computes the chemical potential and occupation numbers.
!!
subroutine elsi_compute_mu_and_occ(elsi_h)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h !< Handle

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
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_compute_mu_and_occ"

   call elsi_get_eval_all(elsi_h)

   ! Determine the smallest and largest eivenvalues
   e_low = elsi_h%eval_all(1,1,1)
   e_high = elsi_h%eval_all(elsi_h%n_states,1,1)

   do i_kpt = 1,elsi_h%n_kpts
      do i_spin = 1,elsi_h%n_spins
         do i_state = 1,elsi_h%n_states
            if(elsi_h%eval_all(i_state,i_spin,i_kpt) < e_low) then
               e_low = elsi_h%eval_all(i_state,i_spin,i_kpt)
            endif
            if(elsi_h%eval_all(i_state,i_spin,i_kpt) > e_high) then
               e_high = elsi_h%eval_all(i_state,i_spin,i_kpt)
            endif
         enddo
      enddo
   enddo

   ! Determine the upper and lower bounds for the chemical potential search
   mu_lower = e_low

   if(e_low == e_high) then
      mu_upper = 0.0_r8
   else
      mu_upper = e_high
   endif

   elsi_h%occ_num = 0.0_r8

   ! Compute the error in electron count
   call elsi_check_electrons(elsi_h,mu_lower,diff_ne_lower)
   call elsi_check_electrons(elsi_h,mu_upper,diff_ne_upper)

   ! If diff_ne_lower*diff_ne_upper > 0, the solution is not in this interval.
   ! Enlarge the interval towards both sides, then recheck.
   n_steps = 0
   do while(diff_ne_lower*diff_ne_upper > 0)
      n_steps = n_steps+1
      if(n_steps > elsi_h%max_mu_steps) then
         write(info_str,"(A,I13,A)") " Chemical potential not found in ",&
            elsi_h%max_mu_steps," iterations! Exiting..."
         call elsi_stop(info_str,elsi_h,caller)
      endif

      mu_lower = mu_lower-0.5_r8*abs(e_high-e_low)
      mu_upper = mu_upper+0.5_r8*abs(e_high-e_low)

      call elsi_check_electrons(elsi_h,mu_lower,diff_ne_lower)
      call elsi_check_electrons(elsi_h,mu_upper,diff_ne_upper)
   enddo

   ! Now the solution should lie in the interval. Use bisection to find it.
   call elsi_find_mu(elsi_h,mu_lower,mu_upper)

end subroutine

!>
!! This routine computes the number of electrons for a given chemical potential,
!! and returns the error in the number of electrons. The occupation numbers will
!! be updated as well.
!!
subroutine elsi_check_electrons(elsi_h,mu_in,diff_ne_out)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h      !< Handle
   real(kind=r8),     intent(in)    :: mu_in       !< Input chemical potential
   real(kind=r8),     intent(out)   :: diff_ne_out !< Difference in number of electrons

   real(kind=r8) :: invert_width ! 1/broaden_width
   real(kind=r8) :: max_exp ! Maximum possible exponent
   real(kind=r8) :: this_exp
   real(kind=r8) :: this_hermite

   integer(kind=i4) :: n_steps
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin

   character*40, parameter :: caller = "elsi_check_electrons"

   if(elsi_h%broadening_width .le. 0.0_r8) then
      call elsi_stop(" Broadening width in chemical potential determination must"//&
              " be a positive number. Exiting...",elsi_h,caller)
   endif

   invert_width = 1.0_r8/elsi_h%broadening_width
   diff_ne_out = 0.0_r8

   if((elsi_h%spin_degen /= 1.0_r8) .and. (elsi_h%spin_degen /= 2.0_r8)) then ! Not set
      if(elsi_h%n_spins == 2) then
         elsi_h%spin_degen = 1.0_r8
      else
         elsi_h%spin_degen = 2.0_r8
      endif
   endif

   select case(elsi_h%broadening_scheme)
   case(GAUSSIAN)
      do i_kpt = 1,elsi_h%n_kpts
         do i_spin = 1,elsi_h%n_spins
            do i_state = 1,elsi_h%n_states
               elsi_h%occ_num(i_state,i_spin,i_kpt) = elsi_h%spin_degen*0.5_r8*&
                  (1.0_r8-erf((elsi_h%eval_all(i_state,i_spin,i_kpt)-mu_in)*invert_width))

               diff_ne_out = elsi_h%occ_num(i_state,i_spin,i_kpt)*elsi_h%k_weight(i_kpt)+&
                                diff_ne_out
            enddo
         enddo
      enddo

   case(FERMI)
      max_exp = maxexponent(mu_in)*log(2.0_r8)

      do i_kpt = 1,elsi_h%n_kpts
         do i_spin = 1,elsi_h%n_spins
            do i_state = 1,elsi_h%n_states
               this_exp = (elsi_h%eval_all(i_state,i_spin,i_kpt)-mu_in)*invert_width

               if(this_exp < max_exp) then
                  elsi_h%occ_num(i_state,i_spin,i_kpt) = elsi_h%spin_degen/(1.0_r8+exp(this_exp))

                  diff_ne_out = elsi_h%occ_num(i_state,i_spin,i_kpt)*elsi_h%k_weight(i_kpt)+&
                                   diff_ne_out
               else ! Exponent in this step is larger than the largest possible exponent
                  elsi_h%occ_num(i_state,i_spin,i_kpt) = 0.0_r8
               endif
            enddo
         enddo
      enddo

   case(METHFESSEL_PAXTON_0)
      do i_kpt = 1,elsi_h%n_kpts
         do i_spin = 1,elsi_h%n_spins
            do i_state = 1,elsi_h%n_states
               elsi_h%occ_num(i_state,i_spin,i_kpt) = elsi_h%spin_degen*0.5_r8*&
                  (1.0_r8-erf((elsi_h%eval_all(i_state,i_spin,i_kpt)-mu_in)*invert_width))

               diff_ne_out = elsi_h%occ_num(i_state,i_spin,i_kpt)*elsi_h%k_weight(i_kpt)+&
                                diff_ne_out
            enddo
         enddo
      enddo

   case(METHFESSEL_PAXTON_1)
      do i_kpt = 1,elsi_h%n_kpts
         do i_spin = 1,elsi_h%n_spins
            do i_state = 1,elsi_h%n_states
               this_hermite = (elsi_h%eval_all(i_state,i_spin,i_kpt)-mu_in)*invert_width

               elsi_h%occ_num(i_state,i_spin,i_kpt) = elsi_h%spin_degen*0.5_r8*&
                  (1.0_r8-erf(this_hermite))-0.5_r8*invert_sqrt_pi*this_hermite*&
                  exp(-this_hermite*this_hermite)

               diff_ne_out = elsi_h%occ_num(i_state,i_spin,i_kpt)*elsi_h%k_weight(i_kpt)+&
                                diff_ne_out
            enddo
         enddo
      enddo
   end select

   diff_ne_out = diff_ne_out-elsi_h%n_electrons

end subroutine

!>
!! This routine computes the chemical potential using a bisection algorithm.
!!
subroutine elsi_find_mu(elsi_h,mu_lower_in,mu_upper_in)

   implicit none
 
   type(elsi_handle), intent(inout) :: elsi_h      !< Handle
   real(kind=r8),     intent(in)    :: mu_lower_in !< Lower bound of chemical potential
   real(kind=r8),     intent(in)    :: mu_upper_in !< Upper bound of chemical potential

   real(kind=r8)    :: mu_left    ! Left bound of chemical potential interval
   real(kind=r8)    :: mu_right   ! Right bound of chemical potential interval
   real(kind=r8)    :: mu_mid     ! Middle point of chemical potential interval
   real(kind=r8)    :: diff_left  ! Difference in number of electrons on left bound
   real(kind=r8)    :: diff_right ! Difference in number of electrons on right bound
   real(kind=r8)    :: diff_mid   ! Difference in number of electrons on middle point
   logical          :: found_mu
   integer(kind=i4) :: n_steps
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_find_mu"

   n_steps = 0
   found_mu = .false.

   mu_left = mu_lower_in
   mu_right = mu_upper_in

   do while((.not.found_mu) .and. (n_steps < elsi_h%max_mu_steps))
      call elsi_check_electrons(elsi_h,mu_left,diff_left)
      call elsi_check_electrons(elsi_h,mu_right,diff_right)

      if(abs(diff_left) < elsi_h%occ_tolerance) then
         elsi_h%mu = mu_left
         found_mu = .true.
      elseif(abs(diff_right) < elsi_h%occ_tolerance) then
         elsi_h%mu = mu_right
         found_mu = .true.
      else
         n_steps = n_steps+1

         mu_mid = 0.5_r8*(mu_left+mu_right)

         call elsi_check_electrons(elsi_h,mu_mid,diff_mid)

         if(abs(diff_mid) < elsi_h%occ_tolerance) then
            elsi_h%mu = mu_mid
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
      call elsi_check_electrons(elsi_h,mu_right,diff_right)

      elsi_h%mu = mu_right

      ! ...with adjusted occupation numbers
      call elsi_statement_print("  Chemical potential cannot reach the"//&
              " required accuracy by bisection method. The error will be"//&
              " arbitrarily canceled in the following state(s):",elsi_h)

      call elsi_adjust_occ(elsi_h,diff_right)
   endif

end subroutine

!>
!! This routine cancels the small error in number of electrons.
!!
subroutine elsi_adjust_occ(elsi_h,diff_ne)

   implicit none

   type(elsi_handle), intent(inout) :: elsi_h  !< Handle
   real(kind=r8),     intent(inout) :: diff_ne !< Error in electron count

   real(kind=r8),    allocatable :: eval_aux(:)
   real(kind=r8),    allocatable :: occ_aux(:)
   integer(kind=i4), allocatable :: idx_aux(:)

   real(kind=r8)    :: tmp_real
   integer(kind=i4) :: tmp_int
   integer(kind=i4) :: min_id
   integer(kind=i4) :: i_state
   integer(kind=i4) :: i_kpt
   integer(kind=i4) :: i_spin
   integer(kind=i4) :: i_val
   integer(kind=i4) :: n_total
   character*200    :: info_str

   character*40, parameter :: caller = "elsi_adjust_occ"

   n_total = elsi_h%n_states*elsi_h%n_spins*elsi_h%n_kpts

   call elsi_allocate(elsi_h,eval_aux,n_total,"eval_aux",caller)
   call elsi_allocate(elsi_h,occ_aux,n_total,"occ_aux",caller)
   call elsi_allocate(elsi_h,idx_aux,n_total,"idx_aux",caller)

   ! Put eigenvalues and occupation numbers into 1D arrays
   i_val = 0

   do i_kpt = 1,elsi_h%n_kpts
      do i_spin = 1,elsi_h%n_spins
         do i_state = 1,elsi_h%n_states
            i_val = i_val+1
            eval_aux(i_val) = elsi_h%eval_all(i_state,i_spin,i_kpt)
            occ_aux(i_val)  = elsi_h%occ_num(i_state,i_spin,i_kpt)
            idx_aux(i_val)  = i_val
         enddo
      enddo
   enddo

   ! Sort eval_aux; move occ_aux and idx_aux accordingly
   do i_val = 1,n_total
      min_id = minloc(eval_aux(i_val:n_total),1)+i_val-1

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

   call elsi_deallocate(elsi_h,eval_aux,"eval_aux")

   do i_val = n_total,1,-1
      if(occ_aux(i_val) > 0.0_r8) then
         i_kpt   = (idx_aux(i_val)-1)/(elsi_h%n_spins*elsi_h%n_states)+1
         i_spin  = mod((idx_aux(i_val)-1)/6,2)+1
         i_state = mod(idx_aux(i_val)-1,elsi_h%n_states)+1

         write(info_str,"(A,I5,A,I5,A,I5)") "  | k-point ",i_kpt,&
            ", spin channel ",i_spin,", state ",i_state
         call elsi_statement_print(info_str,elsi_h)

         if(elsi_h%k_weight(i_kpt)*occ_aux(i_val) > diff_ne) then
            occ_aux(i_val) = occ_aux(i_val)-diff_ne/elsi_h%k_weight(i_kpt)
            diff_ne = 0.0_r8
         else
            diff_ne = diff_ne-elsi_h%k_weight(i_kpt)*occ_aux(i_val)
            occ_aux(i_val) = 0.0_r8
         endif

         elsi_h%occ_num(i_state,i_spin,i_kpt) = occ_aux(i_val)
      endif

      if((diff_ne < elsi_h%occ_tolerance) .or. (diff_ne == 0.0_r8)) exit
   enddo

   call elsi_deallocate(elsi_h,occ_aux,"occ_aux")
   call elsi_deallocate(elsi_h,idx_aux,"idx_aux")

end subroutine

!>
!! This routine collects the eigenvalues of all spins and k-points.
!!
subroutine elsi_get_eval_all(elsi_h)

   implicit none
   include "mpif.h"

   type(elsi_handle), intent(inout) :: elsi_h  !< Handle

   real(kind=r8), allocatable :: tmp_real1(:,:,:)
   real(kind=r8), allocatable :: tmp_real2(:)
   integer(kind=i4)           :: mpierr

   character*40, parameter :: caller = "elsi_get_eval_all"

   if(.not. allocated(elsi_h%eval_all)) then
      call elsi_allocate(elsi_h,elsi_h%eval_all,elsi_h%n_states,&
              elsi_h%n_spins,elsi_h%n_kpts,"eval_all",caller)
   endif

   elsi_h%eval_all = 0.0_r8
   elsi_h%eval_all(1:elsi_h%n_states,elsi_h%i_spin,elsi_h%i_kpt) = &
      elsi_h%eval_elpa(1:elsi_h%n_states)

   if(.not. allocated(elsi_h%occ_num)) then
      call elsi_allocate(elsi_h,elsi_h%occ_num,elsi_h%n_states,&
              elsi_h%n_spins,elsi_h%n_kpts,"occ_num",caller)
   endif

   elsi_h%occ_num = 0.0_r8

   if(.not. allocated(elsi_h%k_weight)) then
      call elsi_allocate(elsi_h,elsi_h%k_weight,&
              elsi_h%n_kpts,"k_weight",caller)
   endif

   elsi_h%k_weight = 0.0_r8
   elsi_h%k_weight(elsi_h%i_kpt) = elsi_h%i_weight

   if(elsi_h%n_spins*elsi_h%n_kpts > 1) then
      call elsi_allocate(elsi_h,tmp_real1,elsi_h%n_states,&
              elsi_h%n_spins,elsi_h%n_kpts,"tmp_real",caller)

      call MPI_Allreduce(elsi_h%eval_all,tmp_real1,&
              elsi_h%n_states*elsi_h%n_spins*elsi_h%n_kpts,&
              mpi_real8,mpi_max,elsi_h%mpi_comm_all,mpierr)

      elsi_h%eval_all = tmp_real1

      call elsi_deallocate(elsi_h,tmp_real1,"tmp_real")

      call elsi_allocate(elsi_h,tmp_real2,elsi_h%n_kpts,"tmp_real",caller)

      call MPI_Allreduce(elsi_h%k_weight,tmp_real2,elsi_h%n_kpts,&
              mpi_real8,mpi_max,elsi_h%mpi_comm_all,mpierr)

      elsi_h%k_weight = tmp_real2

      call elsi_deallocate(elsi_h,tmp_real2,"tmp_real")
   endif

end subroutine

end module ELSI_MU
