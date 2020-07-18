subroutine single_hh_trafo_real_cpu_single(q, hh, nb, nq, ldq)

   use elpa_abstract_impl
   use precision
   ! Perform single real Householder transformation.
   ! This routine is not performance critical and thus it is coded here in Fortran

   implicit none

   integer(kind=ik), intent(in)   :: nb, nq, ldq
   real(kind=rk4), intent(inout)   :: q(ldq, *)
   real(kind=rk4), intent(in)      :: hh(*)
   integer(kind=ik)               :: i
   real(kind=rk4)                  :: v(nq)

   ! v = q * hh
   v(:) = q(1:nq,1)
   do i=2,nb
      v(:) = v(:) + q(1:nq,i) * hh(i)
   enddo

   ! v = v * tau
   v(:) = v(:) * hh(1)

   ! q = q - v * hh**T
   q(1:nq,1) = q(1:nq,1) - v(:)
   do i=2,nb
      q(1:nq,i) = q(1:nq,i) - v(:) * hh(i)
   enddo

end subroutine
