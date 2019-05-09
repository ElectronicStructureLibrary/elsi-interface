module single_hh_trafo_real

  implicit none

  public :: single_hh_trafo_real_cpu_double

  contains

! Perform single real Householder transformation.
! This routine is not performance critical and thus it is coded here in Fortran
  subroutine single_hh_trafo_real_cpu_double(q, hh, nb, nq, ldq)

    use precision

    implicit none

    integer(kind=ik), intent(in) :: nb, nq, ldq
    real(kind=rk8), intent(inout) :: q(1:ldq, 1:nb)
    real(kind=rk8), intent(in) :: hh(1:nb)
    integer(kind=ik) :: i
    real(kind=rk8) :: v(nq)

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

end module
