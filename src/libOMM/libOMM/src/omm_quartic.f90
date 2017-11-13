!================================================!
! calculate the coeffs. of the quartic line      !
! search equation from three energy points and   !
! two gradient points                            !
!================================================!
subroutine omm_fit_quartic(x,y,g,c)
  use omm_params, only : dp

  implicit none

  !**** INPUT ***********************************!

  real(dp), intent(in) :: x(1:3) ! three x-points {x_i}
  real(dp), intent(in) :: y(1:3) ! y(x_i) at the three points
  real(dp), intent(in) :: g(1:2) ! (dy/dx)|x_i at the first two points

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: c(0:4) ! coeffs. of the quartic equation

  !**********************************************!

  ! the following expressions for the coeffs. were produced automatically using Maple 12
  c(4)=(x(3)**3*x(2)*g(1)-3*x(1)*x(2)**2*y(1)+3*y(3)*x(1)*x(2)**2+x(1)**2*x(2)**2*g(1)+x(3)*x(2)**3*&
       g(1)+2*x(1)**2*x(3)**2*g(2)-3*x(2)*x(3)**2*y(1)+3*y(2)*x(1)**2*x(2)-x(3)**3*x(1)*g(1)+x(3)**3*&
       x(2)*g(2)-x(2)**2*x(3)**2*g(2)-x(1)**2*x(2)**2*g(2)-2*x(2)**2*x(3)**2*g(1)+3*x(2)*x(3)**2*&
       y(2)+x(1)**2*x(3)**2*g(1)-x(1)*x(2)**3*g(1)-3*x(1)*x(3)**2*y(1)-x(3)**3*x(1)*g(2)+3*x(1)*&
       x(3)**2*y(2)+x(2)*g(2)*x(1)**3-3*y(3)*x(1)**2*x(2)-x(3)*g(2)*x(1)**3+6*x(1)*x(3)*x(2)*y(1)+2*&
       x(1)*x(3)*x(2)**2*g(2)+x(1)*x(3)**2*x(2)*g(1)-x(1)*x(3)**2*x(2)*g(2)-2*x(1)**2*x(3)*x(2)*g(1)-&
       x(1)**2*x(3)*x(2)*g(2)+x(1)*x(3)*x(2)**2*g(1)-6*x(1)*x(3)*x(2)*y(2)+2*x(3)**3*y(1)-2*x(3)**3*&
       y(2)+x(2)**3*y(1)+y(3)*x(1)**3-y(3)*x(2)**3-y(2)*x(1)**3)/(-2*x(3)**3*x(1)**4+x(3)**4*x(1)**3-&
       x(1)**2*x(2)**5-x(3)**4*x(2)**3-x(2)**5*x(3)**2-3*x(1)**4*x(2)**3+2*x(3)**3*x(2)**4+x(1)**5*&
       x(3)**2+3*x(2)**4*x(1)**3+4*x(3)**3*x(2)*x(1)**3-4*x(3)**3*x(1)*x(2)**3+2*x(1)*x(3)*x(2)**5+4*&
       x(1)**4*x(3)*x(2)**2+8*x(1)**2*x(3)**2*x(2)**3+x(1)**4*x(3)**2*x(2)-x(1)*x(3)**2*x(2)**4-2*&
       x(1)**5*x(3)*x(2)-4*x(1)**2*x(3)*x(2)**4+x(1)**5*x(2)**2-8*x(2)**2*x(3)**2*x(1)**3-3*x(3)**4*&
       x(1)**2*x(2)+3*x(3)**4*x(1)*x(2)**2)

  c(3)=-(-x(1)*g(1)+2*c(4)*x(1)**4-x(1)*g(2)+4*x(1)*c(4)*x(2)**3+x(2)*g(1)-4*x(2)*c(4)*x(1)**3+x(2)*&
       g(2)-2*c(4)*x(2)**4+2*y(1)-2*y(2))/(x(1)**3+3*x(1)*x(2)**2-3*x(2)*x(1)**2-x(2)**3)

  c(2)=-(-y(2)+c(4)*x(2)**4+c(3)*x(2)**3+x(2)*g(1)-4*x(2)*c(4)*x(1)**3-3*x(2)*c(3)*x(1)**2+y(1)+3*&
       c(4)*x(1)**4+2*c(3)*x(1)**3-x(1)*g(1))/(x(1)**2-2*x(1)*x(2)+x(2)**2)

  c(1)=g(1)-4*c(4)*x(1)**3-3*c(3)*x(1)**2-2*c(2)*x(1)

  c(0)=y(1)-c(4)*x(1)**4-c(3)*x(1)**3-c(2)*x(1)**2-c(1)*x(1)

  !if (Node==0) print*, 'f(x)=',c(4),'*x**4+',c(3),'*x**3+',c(2),'*x**2+',c(1),'*x+',c(0)

end subroutine omm_fit_quartic

!================================================!
! find the minimum for the quartic line search   !
! equation                                       !
!================================================!
subroutine omm_solve_quartic(c,x_min,fail)
  use omm_params, only : dp, Pi

  implicit none

  !**** INPUT ***********************************!

  real(dp), intent(in) :: c(0:4) ! coeffs. of the quartic equation


  !**** OUTPUT **********************************!

  logical, intent(out) :: fail ! did we fail to find a minimum?

  real(dp), intent(out) :: x_min ! position of minimum

  !**** LOCAL ***********************************!

  integer :: x_order(1:3)

  real(dp) :: t(1:3)
  real(dp) :: z(1:3)
  real(dp) :: a
  real(dp) :: b
  real(dp) :: d
  real(dp) :: Q
  real(dp) :: R
  real(dp) :: theta
  real(dp) :: S
  real(dp) :: U

  !**********************************************!

  fail=.false.

  ! in order to find the minimum of the quartic equation, we have to solve a cubic equation; the
  ! following method is taken from Numerical Recipes
  a=3.0_dp*c(3)/(4.0_dp*c(4))
  b=2.0_dp*c(2)/(4.0_dp*c(4))
  if ((abs(b)>=1.0d11) .or. (abs(c(4))<=1.0d-11)) then
    x_min=-0.5_dp*c(1)/c(2)
    return
  end if
  d=c(1)/(4.0_dp*c(4))

  Q=(a**2-3.0_dp*b)/9.0_dp
  R=(2.0_dp*a**3-9.0_dp*a*b+27.0_dp*d)/54.0_dp
  if (R**2<Q**3) then
    theta=acos(R/sqrt(Q**3))
    t(1)=-2.0_dp*sqrt(Q)*cos(theta/3.0_dp)-a/3.0_dp
    t(2)=-2.0_dp*sqrt(Q)*cos((theta+2.0_dp*Pi)/3.0_dp)-a/3.0_dp
    t(3)=-2.0_dp*sqrt(Q)*cos((theta-2.0_dp*Pi)/3.0_dp)-a/3.0_dp
    z(1:3)=c(4)*t(1:3)**4+c(3)*t(1:3)**3+c(2)*t(1:3)**2+c(1)*t(1:3)+c(0)
    if (c(4)>0.0_dp) then
      if (all(z(1)>=z(2:3))) then
        x_order(1:3)=(/1,2,3/)
      else if (z(2)>z(3)) then
        x_order(1:3)=(/2,3,1/)
      else
        x_order(1:3)=(/3,1,2/)
      end if
      if ((0.0_dp<=t(x_order(1))) .and. (t(x_order(2))<=t(x_order(1)))) then
        x_min=t(x_order(2))
      else
        x_min=t(x_order(3))
      end if
    else
      if (all(z(1)<=z(2:3))) then
        x_min=t(1)
      else if (z(2)<z(3)) then
        x_min=t(2)
      else
        x_min=t(3)
      end if
    end if
  else
    S=-sign(1.0_dp,R)*(abs(R)+sqrt(R**2-Q**3))**(1.0_dp/3.0_dp)
    if (S==0.0_dp) then
      U=0.0_dp
    else
      U=Q/S
    end if
    x_min=(S+U)-(a/3.0_dp)
    if (c(4)<0.0_dp) fail=.true.
  end if

end subroutine omm_solve_quartic
