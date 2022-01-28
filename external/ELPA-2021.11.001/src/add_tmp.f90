











module add_tmp
  use precision
  implicit none
  private

  public :: add_tmp_double
  public :: add_tmp_single

  contains

! real double precision first
















subroutine add_tmp_&
&double&
&(obj, d1, dbase, ddiff, z, ev_scale_value, na1,i)
  use precision
  use v_add_s
  use elpa_abstract_impl
  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in) :: na1, i

  real(kind=rk8), intent(in)    :: d1(:), dbase(:), ddiff(:), z(:)
  real(kind=rk8), intent(inout) :: ev_scale_value
  real(kind=rk8)                :: tmp(1:na1)

  ! tmp(1:na1) = z(1:na1) / delta(1:na1,i)  ! original code
  ! tmp(1:na1) = z(1:na1) / (d1(1:na1)-d(i))! bad results

  ! All we want to calculate is tmp = (d1(1:na1)-dbase(i))+ddiff(i)
  ! in exactly this order, but we want to prevent compiler optimization

  tmp(1:na1) = d1(1:na1) -dbase(i)
  call v_add_s_&
  &double&
  &(obj, tmp(1:na1),na1,ddiff(i))
  tmp(1:na1) = z(1:na1) / tmp(1:na1)
  ev_scale_value = 1.0_c_double/sqrt(dot_product(tmp(1:na1),tmp(1:na1)))

end subroutine add_tmp_&
&double


! real single precision first


















subroutine add_tmp_&
&single&
&(obj, d1, dbase, ddiff, z, ev_scale_value, na1,i)
  use precision
  use v_add_s
  use elpa_abstract_impl
  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in) :: na1, i

  real(kind=rk4), intent(in)    :: d1(:), dbase(:), ddiff(:), z(:)
  real(kind=rk4), intent(inout) :: ev_scale_value
  real(kind=rk4)                :: tmp(1:na1)

  ! tmp(1:na1) = z(1:na1) / delta(1:na1,i)  ! original code
  ! tmp(1:na1) = z(1:na1) / (d1(1:na1)-d(i))! bad results

  ! All we want to calculate is tmp = (d1(1:na1)-dbase(i))+ddiff(i)
  ! in exactly this order, but we want to prevent compiler optimization

  tmp(1:na1) = d1(1:na1) -dbase(i)
  call v_add_s_&
  &single&
  &(obj, tmp(1:na1),na1,ddiff(i))
  tmp(1:na1) = z(1:na1) / tmp(1:na1)
  ev_scale_value = 1.0_c_float/sqrt(dot_product(tmp(1:na1),tmp(1:na1)))

end subroutine add_tmp_&
&single


end module
