module omm_params

implicit none

public

!**** PARAMS ************************************!

integer, parameter :: dp=selected_real_kind(15,300)

real(dp), parameter :: Pi=3.141592653589793238462643383279502884197_dp

complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

integer, parameter :: i64 = selected_int_kind(18)

!************************************************!

end module omm_params
