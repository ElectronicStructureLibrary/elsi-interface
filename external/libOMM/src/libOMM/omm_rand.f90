module omm_rand

implicit none

contains

    !================================================!
    ! random number generator                        !
    ! -initialize with:                              !
    !  call rand_init()                              !
    ! -generate new number with:                     !
    !  call random_number(rn)                        !
    !  where where rn is a real(dp) variable         !
    !================================================!
    function omm_rand_seed() result(seed)

        implicit none
        integer :: seed
        seed=123456

    end function omm_rand_seed

    subroutine omm_bsd_lcg(x, r)

        ! x_{n+1} = (a * x_{n} + c) mod m
        ! r = x_{n+1) / m

        use omm_params, only: i64, dp

        implicit none

        integer, intent(inout) :: x
        real(dp), intent(out) :: r
        integer(i64), parameter :: a = 1103515245_i64
        integer, parameter :: c = 12345
        integer(i64), parameter :: m = 2_i64**31

        x = int(mod(a*x+c,m))
        r = real(x,dp)/m

    end subroutine omm_bsd_lcg

end module omm_rand
