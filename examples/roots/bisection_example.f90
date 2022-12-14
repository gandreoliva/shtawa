program example_bisection
    use iso_fortran_env, only: dp => real64
    use shtawa, only: bisection
    real(dp) :: x_sol, interv_beg, interv_end
    integer :: istat
    real(dp) :: x
    integer :: i,n,prev_sign

    interv_beg = -10d0
    interv_end = 10d0

    write(*,*) "Checking for sign changes in the interval [", interv_beg, ",", interv_end, "]"
    n = 100 ! subdivisions of the interval to check
    prev_sign = int(sign(1d0,f(interv_beg)))
    do i = 1,n
        x = interv_beg + i*(interv_end-interv_beg)/n
        if ( int(sign(1d0,f(x))) /= prev_sign ) then
            print*, "f(x) changes sign at x ~ ", x
        end if
        prev_sign = int(sign(1d0,f(x)))
    end do


    

    write(*,*) "Example of the bisection subroutine for root finding"

    call bisection(f=f, interv_beg=1.2d0, interv_end=1.8d0, x_sol=x_sol, &
        & tol=1d-8, maxiter=50, istat=istat)
    
    if (istat /= 0) error stop "Bisection method failed!"

    write(*,*) "Solution: x_sol = ", x_sol

contains
    function f(x)
        real(dp), intent(in) :: x
        real(dp) :: f
        f = x**3 - 4d0*x**2+ 0.5d0*x + 4
        ! f = x**2 + 1
        ! f = x - 3
    end function
end program