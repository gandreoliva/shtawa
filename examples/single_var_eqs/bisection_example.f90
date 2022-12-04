program example_bisection
    use iso_fortran_env, only: dp => real64
    use shtawa, only: bisection
    real(dp) :: x_sol
    integer :: istat

    write(*,*) "Test of the function call: ", f(0d0)

    write(*,*) "Example of the bisection subroutine for root finding"

    call bisection(f=f, interv_beg=-10d0, interv_end=10d0, x_sol=x_sol, &
        & tol=0.01d0, maxiter=1000, istat=istat)
    if (istat == 0) then
        write(*,*) "Solution: x_sol = ", x_sol
    end if

contains
    function f(x)
        real(dp), intent(in) :: x
        real(dp) :: f
        f = x**3 - 4d0*x**2+ 0.5d0*x + 4
        ! f = x**2 + 1
        ! f = x - 3
    end function
end program