program example_bisection
    use iso_fortran_env, only: dp => real64
    use shtawa, only: newton
    real(dp) :: x_sol
    real(dp), parameter :: pi = 3.14159265359d0
    integer :: istat

    write(*,*) "Example of the Newton-Raphson subroutine for root finding"

    call newton(f, f_deriv, x_approx=pi/4, x_sol=x_sol, &
        & tol=0.0001d0, maxiter=1000, istat=istat)

    if (istat == 0) then
        write(*,*) "Solution: x_sol = ", x_sol
    else
        write(*,*) "Failed with status code ", istat
    end if

contains
    real(dp) function f(x)
        real(dp), intent(in) :: x
        f = cos(x) - x
    end function
    real(dp) function f_deriv(x)
        real(dp), intent(in) :: x
        f_deriv = -sin(x) - 1
    end function
end program