subroutine secant(f,x_approx0,x_approx1,x_sol,tol,maxiter,istat)
    !! Secant method for root finding.
    !! * Theory
    !! based on Newton's method
    !!  $ p_n = p_{n-1} - \frac{p_{n-1}}{f'(p_{n-1})} $,
    !! but numerically computing the derivative of f, that is
    !!  $f'(p_{n-1}) = \frac{f(p_{n-1}) - f(p_{n-2})}{p_{n-1} - p_{n-2}}$
    !! so that it's no longer necessary to provide f_deriv, but in turn it's necessary
    !! to provide two approximations instead of one.

    real(dp), intent(in) :: x_approx0, x_approx1
        !! initial approximations to the solution
    real(dp), intent(out) :: x_sol
        !! approximation to the solution
    real(dp), intent(in) :: tol
        !! tolerance to the solution
    integer, intent(in) :: maxiter
        !! maximum number of iterations allowed
    interface
        real(dp) function single_var_func_signature(x)
            !! Signature that f must have
            import :: dp
            real(dp), intent(in) :: x
        end function
    end interface
    procedure(single_var_func_signature) :: f
        !! single variable function f(x) for finding the roots
    integer, intent(out) :: istat
    !! end state of the algorithm. 0: success (x_sol valid). 1: maximum number of iterations
    !!  reached without a solution; 2: f'(p_0) = 0, i.e., requirement not satisfied
    integer :: i
    real(dp) :: p,p0,p1,f_deriv_recipr

    p0 = x_approx0
    p1 = x_approx1

    !  check requirements
    if (abs(f(p1) - f(p0)) < tol) then
        istat = 2
        return
        ! requirement f'(p_0) = 0 not satisfied!
    end if

    istat = -1
    i = 1
    do while(i <= maxiter)
        ! approximation to the solution
        f_deriv_recipr = (p1 - p0)/(f(p1) - f(p0))
        p = p0 - f(p0)*f_deriv_recipr
        ! check whether the solution is good enough
        if (abs(p-p0) < tol) then
            !  solution successfully found
            x_sol = p
            istat = 0
            return
        else
            ! next improvement to the ansatz
            i = i + 1
            p0 = p1
            p1 = p
        end if
    end do


    if (istat /= 0) then
       istat = 1
        ! No solution found after max. num. of iterations
    end if

end subroutine