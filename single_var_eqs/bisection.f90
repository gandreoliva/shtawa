subroutine bisection(f,interv_beg,interv_end,x_sol,tol,maxiter,istat)
    use iso_fortran_env, only: dp => real64
    !! Solves an equation of the form f(x) == 0, for the variable x (i.e., finds x_sol).
    !! * Theory: Based on the intermediate value theorem:
    !!    Consider a function defined on the interval [a,b], such that f(a) and f(b)
    !!    have opposite signs. Then, there must be at least a value 
    !!    p for which f(p) == 0 (where the sign changes).
    !! * Algorithm: the interval is halved (=bisected) until the sign change (=0) is found.
    real(dp), intent(in) :: interv_beg, interv_end
        !! endpoints of the first interval to try
    real(dp), intent(out) :: x_sol
        !! approximated solution
    real(dp), intent(in) :: tol
        !! tolerance
    integer, intent(in) :: maxiter
        !! maximum number of iterations
    integer, intent(out) :: istat
        !! integer that indicates the end state of the method:
        !! 0: success (x_sol is valid, otherwise it's invalid)
        !! 1: max. iterations reached without solution; 2: f(a) and f(b) have the same sign.
    interface
        real(dp) function f_signature(x)
            !! Signature that f must have
            import :: dp
            real(dp), intent(in) :: x
        end function
    end interface
    procedure(f_signature) :: f
        !! single variable function f(x) for finding the roots

    integer :: i
    real(dp) :: f_at_a, f_at_xsol, a, b

    a = interv_beg
    b = interv_end

    i = 1
    f_at_a = f(a)
    istat = 1 ! not computed yet

    ! check assumptions on the function
    if (f(a)*f(b) > 0) then
        write(*,*) "(X) Bisection: f(a) and f(b) have the same sign!"
        istat = 2
        return
    end if

    do while(i < maxiter)
        !  approximate solution in the midpoint
        x_sol = a + (b-a)/2
        f_at_xsol = f(x_sol)
        !  if the required tolerance is reached, finish
        if ( (f_at_xsol == 0) .or. ((b-a)/2 < tol) ) then
            istat = 0
            return
        end if
        i = i + 1
        
        !  for next iteration, compute new endpoints
        if ( f_at_a*f_at_xsol > 0 ) then
            a = x_sol
            f_at_a = f_at_xsol
        else
            b = x_sol
        end if

    end do
    ! if the maximum number of iterations is reached, the method failed to find a solution
    write(*,*) "(X) Bisection: No solution found (max. num. of iterations reached)"
    istat = 1
end subroutine