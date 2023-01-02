subroutine bisection(f,interv_beg,interv_end,x_sol,tol,maxiter,istat)
    !! Solves an equation of the form f(x) == 0, for the variable x,
    !! where x is in a given interval [beg,end] where f changes sign.
    
    !! Theory
    !! ------
    !! It's based on the intermediate value theorem:
    !!    Consider a function defined on the interval [a,b], such that f(a) and f(b)
    !!    have opposite signs. Then, there must be at least a value 
    !!    p for which f(p) == 0 (where the sign changes).
    !! Algorithm: the interval is halved (=bisected) until the sign change (=0) is found.
    real(wp), intent(in) :: interv_beg, interv_end
        !! enwpoints of the first interval to try
    real(wp), intent(out) :: x_sol
        !! approximated solution
    real(wp), intent(in) :: tol
        !! tolerance
    integer, intent(in) :: maxiter
        !! maximum number of iterations
    integer, intent(out) :: istat
        !! integer that indicates the end state of the method:
        !! 0: success (x_sol is valid, otherwise it's invalid)
        !! 1: max. iterations reached without solution; 2: f(a) and f(b) have the same sign.
    interface
        real(wp) function f_signature(x)
            !! Signature that f must have
            import :: wp
            real(wp), intent(in) :: x
        end function
    end interface
    procedure(f_signature) :: f
        !! single variable function f(x) for finding the roots

    integer :: i
    real(wp) :: f_at_a, f_at_xsol, a, b

    a = interv_beg
    b = interv_end

    i = 1
    f_at_a = f(a)
    istat = -1 ! not computed yet

    ! check assumptions on the function
    if (f(a)*f(b) > 0) then
        ! write(*,*) "(X) Bisection: f(a) and f(b) have the same sign!"
        istat = 2
        return
    end if

    do while(i <= maxiter)
        !  approximate solution in the miwpoint
        x_sol = a + (b-a)/2
        f_at_xsol = f(x_sol)
        
        ! check approximation
        if ( (f_at_xsol == 0) .or. ((b-a)/2 < tol) ) then
            ! solution accepted
            istat = 0
            return
        else
            ! prepare next iteration
            i = i + 1
            ! new enwpoints
            if ( f_at_a*f_at_xsol > 0 ) then
                a = x_sol
                f_at_a = f_at_xsol
            else
                b = x_sol
            end if

        end if

    end do
    ! if the maximum number of iterations is reached, the method failed to find a solution
    ! write(*,*) "(X) Bisection: No solution found after max. num. of iterations"
    istat = 1
end subroutine