subroutine newton(f,f_deriv,x_approx,x_sol,tol,maxiter,istat)
	!! Newton(-Raphson) method for root finding.
	
	!! Theory
	!! ------
	!! Newton's method aims to solve the one-variable equation f(x) == 0, if we know the
	!! derivative of f, f'(x), and an initial guess for the solution.
	
	!! It's based on a Taylor series expansion around the solution x_sol and a given initial
	!! approximation x_approx. Consider the problem $f(x) = 0$. (f(x) and its derivative
	!! should be continuous and differentiable over an interval). Let $p_0$ be an approximation
	!! to the solution x_sol:= $p$ such that $f'($p_0$) \neq 0$. Now, consider the Taylor expansion
	!!      $ 0 = f(p) = f(p_0) + (p - p_0) f'(p_0) + O(2) $
	!! which means that we can solve for p:
	!!      $ p = p_0 - \frac{p_0}{f'(p_0)} $.
	!! $|p-p_0|$ is assumed to be small (too off approx. might not converge), which means that
	!! if we start with an approximation p_0, it can be improved iteratively.
	real(wp), intent(in) :: x_approx
		!! initial approximation to the solution
	real(wp), intent(out) :: x_sol
		!! approximation to the solution
	real(wp), intent(in) :: tol
		!! tolerance to the solution
	integer, intent(in) :: maxiter
		!! maximum number of iterations allowed
	interface
		real(wp) function single_var_func_signature(x)
			!! Signature that f must have
			import :: wp
			real(wp), intent(in) :: x
		end function
	end interface
	procedure(single_var_func_signature) :: f
		!! single variable function f(x) for finding the roots
	procedure(single_var_func_signature) :: f_deriv
		!! derivative of f(x), i.e., f'(x)
	integer, intent(out) :: istat
	!! end state of the algorithm. 0: success (x_sol valid). 1: maximum number of iterations
	!!  reached without a solution; 2: f'(p_0) = 0, i.e., requirement not satisfied
	integer :: i
	real(wp) :: p,p0

	p0 = x_approx

	!  check requirements
	if (abs(f_deriv(p0)) < tol) then
		istat = 2
		return
		! requirement f'(p_0) = 0 not satisfied!
	end if

	istat = -1
	i = 1
	do while(i <= maxiter)
		! approximation to the solution
		p = p0 - f(p0)/f_deriv(p0)
		! check whether the solution is good enough
		if (abs(p-p0) < tol) then
			!  solution successfully found
			x_sol = p
			istat = 0
			return
		else
			! next improvement to the ansatz
			i = i + 1
			p0 = p
		end if
	end do


	if (istat /= 0) then
	   istat = 1
		! No solution found after max. num. of iterations
	end if


end subroutine