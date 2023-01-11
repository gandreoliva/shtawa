program secant_example
	use iso_fortran_env, only: dp => real64
	use shtawa, only: secant
	real(dp) :: x_sol
	real(dp), parameter :: pi = 3.14159265359d0
	integer :: istat

	write(*,*) "Example of the Secant subroutine for root finding"

	call secant(f, x_approx0=pi/4, x_approx1=0d0, x_sol=x_sol, &
		& tol=1d-8, maxiter=50, istat=istat)

	if (istat /= 0) error stop "Secant method failed"

	write(*,*) "Solution: x_sol = ", x_sol
	

contains
	real(dp) function f(x)
		real(dp), intent(in) :: x
		f = cos(x) - x
	end function
end program