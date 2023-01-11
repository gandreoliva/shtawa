subroutine build_cubic_spline(x_grid,f_grid,spline_coefs,fp_1,fp_n)
	!! Builds a cubic spline for a 1d grid. Word origin: 'spline' is a kind of flexible
	!! ruler used in the past to draw curves.

	!! Theory
	!! -----
	!! The basic idea of cubic interpolation is the following: we take a grid of positions
	!! x_grid, and field values f_grid and fit cubic polinomials in each grid point
	!!   f|    .o·`·o
	!!    | o.'           o: point in the grid
	!!    -------------x
	!! Each one of the polynomials has an equation S_i = k3*x**3 + k2*x**2 + k1*x + k0,
	!! and it is continuous and differentiable at every grid point. For the differentia-
	!! bility condition, we need to know the values of the derivatives of the polynomial
	!! at each boundary of the grid. A condition of f'(1) = f'(n) = 0 is called
	!! "natural spline" and it means we assume the polynomial has zero gradient (flat) on
	!! the extremes. To build the spline, we notice that the second derivative of a
	!! third-degree polynomial is a linear function (=> to find each S''_i we can use
	!! the solution of a linear system of equations). Upon integration, S'_i and S_i
	!! are found; the integration constants are set by the conditions. This subroutine
	!! finds S''_i.

	!! Implementation
	!! --------------
	!! This function is meant to be only called once to build the spline. The function
	!! spline_interpolate does the actual interpolation and should be called after
	!! calling the spline coefficients.


	!! References
	!! ----------
	!! * Stoer (1976) Einfuehrung in die Numerische Mathematik I, Springer. (S. 2.4.2)

	real(wp), dimension(:), intent(in) :: x_grid
		!! grid of x values
	real(wp), dimension(:), intent(in) :: f_grid
		!! Function f evaluated at every grid point
	real(wp), dimension(:), intent(out) :: spline_coefs
		!! Second derivatives of the spline. To feed into spline_interpolate
	real(wp), intent(in), optional :: fp_1
		!! If present, value of f'(1). If not present, a natural spline is used for x(1)
	real(wp), intent(in), optional :: fp_n
		!! If present, value of f'(n). If not present, a natural spline is used for x(n)
	integer :: i,n
	real(wp), dimension(:), allocatable :: a_l, a_d, a_r, b

	n = size(x_grid)
	allocate(a_l(2:n),a_d(1:n),a_r(1:n-1),b(1:n))

	!! Tridiagonal system of linear equations
	!!   a_l(i) * S''(i-1) + a_d(i) * S''(i)  + a_r(i) * S''(i+1) = b(i)
	!! Let h(i) = x(i) - x(i-1). For i=2,n-1, the matrix elements are given by
	!!   a_l(i) = h(i)
	!!   a_r(i) = h(i+1) = a_l(i+1)
	!!   a_d(i) = 2*( h(i) + h(i+1) ) = 2*(a_l(i) + a_r(i))
	!!   b(i) = 6 * ( (f(i+1) - f(i))/h(i+1) - (f(i) - f(i-1))/h(i) ) = b`(i) - b`(i-1)
	
	!! For i=2,..,n-1 (values of i=1 will be overwritten ahead)
	a_r(1:n-1) = x_grid(2:n) - x_grid(1:n-1)
	!! a_l(i) = a_r(i-1)
	a_l(2:n-1) = a_r(1:n-2)
	b(1:n-1) = 6d0*( f_grid(2:n) - f_grid(1:n-1) )/a_r(1:n-1)
	b(2:n-1) = b(2:n-1) - b(1:n-2)
	a_d(2:n-1) = 2d0*( a_l(2:n-1) + a_r(2:n-1) )

	!! Boundaries: i=1, i=n

	!! For a natural spline, a_r(1) = 0;  d(1) = 0; ==> 1*S''(1) = 0
	if (.not. present(fp_1)) then
		a_d(1) = 1
		a_r(1) = 0
		b(1) = 0
	else
		!! For specified derivatives, i=1
		!!      b(1) = 6*( (f(2) - f(1))/h(2) - f'(1) )/h(2)
		a_d(1) = 2
		a_r(1) = 1
		b(1) = 6d0*( (f_grid(2)-f_grid(1))/(x_grid(2)-x_grid(1)) - fp_1 )&
				&/(x_grid(2)-x_grid(1))
	end if

	!!   a_l(n) = 0; b(n) = 0; ==> 1*S''(n) = 0
	if (.not. present(fp_n)) then
		a_l(n) = 0
		a_d(n) = 1
		b(n) = 0
	else
		!! For specified derivatives, i=n
		!!      b(n) = 6*( f'(n) - (f(n) - f(n-1))/(x(n)-x(n-1)) )/(x(n)-x(n-1))
		a_l(n) = 1
		a_d(n) = 2
		b(n) = 6d0*(fp_n - (f_grid(n) - f_grid(n-1))/(x_grid(n)-x_grid(n-1)) )&
						 &/(x_grid(n)-x_grid(n-1))
	end if

	!! Solve for the spline coefficients
	call tridiag(a_l,a_d,a_r,b,spline_coefs)

end subroutine


function spline_interpolate(x,x_grid,f_grid,spline_coefs) result(f)
	!! Interpolates a value by evaluating a cubic spline, i.e., approximates f(x)
	!! given tabulated values of x and f in a grid, using the method of cubic spline
	!! interpolation.
	real(wp), intent(in) :: x
		!! Value of x for which f(x) is desired
	real(wp), dimension(:), intent(in) :: x_grid
		!! Grid with x values
	real(wp), dimension(:), intent(in) :: f_grid
		!! Values of f evaluated/known at each grid point
	real(wp), dimension(:), intent(in) :: spline_coefs
		!! Coefficients of the cubic spline obtained with build_cubic_spline
	real(wp) :: f
		!! Interpolated value of x
	integer :: il,n
	integer, dimension(1:1) :: ind_left
	real(wp) :: alp,bet,gam,del,h_right,x_diff

	n = size(x_grid)

	call locate_points_in_grid([x],x_grid,ind_left)
	il = ind_left(1)
	if ((il == 0) .or. (il == n)) error stop "Point not found in grid!"

	!! Evaluation of the spline
	h_right = x_grid(il+1) - x_grid(il)
	alp = f_grid(il)
	bet = (f_grid(il+1)-f_grid(il))/h_right - &
		& h_right*(2*spline_coefs(il) + spline_coefs(il+1))/6d0
	gam = spline_coefs(il)/2d0
	del = (spline_coefs(il+1) - spline_coefs(il))/(6*h_right)

	x_diff = x - x_grid(il)
	f = alp + bet*x_diff + gam*x_diff**2 + del*x_diff**3

end function