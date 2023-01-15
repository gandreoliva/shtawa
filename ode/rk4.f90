subroutine rk4(f,y0,t_initial,t_final,niter,outfile,stopping_cond)
	!! Approximates the initial-value problem y'(t) = f(t,y(t)) with the initial
	!! condition y(t=t_initial) = y0

	!! Theory
	!! ------
	!! * Introduction to Runge-Kutta methods:
	!! Consider the Taylor expansion of the function y(t) up to second order:
	!!   y(t_{i+1}) = y(t_i) + h * y'(t_i) + (h^2/2)*y''(t(i)) + O(3)
	!!             = y(t_i) + h * [ f(t_i,y(t_i)) + (h/2)*f'(t_i,y(t_i)) ] + O(3)
	!!  By the chain rule,
	!!  f' = df/dt + df/dy*y'  = df/dt + df/dy*f
	!! (d: partial derivative, depend. omitted). We substitute in the sq. brackets []:
	!!    f + (h/2)*(df/dt) + (h/2)*(df/dy)*f ,  (dep. with (t_i,y(t_i)) omitted)
	!! Given that the Taylor expansion for a function g of two variables is
	!!   g(t+p,y+q) = f(t,y) + p*(df/dt) + q*(df/dy),
	!! we compare and find that one can express the term in sq. brackets as
	!!  f( t + h/2, y + (h/2)*f(t,y(t)) )
	!! So that
	!! y(t_{t+1}) = y(t_i) + h * [ f( t + h/2, y + (h/2)*f(t,y(t)) ) ].
	!! That is called the "Midpoint method".
	!! 
	!! * Runge-Kutta 4 method:
	!! From the previous discussion, we see that a) expansions of further order are possible,
	!! b) p and q are in principle arbitrary, so that different combinations give different
	!! methods. The Runge-Kutta 4 method is has a truncation error O(4).

	real(wp), dimension(:), intent(in) :: y0
		!! Initial condition y(t=t_initial) = y0
	real(wp), intent(in) :: t_initial, t_final
		!! Initial and final values of t
	integer, intent(in) :: niter
		!! Total number of iterations
	interface
		function f(t,y)
			!! Right-hand-side of the (system of) ODEs
			import wp
			real(wp), intent(in) :: t
			real(wp), dimension(:), intent(in) :: y
			real(wp), dimension(size(y)) :: f
		end function
		function stopping_cond_signature(t,y)
			!! Function that determines when to stop integrating
			import wp
			real(wp), intent(in) :: t
			real(wp), dimension(:), intent(in) :: y
			logical :: stopping_cond_signature
		end function
	end interface
	procedure(stopping_cond_signature), optional :: stopping_cond
		!! If present, this function provides a condition to be satisfied
		!! by the solution (when it is satisfied, the integration stops).
		!! If it's not present, the integration continues until the maximum
		!! number of iterations is reached.
	integer, intent(in) :: outfile
		!! Output file unit (default formatting)
		!! If the condition is reached, it prints a line with a space.

	integer :: i
	real(wp) :: h, t
	real(wp), dimension(:,:), allocatable :: k
	real(wp), dimension(:), allocatable :: y

	allocate(y(size(y0)))
	allocate(k(1:4,size(y0)))


	!! Step and initialization
	h = (t_final - t_initial)/niter
	t = t_initial
	y = y0

	!! Output initial step
	write(*,*) t, y


	evolve: do i = 1,niter
		k(1,:) = h*f(t,y)
		k(2,:) = h*f(t+h/2d0, y + k(1,:)/2d0)
		k(3,:) = h*f(t+h/2d0, y + k(2,:)/2d0)
		k(4,:) = h*f(t+h/2d0, y + k(3,:))

		y = y + (k(1,:) + 2*k(2,:) + 2*k(3,:) + k(4,:))/6d0
		t = t_initial + i*h

		if (present(stopping_cond)) then
			if (stopping_cond(t,y) .eqv. .true.) then
				write(*,*) " "
				exit evolve
			end if
		end if

		!! Output every step
		write(*,*) t, y

	end do evolve


end subroutine