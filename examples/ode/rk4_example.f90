program rk4_example
	use iso_fortran_env, only: output_unit, dp=> real64
	use shtawa
	implicit none
	real(dp) :: c

	c = 0.8

	!! y0 = [ y(1,t=0), y(2,t=0) ] = [ x'(t=0), x(t=0) ]
	call rk4(f,y0=[0d0,0.3d0],t_initial=0d0,t_final=15d0,&
				& niter=100,outfile=output_unit,stopping_cond=stopping_cond)

contains
	function f(t,y)
		real(dp), intent(in) :: t
		real(dp), dimension(:), intent(in) :: y
		real(dp), dimension(size(y)) :: f
		!! System of equations (oscillations)
		!! 	x''(t) = -c*x
		!! ==> v'(t) = -c*x,  x'(t) = v
		!! ==> y'(1) = -c*y(2), y'(2) = y(1)
		
		f(1) = -c*y(2)
		f(2) = y(1)
	end function

	function stopping_cond(t,y)
		real(dp), intent(in) :: t
		real(dp), dimension(:), intent(in) :: y
		logical :: stopping_cond
		logical :: conditions
		stopping_cond = .false.
		conditions = t > 4d0
		conditions = conditions .and. y(2) > 0d0  ! position
		conditions = conditions .and. y(1) < 0d0  ! velocity
		if ( conditions .eqv. .true. ) stopping_cond = .true.
	end function

end program