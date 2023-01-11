program multilinear_interpolate_example
	use iso_fortran_env, only: dp => real64
	use shtawa
	implicit none
	integer, parameter :: nx=50,ny=40,nz=10
	type(grid), dimension(3) :: grids
	real(dp), dimension(nx) :: x_grid
	real(dp), dimension(ny) :: y_grid
	real(dp), dimension(nz) :: z_grid
	real(dp), dimension(nx,ny,nz) :: f_grid
	integer :: i
	real(dp) :: f
	real(dp), dimension(3) :: point

	x_grid = [(i*0.05d0,i=1,nx)]
	y_grid = [(i*0.25d0,i=1,ny)]
	z_grid = [(-0.5d0+(i-1)*0.1d0,i=1,nz)]

	grids(1)%points = x_grid
	grids(2)%points = y_grid
	grids(3)%points = z_grid

	generate_data: block
		integer :: j,k
		do i=1,nx
			do j=1,ny
				do k=1,nz
					f_grid(i,j,k) = sin(x_grid(i))*y_grid(j)+z_grid(k)**2
				end do
			end do
		end do
	end block generate_data
	
	point = [0.85d0,0.413d0,0.284d0]
	call multilinear_interpolate_point(point,grids,f_grid,f)

	print*, "point =",point
	print*, "Interpolated value:",f

	print*, "True value:",sin(point(1))*point(2)+point(3)**2

end program