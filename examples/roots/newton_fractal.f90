program complex_newton_example
    use iso_fortran_env, only: dp => real64
    use shtawa
    implicit none
    integer, parameter :: nx = 500, ny = nx
    complex(dp), parameter :: iimag = (0,1) ! z = x + y*iimag = r*exp(iimag*th)
    complex(dp), dimension(0:nx,0:ny) :: z_grid
    complex(dp), dimension(:), allocatable :: est_zeros
    real(dp), parameter :: epsilon = 1d-1, l = 10d0
    integer, dimension(0:nx,0:ny) :: colors
    integer :: i,j,istat
    complex(dp) :: z_sol
    real(dp), dimension(:), allocatable :: dist_sol_zeros
    integer :: file_colors
    character(100), parameter :: filename_colors = "data/newton_fractal.dat"
    character(1000) :: pycmd

    ! Part I: Estimate the zeros of the function in a regular grid by minimizing the norm
    ! of f(z) (small norm => f(z) ~ 0). This algorithm only finds the minimum values within a
    ! grid's precision. It doesn't guarantee that there is a root near an estimated zero.
    call estimate_zeros(epsilon,-l,l,-l,l,z_grid,est_zeros)
    colors = -1
    print*, "Roots estimated at", est_zeros

    open(newunit=file_colors, file=filename_colors, form="unformatted",&
        & access="stream", status="replace")

    ! Part II: generation of the fractal.
    ! Solve f(z) == 0 with the Newton's method for all cells in the grid.
    ! In general, Newton's method doesn't tell us to which of the many roots it converges to.
    ! The following algorithm classifies the roots yielded by Newton's method by min. distance to
    ! the estimated solutions from Part I (which were generated on a regular grid).
    ! The goal is to produce a map of the complex plane such that, if we start Newton's method at
    ! a given z_approx, we know to which of the multiple solutions of f(z) == 0 
    ! Newton's method converges to.

    do i=0,nx
        do j=0,ny

            call complex_newton(f=f, f_deriv=f_deriv, z_approx=z_grid(i,j), z_sol=z_sol,&
                & tol=1d-2,maxiter=1000,istat=istat)
            
            ! if the solution is valid and there are estimated solutions from Part I...
            if ((istat == 0) .and. (size(est_zeros)>0)) then
                ! decide which color corresponds to each point in the grid, i.e.,
                ! classify each grid point according to which solution the Newton's method converges
                dist_sol_zeros = sqrt((z_sol-est_zeros)*conjg(z_sol-est_zeros))
                colors(i,j) = minloc(dist_sol_zeros,dim=1)
            end if

        end do
    end do

    write(file_colors) colors
    close(file_colors)

    ! plot the fractal
    write(pycmd,*) "python newton_fractal_plot.py ", filename_colors, " ", nx+1, " ", ny+1
    ! print*, pycmd
    call execute_command_line(pycmd)

contains

    ! Function for finding the roots (don't forget to update derivative!)
    complex(dp) function f(z)
        complex(dp), intent(in) :: z
        f = z**3 - z**2 + 2
    end function

    ! Derivative of f(z)
    complex(dp) function f_deriv(z)
        complex(dp), intent(in) :: z
        f_deriv = 3*z**2 - 2*z
    end function

    subroutine estimate_zeros(epsilon,x_beg, x_end, y_beg, y_end,z_grid,est_zeros)
            !! Finds points where the function f is close to zero, assuming that they are
            !! point-like in the complex plane.
            !! This is true for polynomials or functions like sin(z) or cos(z).
            !! epsilon: threshold below which we consider a minimum to be a zero
        real(dp), intent(in) :: epsilon, x_beg, x_end, y_beg, y_end
        complex(dp), dimension(0:nx,0:ny), intent(out) :: z_grid
        complex(dp), dimension(:), allocatable, intent(out) :: est_zeros
        real(dp), dimension(0:nx,0:ny) :: f_r_grid
        complex(dp), dimension(0:nx,0:ny) :: f_grid
        integer :: i,j
        real(dp) :: dx,dy

        dx = (x_end-x_beg)/nx
        dy = (y_end-y_beg)/ny

        ! Generate grid
        z_grid = reshape(  [( ( (x_beg + i*dx) + iimag*(y_beg + j*dy), j=0,ny ), i=0,nx )]   &
                        &, shape=[nx+1,ny+1]  )

        ! Evaluate the function and compute r (sqrt of norm). f(z) ~ 0 has f_r ~ 0.
        f_grid = reshape( [ (( f(z_grid(i,j)) ,j=0,ny), i=0,nx) ] , shape = [nx+1,ny+1])
        f_r_grid = sqrt(f_grid*conjg(f_grid))

        est_zeros = [0d0] ! temporary initialization of the list with a first element
        ! Find local minima, i.e., where f_r ~ 0 but also smaller than its neighbors in cross
        do i=1,nx-1
            do j=1,ny-1
                if ((f_r_grid(i,j) <= epsilon) .and.        &
                &   (f_r_grid(i,j) < f_r_grid(i-1,j)) .and. &
                &   (f_r_grid(i,j) < f_r_grid(i+1,j)) .and. &
                &   (f_r_grid(i,j) < f_r_grid(i,j-1)) .and. &
                &   (f_r_grid(i,j) < f_r_grid(i,j+1))       &
                & ) then
                    ! add to the list
                    est_zeros = [est_zeros, z_grid(i,j)]
                end if
            end do
        end do
        est_zeros = est_zeros(2:) ! get rid of first element
    end subroutine

end program

