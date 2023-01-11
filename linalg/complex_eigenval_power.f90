subroutine complex_eigenval_power(a,x,eigenval,tol,niter,istat)
	!! Power method to approximate the dominant eigenvalue and associated eigenvector
	!! of the matrix a.

	!! Theory
	!! ------
	!! In general, the eigenvalues l_i and eigenvectors v of a matrix 'a' are defined
	!! via  a @ v == l_i * v. Seeing matrices as transformations of vectors, we say that 
	!! the matrix 'a' simply scales v by l_i, without changing its direction.
	!! 
	!! Let 'a' be a n x n diagonalizable matrix with eigenvalues l_i. For some matrices,
	!! there is a 'dominant' eigenvalue l1 larger than the rest, i.e.,
	!!  |l1| > |l2| >= |l3| >= ... >= |l_n|
	!! The power method aims to find l1 and v1 by succesively applying 'a' to an inital
	!! vector x, which it must not be orthogonal to v1. In practice, we can choose x
	!! as a random vector to reduce the chances that the choice is orthogonal to v1.
	!! The set of eigenvectors is an orthogonal basis to the matrix's linear space.
	!! this means that any vector x can be decomposed as
	!!          x = c1*v1 + c2*v2 + ... + c_n*v_n
	!! applying 'a' one time,
	!!      a @ x = c1*l1*v1 + c2*l2*v2 + ... + c_n*l_n*v_n .
	!! After successive applications of 'a', the first term with l1 will grow faster
	!! because l1 is larger than the rest of eigenvalues. Because of this as well,
	!! successive applications of 'a' turn 'x' into something similar to v1. To avoid
	!! both sides of the equation to become large after succesive applications of 'a',
	!! we normalize the vector. Customary, the L_inf norm is used for general matrices
	!! and the L2 norm is used for symmetric matrices. The L_inf norm of x is essentially
	!! the maximum value of the absolute value (or complex norm) of its entries.
	complex(wp), dimension(:,:), intent(in) :: a
		!! Matrix for which the eigenvalue and eigenvector is to be computed. Shape (1:n,1:n)
	complex(wp), dimension(:,:), intent(inout) :: x
		!! Initial vector when it goes in, approximation to the dominant eigenvector 
		!! when it goes out. It should have shape (1:n,1) (column vector)
	complex(wp), intent(out) :: eigenval
		!! Approximation to the dominant eigenvalue
	real(wp), intent(in) :: tol
		!! tolerance of the solution
	integer, intent(in) :: niter
		!! Maximum number of iterations
	integer, intent(out) :: istat
		!! 0: sucess; 1: maximum number of iterations reached; 2: eigenval = 0, another vector
		!! should be selected.
	complex(wp), parameter :: iimag = (0d0,1d0)
	integer :: k,n,normloc
	complex(wp), dimension(:,:), allocatable :: x_next !! next approx for x
	real(wp) :: diff_convergence

	n = size(a,1)
	allocate(x_next(n,1))

	k = 1
	istat = -1

	!! Find which entry contains the L_inf norm of x, and normalize x
	normloc = maxloc(abs(x(:,1)),dim=1)
	x = x/x(normloc,1)
	
	do while(k <= niter)
		x_next = matmul(a,x)
		eigenval = x_next(normloc,1)

		!! Normalize x_next
		normloc = maxloc(abs(x_next(:,1)),dim=1)

		if (x_next(normloc,1) == 0d0+0*iimag) then
			istat = 2
			!! The eigenvalue is zero. Another x should be selected.
		end if

		diff_convergence = maxval(abs( x(:,1) - x_next(:,1)/x_next(normloc,1) ))
		x = x_next/x_next(normloc,1)

		if (diff_convergence < tol) then
			istat = 0
			!! Success
			return
		else
			!! Prepare next iteration
			k = k + 1
		end if
	end do

	istat = 1
	!! If this point is reached, there was no convergence after the max num. of iterations.

end subroutine