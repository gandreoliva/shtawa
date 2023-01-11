subroutine linsys_gaussjordan(a,x,istat)
	!! Solves the linear systems of equations represented by
	!!  a(1,1)*x(1)+a(1,2)*x(2)+...= a(1,n+1), a(1,n+2),..., a(1,n+m)
	!!  a(2,1)*x(1)+a(2,1)*x(2)+...= a(2,n+1), a(2,n+2),..., a(2,n+m)
	!!  ...
	!!  a(1,n)*x(1)+a(2,n)*x(2)+...= a(n,n+1), a(n,n+2),..., a(n,n+m)
	!! by means of Gauss Jordan elimination. This is, the systems
	!! a*x == b1, a*x == b2,..., a*x == bm are solved simultaneously.
	!! This subroutine generalizes x to be a 2d array of shape (n,m),
	!! where m is given in the shape of the matrix a. For m = 1, x is
	!! a column vector shape(n,1).
	!!
	!!  Theory
	!!  -----
	!!  The same as Gauss's method but Gauss elimination is also made on the rows
	!!  above the diagonal, not only below. This results in a unit echelon matrix
	!!  | a11   0       0       b1... |
	!!  | 0     a22     0       b2... |
	!!  | 0     0       a33     b3... |
	!! after which the solution is easily read in col. n+1 by dividing each row 
	!! over the diagonal.
	real(wp), dimension(:,:), allocatable, intent(out) :: x
		!! 2d array, shape (n,m) with the solution vectors of each one of the m equations.
	real(wp), dimension(1:,1:), intent(inout) :: a
		!! augmented 2d matrix with the coefficients (columns) of each equation (rows).
		!! It has a shape (n,n+m).
	integer, intent(out) :: istat
		!! 0: success, 1: no unique solution exists
	integer :: i,k,n,m,l
	real(wp) :: multiplier

	n = size(a,1)
	m = size(a,2) - n
	allocate(x(n,m))

	istat = -1

	!! For debugging: print matrix
	! print*, "---"
	! do l = 1,n
	!     write(*,"(*(f8.2))") a(l,:)
	! end do

	! for all the rows except the last one...
	do i=1,n-1
		!! Checking whether a row exchange is necessary (diag. elem. must be nonzero)
		!! Find a non-zero element in the remaining rows of the same column

		! !! With partial pivoting: exchange rows with the maximum non-zero value.
		! !! This reduces some numerical rounding errors.
		! !! Look for maxima(abs) in this and the remaining rows of the same column.
		! !! the (i-1) converts back to the row index of the full matrix.
		k = maxloc(abs(a(i:,i)),dim=1) + (i-1)
		if (a(k,i) == 0) then ! if the maximum is zero => all coefficients are zero
			istat = 1
			! No unique solution exists
			return
		end if

		if (k /= i) then ! a(i,i) == 0, so a row switch is necessary
			call swap(a(k,:),a(i,:))
		end if

		! Gaussian elimination: apply row operation to rows below the diagonal
		do k=i+1,n
			multiplier = a(k,i)/a(i,i)
			a(k,:) = a(k,:) - multiplier*a(i,:)
		end do
	end do



	i = n
	do while(i > 1)
		! Gaussian elimination: apply row operation to rows above the diagonal
		do k=1,i-1
			multiplier = a(k,i)/a(i,i)
			a(k,:) = a(k,:) - multiplier*a(i,:)
		end do
		i = i - 1
	end do


	if (a(n,n) == 0) then
		istat = 1
		! No unique solution exists
		return
	end if

	! Create unit matrix in a(:,:n), to read the solution in a(:,n+m), by rescaling each row
	do i = 1,n
		a(i,:) = a(i,:)/a(i,i)
	end do

	x = a(:,n+1:n+m)

	istat = 0 ! success

contains
	elemental subroutine swap(a1,a2)
		!! swaps the values of two reals (or arrays)
		real(wp), intent(inout) :: a1,a2
		real(wp) :: temp
		temp = a1
		a1 = a2
		a2 = temp
	end subroutine
end subroutine