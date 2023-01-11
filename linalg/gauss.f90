subroutine linsys_gauss(a,x,istat)
	!! Solves the linear system represented by a(:,:n)*x == a(:,n+1),
	!!  a(1,1)*x(1) + a(1,2)*x(2) + ... = a(1,n+1)
	!!  a(2,1)*x(1) + a(2,1)*x(2) + ... = a(2,n+1)
	!!  ...
	!!  a(1,n)*x(1) + a(2,n)*x(2) + ... = a(n,n+1)
	!! by means of Gauss elimination with backward substitution.
	!! Elemental operations on a matrix m:
	!!      1. Multiply row i by a constant k: m(i,:) = k*m(i,:)
	!!      2. Sum two rows: m(i,:) = m(i,:) + k*m(j,:)
	!!      3. Swap two rows: swap(m(i,:),m(j,:))

	!! Theory
	!! ------
	!! ```
	!! Example 1:
	!!  6*x1 - 3*x2 = 0    =>  a = | 6  -3  0   |     n = 2
	!!  5*x1 + x2   = 7            | 5  1   7   |
	!!  The goal of Gaussian elimination is to find a diagonal matrix,
	!!  i.e., modify 'a' (using only elemental operations)
	!!  so that it becomes a matrix of the form
	!!      | a11   a12    b1  |    =>    a11*x1 + a12*x2 = b1
	!!      | 0     a22    b2  |                   a22*x2 = b2
	!!  because now, we can easily read x2 = a22*b2, and
	!!      x1 = ( b1 - a12*x2 )/a11 = ( b1 - a12*(a22*b2) )/a11 .
	!!  This last operation is called "backwards substitution".
	!!  To diagonalize the original matrix in the example, we start with row 1.
	!!  All remaining k rows (in this case only row 2) must have ak1 = 0.
	!!  We can achieve this by performing the operation row2 - a21/a11*row1 --> row2
	!!  because the first term of row k will be ak1 - ak1/a11*a11 = 0. In
	!!  our example, after doing this operation we obtain
	!!     | 6  -3      0 |
	!!     | 0  7/2     7 |
	!!  which is already diagonal. Then backwards substitution gives us x2=2, x1=1.
	!!
	!! Example 2:
	!!        3*y = 6   =>  a = |0  3   6|    n = 2
	!!  2*x +   y = 4           |2  1   4|
	!!  We start with row 1. However, here the element we want to be in the diagonal
	!! (a11) is zero. So, first we find an element in column 1 which is nonzero,
	!!  which is the element a21, and then we exchange rows 2 and 1. After that,
	!! | 2   1   4 |
	!! | 0   3   6 |
	!! which is already diagonal. If we apply the rest of the algorithm for Gaussian
	!! elimination, we see that row2 - a21*a11*row1 = row2 + 0 --> row2  doesn't change
	!! the matrix. Backwards substitution gives us x2=2, x1=1.
	!! ```
	real(wp), dimension(1:), intent(out) :: x
		!! 1d array with the solution to the unknowns: x = [x1,x2,x3,...]
	real(wp), dimension(1:,1:), intent(inout) :: a
		!! augmented 2d matrix with the coefficients (columns) of each equation (rows).
		!! If n is the number of equations, a should have shape (n,n+1), with the last
		!! column being the coefficients not multiplied by an unknown.
	integer, intent(out) :: istat
		!! 0: success, 1: no unique solution exists
	integer :: i,k,n  ,l
	real(wp) :: multiplier

	istat = -1
	n = size(x)

	!! For debugging: print matrix
	! print*, "---"
	! do l = 1,n
	!     write(*,"(*(f8.2))") a(l,:)
	! end do

	! for all the rows except the last one...
	do i=1,n-1
		!! Checking whether a row exchange is necessary (diag. elem. must be nonzero)
		!! Find a non-zero element in the remaining rows of the same column

		!! Without pivoting: simply find the closest non-zero value
		! k = i
		! do while ((a(k,i) == 0) .and. (k < n+1))
		!     k = k + 1
		! end do
		! if ( k > n ) then ! if all the coeficients in the column are zero
		!     istat = 1
		!     ! No unique solution exists
		!     return
		! end if

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

		! Gaussian elimination: apply row operation to remaining rows (see Theory > Example 1)
		! to make matrix diagonal
		do k=i+1,n
			multiplier = a(k,i)/a(i,i)
			a(k,:) = a(k,:) - multiplier*a(i,:)
		end do
	end do

	if (a(n,n) == 0) then
		istat = 1
		! No unique solution exists
		return
	end if

	! Backward substitution
	x(n) = a(n,n+1)/a(n,n)

	i = n-1
	do while( i >= 1 )
		x(i) = (a(i,n+1) - sum(a(i,i+1:n)*x(i+1:)) )/a(i,i)
		i = i-1
	end do

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