recursive function determinant(a) result(res)
    !! Finds the determinant of a matrix a
    !! Algorithm from rosettacode.org
    !! 
    !! Theory
    !! ------
    !! Consider an example 3x3 matrix
    !! a =  | a11   a12     a13 |
    !!      | a21   a22     a23 |
    !!      | a31   a32     a33 |
    !!  The determinant is built recursively by taking the first row
    !!  and building a submatrix b for each element a(1,j). b does not contain
    !!  the first row of a, nor the column j of a. Then, the determinant of b
    !!  is computed. The determinant of a is then computed as
    !!   det(a) = +a(1,1)*det(b(-1,-1)) - a(1,2)*det(b(-1,-2)) + a(1,3)*det(b(-1,-3))
    !!  i.e., with a sum with alternating signs (even columns have a negative sign),
    !!  and where b(-i,-j) means the submatrix without row i and column j.

    real(dp), dimension(1:,1:), intent(in) :: a
        !! (square, 2d) matrix to compute the determinant
    real(dp), dimension(:,:), allocatable :: b
    integer :: sgn, n, j
    real(dp) :: res

    n = size(a,1)
    allocate(b(n-1,n-1))

    if (n == 1) then
        res = a(1,1)
    else
        res = 0
        sgn = 1
        do j = 1,n
            ! Build submatrix: what comes before column j
            b(:, :j-1) = a(2:, :j-1)
            ! Build submatrix: what comes after column j
            b(:, j:) = a(2:, j+1:)
            ! Recursive call
            res = res + sgn * a(1, j) * determinant(b)
            ! Alternation of signs
            sgn = -sgn
        end do
    end if

end function