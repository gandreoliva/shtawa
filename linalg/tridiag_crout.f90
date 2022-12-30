subroutine tridiag_crout(a,x)
    !! Solves the system a(:,:n)*x == a(:,n+1), where a(:,:n) is (band) tridiagonal,
    !! a(:,:n) =    | a11   a12                    0        |
    !!              | a21   a22     a23                     |
    !!              |       a32     a33     a34             |
    !!              |        ...    ...     ...             |
    !!              |                               a(n-1,n)|
    !!              |      0               a(n,n-1)  a(n,n) |
    !!  by using Crout factorization.

    !! Theory
    !! ------
    !! Notation: in this section, a := a(:,:n); b := a(:,n+1); @ is matrix multiplication.
    !! Some matrices, among them tridiagonal matrices, can be LU factorized, i.e., 'a'
    !! be written as  a = l @ u, where l is lower-triangular and u upper-triangular,
    !!
    !!          | l11   0   0   0 |         | 1 u12 0   0 |
    !!   l =    | l12   l22 0   0 |     u = | 0 1   .   0 |
    !!          | 0     .   .   0 |         | 0 0   1   . |
    !!          | 0     0   .   . |         | 0 0   0   1 |
    !!  To find l and u from a, we simply perform by hand the matrix multiplication
    !!  a = l @ u and find that the non-zero elements are
    !!      a(1,1) = l(1,1)
    !!      a(i,i) = l(i,i-1)*u(i-1,i) + l(i,i)    for i = 2,3,...,n
    !!      a(i,i-1) = l(i,i-1)                    for i = 2,3,...,n
    !!      a(i,i+1) = l(i,i)*u(i,i+1)             for i = 1,2,...,n-1 ,
    !! from which we can solve for each element of l and u given a, as done in the
    !! algorithm.
    !! 
    !! The advantage of doing this is that the original system, a @ x == b, i.e.,
    !!  (l @ u) @ x == b can now be solved in two steps with fewer calculations.
    !! First, we define z := u @ x. Then, we first solve l @ z == b for z.
    !! Because l is lower triangular, we can do forward substitution. By performing
    !! the matrix multiplication l @ z == b by hand, we find
    !!  l(1,1)*z(1) = b(1) ==> z(1) = b(1)/l(1,1)
    !!  l(i,i-1)*z(i-1) + l(i,i)*z(i) = b(i) ==> z(i) given z(i-1), etc.
    !! Once z is determined, then we can solve u @ x == z by backwards substitution
    !! and the problem is solved.

    real(dp), dimension(:,:), intent(inout) :: a
        !! augmented tridiagonal matrix to be solved, shape (n,n+1). It is overwritten
        !! after calling the subroutine, so make sure you make a copy beforehand if needed!
    real(dp), dimension(:), intent(out) :: x
        !! solution vector, shape (n)
    integer :: i,n
    real(dp), dimension(:), allocatable :: z

    n = size(a,1)
    allocate(z(1:n))

    !! We note that l and u don't overlap any non-zero and non-one components.
    !! Because of this, we can store their non-overlapping elements in one matrix.
    !! Additionally, we have built the method such that the values of overlapping
    !! elements (0 or 1) have been already substituted in the equations, and we only
    !! use each element of 'a' once. Because of this, we can directly use a to store
    !! l and u instead of using extra memory. For clarity, the 'associate' construct
    !! adds an alias to l and u so that the method can be followed easily, but in reality
    !! l and u point to the same memory space as 'a'. This means that 'a' is overwritten
    !! and shouldn't be used after the subroutine has been called. In case the user wants
    !! to use a, they need to make a copy before calling the subroutine.
    associate (l => a, u => a)
        
        !! Part I: Building l and u, and solving l @ z == b. Remember: b = a(:,n+1)
        
        !! First element
        l(1,1) = a(1,1)
        u(1,2) = a(1,2)/l(1,1)
        z(1) = a(1,n+1)/l(1,1)

        do i=2,n-1
            !! Build the elements of l left to the diagonal
            l(i,i-1) = a(i,i-1)
            !! Build the diagonal elements of l
            l(i,i) = a(i,i) - l(i,i-1)*u(i-1,i)
            !! Build the elements of u on top of the diagonal
            u(i,i+1) = a(i,i+1)/l(i,i)
            !! Forward substitution to solve for z
            z(i) = (a(i,n+1) - l(i,i-1)*z(i-1))/l(i,i)
        end do

        !! Final elements
        l(n,n-1) = a(n,n-1)
        l(n,n) = a(n,n) - l(n,n-1)*u(n-1,n)
        z(n) = (a(n,n+1) - l(n,n-1)*z(n-1))/l(n,n)

        !! Part II: Solution of u @ x == z: backwards substitution
        x(n) = z(n)
        i = n-1
        do while(i >= 1)
            x(i) = z(i) - u(i,i+1)*x(i+1)
            i = i-1
        end do

    end associate

end subroutine