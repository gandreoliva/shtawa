program tridiag_example
    use iso_fortran_env, only: dp => real64
    use shtawa
    implicit none
    integer, parameter :: n = 4
    real(dp), dimension(n) :: a_d, b, x
    real(dp), dimension(n-1) :: a_r, a_l

    ! Solution of the system
    ! 2   -1  0   0   1
    ! -1  2   -1  0   0
    ! 0   -1  2   -1  0
    ! 0   0   -1  2   1

    a_d = 2
    a_l = -1
    a_r = -1
    b = [1d0, 0d0, 0d0, 1d0]

    call tridiag(a_l,a_d,a_r,b,x)

    print*, "Solution by tridiag:"
    print*, x
    
    verification: block
        real, dimension(n,n) :: a_full
        integer i

        ! Remember: the tridiag subroutine overwrites the matrix!
        ! Reinitialization
        a_d = 2
        a_l = -1
        a_r = -1

        a_full = 0
        a_full(1,1) = a_d(1)
        a_full(1,2) = a_r(1)
        do i=2,n-1
            a_full(i,i-1) = a_l(i)
            a_full(i,i) = a_d(i)
            a_full(i,i+1) = a_r(i)
        end do
        a_full(n,n) = a_d(n)
        a_full(n,n-1) = a_l(n-1)

        print*, "a_full"
        do i=1,n
            print*, a_full(i,:)
        end do

        print*, "Verification (a @ x - b == 0)"
        print*, matmul(a_full,reshape(x,[n,1],order=[2,1])) - reshape(b,[n,1],order=[2,1])
    end block verification

end program