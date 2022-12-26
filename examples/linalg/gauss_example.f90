program gauss_example
    use iso_fortran_env, only: dp => real64
    use shtawa
    implicit none
    real(dp), dimension(:,:), allocatable :: a,acopy
    real(dp), dimension(:), allocatable :: x
    integer :: istat, n,i

    ! Examples with static arrays a(1:2,1:3), x(1:2)
    ! a = reshape( [ &
    !         6d0,  -3d0,   0d0, &
    !     &   5d0,  1d0,  7d0 &
    !     &] ,shape=[2,2+1], order=[2,1])
    !
    ! a = reshape( [ &
    !         0d0,  3d0,   6d0, &
    !     &   2d0,  1d0,  4d0 &
    !     &] ,shape=[2,2+1], order=[2,1])


    !! Interactive of the Gaussian elimination method for solving linear systems of equations

    print*, "> Enter number of equations:"
    read(*,*) n
    
    allocate(a(n,n+1))
    allocate(x(n))
    print*, "> Enter augmented matrix a(1:n,1:n+1):"
    !!  An augmented matrix is built such that
    !!  such that each eq. i is 
    !!  a(i,1)*x(1) + ... a(i,n)*x(n) = a(i,n+1)
    do i=1,n
        read(*,*) a(i,:)
    end do

    !! To verify the solution later, we need to make a copy of a. The subroutine
    !! linsys_gauss modifies the value of 'a'.
    allocate(acopy(n,n+1))
    acopy(:,:) = a(:,:)
    
    call linsys_gauss(a,x,istat)

    if (istat /= 0) then
        print*, "Error code: ", istat
        error stop "(X) Error in the solution"
    end if

    print*, "Solution: x = "
    print*, x


    print*, "Verification of the solution: a @ x - b should be 0"
    print*, [(sum(acopy(i,:n)*x(:)) ,i=1,n)] -  acopy(:,n+1)



end program