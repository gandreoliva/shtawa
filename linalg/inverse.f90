function inverse(a)
    !! Computes the inverse of the matrix 'a' using the Gauss-Jordan procedure
    !! Warning: if 'a' is not invertible, the function throws an error and stops
    real(wp), dimension(:,:), intent(in) :: a
    real(wp), dimension(:,:), allocatable :: a_augm, inverse, unit
    integer :: i,n, istat
    n = size(a,1)

    allocate(a_augm(n,2*n))
    allocate(inverse(n,n))
    allocate(unit(n,n))

    a_augm(:,:n) = a

    ! Building unit matrix at the right side of the augmented matrix
    a_augm(:,n+1:2*n) = 0d0
    do i=1,n
        a_augm(i,i+n) = 1d0
    end do

    call linsys_gaussjordan(a_augm,inverse,istat)
    if (istat /= 0) then
        write(*,*) "(X) Error: inverse could not be computed. Status code =",istat
        error stop
    end if

end function