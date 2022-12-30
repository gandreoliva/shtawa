program combinations_example
    use iso_fortran_env, only: dp=> real64
    use shtawa
    implicit none
    integer, dimension(:,:), allocatable :: combinations
    integer :: i,n,k
    
    n = 3
    k = 2
    call get_combinations_slow(n,k,combinations)

    do i=1,size(combinations,2)
        print*, combinations(:,i)
    end do
    
    print*, "Computed", size(combinations,2), "combinations, which should be C(n,k)"
    print*, "C(n,k) =", int(gamma(real(n+1)))/(int(gamma(real(k+1)))*int(gamma(real(n-k+1))))

end program