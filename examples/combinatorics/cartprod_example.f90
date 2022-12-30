program cartprod_example
    use iso_fortran_env, only: dp => real64
    use shtawa
    implicit none
    integer, parameter :: n = 3
    integer, dimension(1:n) :: sets_ubounds,ntuple
    integer, dimension(:,:), allocatable :: cartprod
    integer :: i,nproducts

    sets_ubounds = [2,4,1]
    
    print*, "Recursive general version:"
    call get_ncartprod(sets_ubounds,cartprod)

    nproducts = size(cartprod,2)
    do i=1,nproducts
        print*, cartprod(:,i)
    end do
    print*, "----"
    print*, nproducts, "products in total"


    block
        integer :: i,j,k
        nproducts = 0
        print*, "========="
        print*, "Iterative version:"
        do i=1,sets_ubounds(1)
            do j=1,sets_ubounds(2)
                do k=1,sets_ubounds(3)
                    print*, i,j,k
                    nproducts = nproducts + 1
                end do
            end do
        end do
        print*, nproducts

    end block


end program