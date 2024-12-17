function discrete_fourier_coeff(f)
    ! Computes the coefficients of the Fourier series expansion of the sampled function f
    ! assuming periodicity between -pi and pi.
    !
    ! Theory
    ! ------
    ! Let's assume we know N samples of a function f: Real->Complex with period 2*pi. The samples are
    ! equidistant at known at x_l = l * 2*pi/N, for l = -N/2,...,N/2.
    ! In the code, f is the sampled function f. We want to compute the Fourier coefficients of the
    ! Fourier series for this function.
    ! 
    ! Analytically, the Fourier coefficients are computed as
    !   c_k = 1/T  integral_{-T/2}^{T/2}  f(x) exp( i*k*2*pi/T * x) dx  (for k in Integers, T is the period)
    ! and this forms the series  f(x) = sum_{k=-inf}^{inf}  c_k * e^{i*k*2*pi/T * x}.
    ! 
    ! With our assumptions, the exponential is simplified to exp( -2*pi*i / N )**(k*l). We can approximate
    ! the coefficients for k = -N/2,..,0,..,N/2 with the sum
    !    c_k = (1/N) * sum_{l=-N/2}^{N/2} f(l) * exp( -2*pi*i / N )**(k*l)
    ! This works best with N odd (zero is in the middle of the indices and adds one).
    integer :: k, l
    real(wp), dimension(:) :: f
    real(wp), dimension(:), allocatable :: f_ordered
    complex(wp), dimension(:), allocatable :: discrete_fourier_coeff
    complex(wp), parameter :: iimag = (0d0,1d0)
    real(wp), parameter :: pi = 3.141592653589793
    integer :: n

    n = size(f)
    if (mod(n,2) == 0) then
        n = n-1
        print*, "Warning: even number of samples. Last sample ignored!"
    end if
    allocate(discrete_fourier_coeff(-n/2:n/2))
    allocate(f_ordered(-n/2:n/2))
    
    f_ordered = f(1:n)
    
    do k=-n/2,n/2
        discrete_fourier_coeff(k) = 0
        do l=-n/2,n/2
            discrete_fourier_coeff(k) = discrete_fourier_coeff(k) + f_ordered(l) * exp(-2*pi*iimag/real(n,wp))**(k*l)
        end do
        discrete_fourier_coeff(k) = discrete_fourier_coeff(k)*(1/real(n,wp))
    end do

end function


function discrete_fourier_transform(f) result(fourier)
    ! Computes the discrete Fourier transform of the complex sampled function f assuming
    ! periodicity between 0 and 2*pi.
    ! This does the same as the function discrete_fourier_coeff, but f here is complex and
    ! the indices are ordered starting from zero instead of negative values.
    !
    ! Theory (continuation of discrete_fourier_coeff)
    ! ------
    ! Let's assume we know N samples of a function f: Complex->Complex with period 2*pi. The samples are
    ! now equidistant at known at x_l = l * 2*pi/N, for l = 0,...,N-1.
    ! 
    ! Analytically, the Fourier coefficients defined within 0 and the period T are
    !   c_k = 1/T  integral_{0}^{T}  f(x) exp( i*k*2*pi/T * x) dx  (for k in Integers, T is the period)
    ! We can write the Fourier coefficients  k = 0,...,N-1 with the sum
    !    c_k = (1/N) * sum_{l=0}^{n-1} f(l) * exp( -2*pi*i / N )**(k*l)
    !  Note that one can also write the exponential part of the Fourier coefficients
    ! as a matrix zeta(k,l) = exp( -2*pi*i / N )**(k*l)
    ! and then this whole operation of computing the coefficients is reduced to a matrix product.
    ! Now we notice that the Fourier transform is defined as
    !    F(w) =  Cauchy principal value of integral_{-inf}^{inf}  f(t) * e^{-i*w*t) dt
    ! For a periodic, sampled function we take the integral only within one period, and we notice
    ! that this integral is the same as the Fourier coefficients for the Fourier series. For historical
    ! reasons and numerical convenience, the indices here start at zero.
    
    complex(wp), dimension(0:) :: f
    complex(wp), dimension(:,:), allocatable :: zeta
    complex(wp), dimension(:), allocatable :: fourier
    complex(wp), parameter :: iimag = (0d0,1d0)
    complex(wp) :: zetaconstant
    real(wp), parameter :: pi = 3.141592653589793
    integer :: k,l, n

    n = size(f)
    allocate(zeta(0:n,0:n))
    allocate(fourier(0:n))

    zetaconstant = exp(-2*pi*iimag/n)
    zeta = reshape([((zetaconstant**(k*l), l=0,n-1), k=0,n-1)] ,[n,n])
    
    fourier = 1/real(n,wp) * matmul(zeta,f)

end function




function discrete_inverse_fourier_transform(f) result(invfourier)
    ! Computes the discrete inverse Fourier transform of the complex sampled function f.
    ! This function does the same as discrete_fourier_transform, but one
    ! takes the complex conjugate of zetaconstant (i.e., the exponential with a + instead of a -)
    ! and the normalization factor is also gone.
    
    
    complex(wp), dimension(0:) :: f
    complex(wp), dimension(:,:), allocatable :: zeta
    complex(wp), dimension(:), allocatable :: invfourier
    complex(wp), parameter :: iimag = (0d0,1d0)
    complex(wp) :: zetaconstant
    real(wp), parameter :: pi = 3.141592653589793
    integer :: k,l, n

    n = size(f)
    allocate(zeta(0:n,0:n))
    allocate(invfourier(0:n))

    zetaconstant = exp(2*pi*iimag/n)
    zeta = reshape([((zetaconstant**(k*l), l=0,n-1), k=0,n-1)] ,[n,n])
    
    invfourier = matmul(zeta,f)

end function