program fourier_transform_example
    !! (slow) Fourier transform of a discrete, periodic, sampled function
    !! Usage: ./fourier_transform.bin > data/data.txt
    !! then plot with python fourier_transform_plot.py data/data.txt
    use iso_fortran_env, only: dp => real64
	use shtawa
	implicit none
    integer, parameter :: n = 50
    real(dp), parameter :: pi = 3.141592653589793
    real(dp), dimension(1:n) :: x
    complex(dp), dimension(1:n) :: sampled_function, fourier, invfourier
    integer :: i
    
    x = [( i*2*pi/n, i=0,n-1 )]
    
    !! Case 1: simple sinus function
    ! sampled_function = sin(x)
    ! fourier = discrete_fourier_transform(sampled_function)

    !! Case 2: harmonics
    ! sampled_function = sin(x) + 0.3*sin(4*x)
    ! fourier = discrete_fourier_transform(sampled_function)
    
    !! Case 3: sum of sinus functions
    sampled_function = sin(x) + 0.3*sin(x/2) + sin(2*x-0.5)
    fourier = discrete_fourier_transform(sampled_function)

    !! Case 4: sawtooth function
    ! sampled_function = x
    ! fourier = discrete_fourier_transform(sampled_function)

    !! Testing the direct and inverse transform
    !! (we should recover the sampled function)
    invfourier = discrete_inverse_fourier_transform(fourier)

    do i=1,n
        write(*,*) fourier(i)%re,fourier(i)%im, sampled_function(i)%re, invfourier(i)%re
    end do


end program