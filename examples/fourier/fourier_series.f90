program fourier_series_example
    !! Computes the coefficients of a Fourier series expansion of a discrete, periodic, sampled function
    !! Usage: ./fourier_series > data/data.txt
    !! Then, plot with python fourier_series_plot.py data/data.txt
	use iso_fortran_env, only: dp => real64
	use shtawa
	implicit none
    !! Try n = 100 (same as in plotting) and another number
    integer, parameter :: n = 100
    real(dp), parameter :: pi = 3.141592653589793
    real(dp), dimension(-n/2:n/2) :: sampled_function, x
    complex(dp), dimension(-n/2:n/2) :: coeff
    complex(dp), parameter :: iimag = (0d0,1d0)
    integer :: i
    real(dp), dimension(1:n/2) :: a,b
    real(dp) :: a0
    
    do i=-n/2,n/2
        x(i) = i*2*pi/n
    end do

    !! Test function 1: cosine
    ! sampled_function = cos(x)
    ! coeff = discrete_fourier_coeff(sampled_function)

    !! Test function 2: sum of sinus
    ! sampled_function = sin(x) + 0.3*sin(x/2) + sin(2*x-0.5)
    ! coeff = discrete_fourier_coeff(sampled_function)

    !! Test function 3: sawtooth
    ! sampled_function = x
    ! coeff = discrete_fourier_coeff(sampled_function)

    !! Test function 4: harmonics
    sampled_function = sin(x) + 0.3*sin(4*x)
    coeff = discrete_fourier_coeff(sampled_function)

    !! Transformation to real coefficients
    a0 = 2*real(coeff(0))

    write(*,*) a0, 0d0
    
    do i=1,n/2
        a(i) = real(coeff(i) + coeff(-i))
        b(i) = real(iimag*(coeff(i)-coeff(-i)))
        write(*,*) a(i), b(i)
    end do

    write(*,*) "---"

    do i=-n/2,n/2
        write(*,*) sampled_function(i)
    end do

end program