program rk4_example
    use iso_fortran_env, only: dp=> real64
    use shtawa
    implicit none
    
    call rk4(f,y0=[0.3d0],t_initial=0d0,t_final=2d0,niter=100)

contains
    function f(t,y)
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(in) :: y
        real(dp) :: f
        f = -0.4*t
    end function

end program