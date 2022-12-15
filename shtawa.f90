module shtawa
    use iso_fortran_env, only: dp => real64
    !! Reference for algorithms: Burden & Faires's Numerical Analysis, 9th ed, Cengage Learning.
    implicit none
contains
    include "roots/bisection.f90"
    include "roots/newton.f90"
    include "roots/complex_newton.f90"
    include "roots/secant.f90"
end module