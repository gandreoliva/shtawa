module shtawa
    use iso_fortran_env, only: dp => real64
    !! Reference for algorithms: Burden & Faires's Numerical Analysis, 9th ed, Cengage Learning.
    implicit none
contains
    include "single_var_eqs/bisection.f90"
    include "single_var_eqs/newton.f90"
end module