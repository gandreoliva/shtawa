module shtawa
	use iso_fortran_env, only: wp => real64
	!! References: Burden & Faires's Numerical Analysis, 9th ed, Cengage Learning.
	implicit none
	include "interp/type_grid.f90"
contains
	include "roots/bisection.f90"
	include "roots/newton.f90"
	include "roots/complex_newton.f90"
	include "roots/secant.f90"

	include "linalg/gauss.f90"
	include "linalg/gaussjordan.f90"
	include "linalg/determinant.f90"
	include "linalg/inverse.f90"
	include "linalg/tridiag_crout.f90"
	include "linalg/complex_eigenval_power.f90"

	include "combinatorics/ncartprod.f90"
	include "combinatorics/permutations.f90"
	include "combinatorics/combinations.f90"
	include "combinatorics/quicksort.f90"

	include "interp/linear_interp.f90"
	include "interp/multilinear_interp.f90"
	include "interp/spline_interp.f90"

	include "ode/rk4.f90"


end module