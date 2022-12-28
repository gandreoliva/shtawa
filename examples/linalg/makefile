MAKEFLAGS += --no-builtin-rules --no-builtin-variables

bin := ../../bin

all:
	echo "No target chosen"

gauss_example: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

gaussjordan_example: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

determinant_example: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin
	
inverse_example: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

tridiag_crout_example: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

complex_eigenval_power_example: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

shtawa:
	make -C ../../

clean:
	rm *.bin

clean_data:
	rm data/*