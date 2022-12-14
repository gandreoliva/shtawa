MAKEFLAGS += --no-builtin-rules --no-builtin-variables

bin := ../../bin

all:
	echo "No target chosen"

bisection_example:
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

newton_example:
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

newton_fractal:
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

clean:
	rm *.bin

clean_data:
	rm data/*