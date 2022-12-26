MAKEFLAGS += --no-builtin-rules --no-builtin-variables

bin := ../../bin

all:
	echo "No target chosen"

bisection_example: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

newton_example: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

newton_fractal: shtawa
	mkdir -p data
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

secant_example: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

shtawa:
	make -C ../../

clean:
	rm *.bin

clean_data:
	rm data/*