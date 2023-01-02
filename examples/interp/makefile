MAKEFLAGS += --no-builtin-rules --no-builtin-variables

bin := ../../bin

all:
	echo "No target chosen"

linear_interp_example: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

multilinear_interp_example: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

multilinear_interp_example2: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

shtawa:
	make -C ../../

clean:
	rm *.bin

clean_data:
	rm data/*