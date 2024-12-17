MAKEFLAGS += --no-builtin-rules --no-builtin-variables

bin := ../../bin

all:
	echo "No target chosen"

fourier_series: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

fourier_transform: shtawa
	gfortran -I $(bin) $(bin)/shtawa.o $@.f90 -o $@.bin

shtawa:
	make -C ../../

clean:
	rm *.bin

clean_data:
	rm data/*