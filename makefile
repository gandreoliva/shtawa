MAKEFLAGS += --no-builtin-rules --no-builtin-variables

all:
	mkdir -p bin
	gfortran -c -o bin/shtawa.o -J bin/ shtawa.f90

clean:
	rm bin/*