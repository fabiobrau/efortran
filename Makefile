FC = gfortran

FCFLAGS = -O0 -fbounds-check 
FCFLAGSDB = -g -O0 -Wall -Wextra -pedantic -fbounds-check -fimplicit-none -fbacktrace 
PROGRAMS =main.out 

SOURCE = src

LAPACK = -llapack -lblas

all: $(PROGRAMS)

main.out:
	$(FC) $(FCFLAGS) -c $(SOURCE)/easy_matrix.f08 -o obj/easy_matrix.o $(LAPACK)
	$(FC) $(FCFLAGS) -c $(SOURCE)/iter.f08 -o obj/iter.o $(LAPACK)
	$(FC) $(FCFLAGS) -c $(SOURCE)/main.f08 -o obj/main.o $(LAPACK)

	$(FC) $(FCFLAGS) obj/*.o -o exec/main.out $(LAPACK)
clean:
	-(rm obj/*)
	-(rm *.mod)
	-(rm exec/*.out)
