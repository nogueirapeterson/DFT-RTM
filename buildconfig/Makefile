PROG = main
F90 = nvfortran
#F90 = ifort
FLAGS= -g -O3 -acc -Minfo=accel -ta=multicore #-fopenmp
MAIN= ../app/main.f90
MODS= ../src/mdle_io_utils.f90 ../src/mdle_adds.f90
PYPLOT= ../external/pyplot-fortran/src/pyplot_module.F90
TEST=../test/test.f90
TEST1=../test/test_atenuation.f90


all:
	$(F90) -c $(MODS) $(FLAGS)
	$(F90) -c $(MAIN) $(FLAGS)
	$(F90) *.o -o $(PROG) $(FLAGS)
	rm *.o *.mod
	#qsub run.sh
clean:
	rm $(PROG) *.o *.mod *~ *.dir

test: $(MODS) $(TEST)
	$(F90) -c $(PYPLOT) $(MODS) $(FLAGS)
	$(F90) -c $(PYPLOT) $(FLAGS) -cpp -DREAL32
	$(F90) -c $(TEST) $(FLAGS)
	$(F90) *.o -o test $(FLAGS)
	$(F90) *.o -o test-att $(FLAGS)
	rm *.o *.mod

test1: $(MODS) $(TEST)
	$(F90) -c $(PYPLOT) $(MODS) $(FLAGS)
	$(F90) -c $(PYPLOT) $(FLAGS) -cpp -DREAL32
	$(F90) -c $(TEST1) $(FLAGS)
	$(F90) *.o -o test-att $(FLAGS)
	rm *.o *.mod




