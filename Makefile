PROG= main
F90 = gfortran
#F90 = ifort
FLAGS= -g -O3 #-fopenmp 
MODS= mdle_io_utils.f90 mdle_adds.f90 


all:
	$(F90) -c $(MODS) $(FLAGS)
	$(F90) -c $(PROG).f90 $(FLAGS) 
	$(F90) *.o -o $(PROG) $(FLAGS) 
	rm *.o *.mod
	#qsub run.sh
clean:
	rm $(PROG) *.o *.mod *~ *.dir


