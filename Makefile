objects = angular_hod.o
LINKS = -lm -lgsl -lgslcblas --std=c99
run:$(objects)
	mpicc -o run $^ $(LINKS) 
angular_hod.o:angular_hod.c
	mpicc -c $^ -O3  
clean:
	rm run $(objects)
