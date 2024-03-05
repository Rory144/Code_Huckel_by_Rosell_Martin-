
compile: 
	gfortran -c mod.f90 -o mod.o
	gfortran -c main_program.f90 -o main_program.o -I./
	gfortran -o main_program.exe main_program.o mod.o -llapack 
execute: 
	./main_program.exe

clean: 
	@rm *.o
	@rm *.mod

