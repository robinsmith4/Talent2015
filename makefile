XLF = gfortran 

all_objects = basis11_Robin2.o 

prog.exe : ${all_objects}
	${XLF} -o prog.exe ${all_objects} 

%.o: %.F90 
	${XLF} -c $<

clean:
	rm *.mod *.exe *.o


