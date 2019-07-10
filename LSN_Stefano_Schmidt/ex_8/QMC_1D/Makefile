%.o : %.cpp
	g++ -Wall -c $< 

qmc1d: qmc1d.o
	g++ -O3 -Wall -o $@ $^ 

clean:
	rm *.o qmc1d potential.dat kinetic.dat probability.dat
