CC = g++
CFLAGS = -Wall -O3 --std=c++11 -g
U_F = ../general_code/
#useful folder
PS_F = ../particle_simulator/

EX = ex_4

all : $(EX).o $(U_F)random.o $(U_F)useful_functions.o $(PS_F)particle_simulator.o
	$(CC) $(EX).o $(U_F)random.o $(U_F)useful_functions.o $(PS_F)particle_simulator.o -o $(EX).exe
$(EX).o : $(EX).cpp
	$(CC) -c $(EX).cpp -o $(EX).o $(CFLAGS)
random.o : $(U_F)random.cpp $(U_F)random.h
	$(CC) -c $(U_F)random.cpp -o $(U_F)random.o $(CFLAGS)
useful_functions.o : $(U_F)useful_functions.cpp $(U_F)standard_header.h
	$(CC) -c $(U_F)useful_functions.cpp -o $(U_F)useful_functions.o $(CFLAGS)
particle_simulator.o: $(PS_F)particle_simulator.cpp $(PS_F)particle_simulator.h
	$(CC) -c $(PS_F)particle_simulator.cpp -o $(PS_F)particle_simulator.o $(CFLAGS)
clean:
	rm *.o
	rm *.exe
