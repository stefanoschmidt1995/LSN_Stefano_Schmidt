CC = g++
CFLAGS = -Wall -O3 --std=c++11 -g
U_S = ../general_code/

EX = ex_3

all : $(EX).o $(U_S)random.o $(U_S)useful_functions.o
	$(CC) $(EX).o $(U_S)random.o $(U_S)useful_functions.o -o $(EX).exe
$(EX).o : $(EX).cpp
	$(CC) -c $(EX).cpp -o $(EX).o $(CFLAGS)
random.o : $(U_S)random.cpp $(U_S)random.h
	$(CC) -c $(U_S)random.cpp -o $(U_S)random.o $(CFLAGS)
useful_functions.o : $(U_S)useful_functions.cpp $(U_S)standard_header.h
	$(CC) -c $(U_S)useful_functions.cpp -o $(U_S)useful_functions.o $(CFLAGS)
clean:
	rm *.o
	rm *.exe
