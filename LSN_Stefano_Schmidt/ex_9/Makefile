CC = g++
CFLAGS = -Wall -O3 --std=c++11 -g
GA_F = ../GA_helper/
U_F = ../general_code/

EX = ex_9

all : $(EX).o $(U_F)random.o $(U_F)useful_functions.o $(GA_F)GA_helpers.o
	$(CC) $(EX).o $(U_F)random.o $(U_F)useful_functions.o  $(GA_F)GA_helpers.o -o $(EX).exe
$(EX).o : $(EX).cpp ../standard_header.h
	$(CC) -c $(EX).cpp -o $(EX).o $(CFLAGS)
GA_helpers.o: $(GA_F)GA_helpers.h $(GA_F)GA_helpers.cpp
	$(CC) -c $(GA_F)GA_helpers.cpp -o $(GA_F)GA_helpers.o $(CFLAGS)
random.o : $(U_F)random.cpp $(U_F)random.h
	$(CC) -c $(U_F)random.cpp -o $(U_F)random.o $(CFLAGS)
useful_functions.o : $(U_F)useful_functions.cpp standard_header.h
	$(CC) -c $(U_F)useful_functions.cpp -o $(U_F)useful_functions.o $(CFLAGS)
clean:
	rm *.o
	rm *.exe
