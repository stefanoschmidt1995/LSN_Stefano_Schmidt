CC = g++
CCPAR = mpicxx
CFLAGS = -Wall -O3 --std=c++11 -g
U_F = ../general_code/
GA_F = ../GA_helper/

EX = ex_10
PAR = ex_10_par

all : $(EX).o $(U_F)random.o $(U_F)useful_functions.o  GA_helpers.o
	$(CC) $(EX).o $(U_F)random.o $(U_F)useful_functions.o $(GA_F)GA_helpers.o -o $(EX).exe
$(EX).o : $(EX).cpp $(GA_F)GA_helpers.cpp
	$(CC) -c $(EX).cpp -o $(EX).o $(CFLAGS)
parallel : $(PAR).cpp GA_helpers.o random.o
	$(CCPAR) -c $(PAR).cpp -o $(PAR).o $(CFLAGS)
	$(CCPAR) $(PAR).o $(U_F)random.o $(GA_F)GA_helpers.o $(U_F)useful_functions.o -o $(PAR).exe
GA_helpers.o: $(GA_F)GA_helpers.h $(GA_F)GA_helpers.cpp
	$(CC) -c $(GA_F)GA_helpers.cpp -o $(GA_F)GA_helpers.o $(CFLAGS)
random.o : $(U_F)random.cpp $(U_F)random.h
	$(CC) -c $(U_F)random.cpp -o $(U_F)random.o $(CFLAGS)
useful_functions.o : $(U_F)useful_functions.cpp standard_header.h
	$(CC) -c $(U_F)useful_functions.cpp -o $(U_F)useful_functions.o $(CFLAGS)
clean:
	rm *.o
	rm *.exe
