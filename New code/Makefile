CC=clang++
CFLAGS=-Wall -Ofast -std=c++11
BIN=local/bin
SRC=src
CDEPS=/usr/local/lib/libbifrost.a -lz -pthread

all: $(BIN) $(BIN)/mince $(BIN)/sketch $(BIN)/find_cliques $(BIN)/create_bitmatrix $(BIN)/cliqdist

$(BIN)/mince: $(SRC)/Mince.cpp
	$(CC) $(CFLAGS) $(SRC)/Mince.cpp -o $(BIN)/mince $(CDEPS)

$(BIN)/sketch: $(SRC)/Sketch.cpp $(SRC)/Sketch.hpp
	$(CC) $(CFLAGS) $(SRC)/Sketch.cpp -o $(BIN)/sketch $(CDEPS)

$(BIN)/find_cliques: $(SRC)/CliquesFromHL.cpp
	$(CC) $(CFLAGS) $(SRC)/CliquesFromHL.cpp -o $(BIN)/find_cliques $(CDEPS)

$(BIN)/create_bitmatrix: $(SRC)/bitmats.cpp
	$(CC) $(CFLAGS) $(SRC)/bitmats.cpp -o $(BIN)/create_bitmatrix $(CDEPS)

$(BIN)/cliqdist: $(SRC)/distance_within_cliques.cpp
	$(CC) $(CFLAGS) $(SRC)/distance_within_cliques.cpp -o $(BIN)/cliqdist $(CDEPS)

clean:
	rm -rf $(BIN)/sketch $(BIN)/mince $(BIN)/find_cliques $(BIN)/create_bitmatrix $(BIN)/cliqdist

	rm -f *.sketch fastx/*.xsketch fastx/*.msh fastx/*.sketch
