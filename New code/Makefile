CC=clang++
CFLAGS=-Wall -Ofast -std=c++11
BIN=local/bin
SRC=src
CDEPS=/usr/local/lib/libbifrost.a -lz -pthread

all: $(BIN) $(BIN)/mince $(BIN)/sketch $(BIN)/find_cliques $(BIN)/create_bitmatrix

$(BIN)/mince: $(SRC)/Mince.cpp
	$(CC) $(CFLAGS) $(SRC)/Mince.cpp -o $(BIN)/mince $(CDEPS)

$(BIN)/sketch: $(SRC)/Sketch.cpp $(SRC)/Sketch.hpp
	$(CC) $(CFLAGS) $(SRC)/Sketch.cpp -o $(BIN)/sketch $(CDEPS)

$(BIN)/find_cliques: $(SRC)/FindCliques.cpp
	$(CC) $(CFLAGS) $(SRC)/FindCliques.cpp -o $(BIN)/find_cliques $(CDEPS)

$(BIN)/create_bitmatrix: $(SRC)/CreateBitmatrix.cpp
	$(CC) $(CFLAGS) $(SRC)/CreateBitmatrix.cpp -o $(BIN)/create_bitmatrix $(CDEPS)

clean:
	rm -rf $(BIN)/sketch $(BIN)/mince $(BIN)/find_cliques $(BIN)/create_bitmatrix

	rm -f *.sketch fastx/*.xsketch fastx/*.msh fastx/*.sketch
