CC=clang++
CFLAGS=-Wall -Ofast -std=c++11
BIN=local/bin
SRC=src
CDEPS=/usr/local/lib/libbifrost.a -lz -pthread

all: $(BIN) $(BIN)/mince $(BIN)/sketch $(BIN)/cluster $(BIN)/distinguish $(BIN)/sift

$(BIN)/mince: $(SRC)/Mince.cpp
	$(CC) $(CFLAGS) $(SRC)/Mince.cpp -o $(BIN)/mince $(CDEPS)

$(BIN)/sketch: $(SRC)/Sketch.cpp $(SRC)/Sketch.hpp
	$(CC) $(CFLAGS) $(SRC)/Sketch.cpp -o $(BIN)/sketch $(CDEPS)

$(BIN)/cluster: $(SRC)/UFClusters.cpp
	$(CC) $(CFLAGS) $(SRC)/UFClusters.cpp -o $(BIN)/cluster $(CDEPS)

$(BIN)/distinguish: $(SRC)/Distinguish.cpp
	$(CC) $(CFLAGS) $(SRC)/Distinguish.cpp -o $(BIN)/distinguish $(CDEPS)

$(BIN)/sift: $(SRC)/KmerSieve.cpp
	$(CC) $(CFLAGS) $(SRC)/KmerSieve.cpp -o $(BIN)/sift $(CDEPS)

clean:
	rm -rf $(BIN)/sketch $(BIN)/mince $(BIN)/cluster $(BIN)/distinguish $(BIN)/sift

	rm -f *.sketch fastx/*.xsketch fastx/*.msh fastx/*.sketch
