CC=clang++
CFLAGS=-Wall -Ofast -std=c++11
BIN=local/bin
SRC=src
CDEPS=/usr/local/lib/libbifrost.a -lz -pthread

all: $(BIN) $(BIN)/mince $(BIN)/sketch $(BIN)/cluster $(BIN)/sequences $(BIN)/hashSequences $(BIN)/inspect

$(BIN)/mince: $(SRC)/Mince.cpp
	$(CC) $(CFLAGS) $(SRC)/Mince.cpp -o $(BIN)/mince $(CDEPS)

$(BIN)/sketch: $(SRC)/Sketch.cpp $(SRC)/Sketch.hpp
	$(CC) $(CFLAGS) $(SRC)/Sketch.cpp -o $(BIN)/sketch $(CDEPS)

$(BIN)/cluster: $(SRC)/UFClusters.cpp
	$(CC) $(CFLAGS) $(SRC)/UFClusters.cpp -o $(BIN)/cluster $(CDEPS)

$(BIN)/sequences: $(SRC)/Sequences.cpp
	$(CC) $(CFLAGS) $(SRC)/Sequences.cpp -o $(BIN)/sequences $(CDEPS)

$(BIN)/hashSequences: $(SRC)/makeSeqHashmap.cpp
	$(CC) $(CFLAGS) $(SRC)/makeSeqHashmap.cpp -o $(BIN)/hashSequences $(CDEPS)

$(BIN)/inspect: $(SRC)/Inspect.cpp
	$(CC) $(CFLAGS) $(SRC)/Inspect.cpp -o $(BIN)/inspect $(CDEPS)

clean:
	rm -rf $(BIN)/mince $(BIN)/sketch  $(BIN)/cluster $(BIN)/sequences $(BIN)/hashSequences $(BIN)/inspect

	rm -f *.sketch fastx/*.xsketch fastx/*.msh fastx/*.sketch
