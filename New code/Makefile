CC=clang++
CFLAGS=-Wall -Ofast -std=c++11
BIN=local/bin
SRC=src
CDEPS=/usr/local/lib/libbifrost.a -lz -pthread

all: $(BIN) $(BIN)/mince $(BIN)/sketch $(BIN)/cluster $(BIN)/sequences $(BIN)/trialSequences $(BIN)/hashSequences

$(BIN)/mince: $(SRC)/Mince.cpp
	$(CC) $(CFLAGS) $(SRC)/Mince.cpp -o $(BIN)/mince $(CDEPS)

$(BIN)/sketch: $(SRC)/Sketch.cpp $(SRC)/Sketch.hpp
	$(CC) $(CFLAGS) $(SRC)/Sketch.cpp -o $(BIN)/sketch $(CDEPS)

$(BIN)/cluster: $(SRC)/UFClusters.cpp
	$(CC) $(CFLAGS) $(SRC)/UFClusters.cpp -o $(BIN)/cluster $(CDEPS)

$(BIN)/trialSequences: $(SRC)/TrialSequences.cpp
	$(CC) $(CFLAGS) $(SRC)/TrialSequences.cpp -o $(BIN)/trialSequences $(CDEPS)

$(BIN)/hashSequences: $(SRC)/makeSeqHashmap.cpp
	$(CC) $(CFLAGS) $(SRC)/makeSeqHashmap.cpp -o $(BIN)/hashSequences $(CDEPS)

clean:
	rm -rf $(BIN)/sketch $(BIN)/mince $(BIN)/cluster $(BIN)/sequences $(BIN)/trialSequences $(BIN)/hashSequences

	rm -f *.sketch fastx/*.xsketch fastx/*.msh fastx/*.sketch
