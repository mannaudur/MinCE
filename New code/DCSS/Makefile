CC=clang++
CFLAGS=-Wall -Ofast -std=c++11
BIN=.
SRC=.
CDEPS=/usr/local/lib/libbifrost.a -lz -pthread

all: $(BIN) $(BIN)/testFastX

$(BIN)/testFastX: $(SRC)/fastXparserTest.cpp
	$(CC) $(CFLAGS) $(SRC)/fastXparserTest.cpp -o $(BIN)/testFastX $(CDEPS)

clean:
	rm -rf $(BIN)/testFastX
