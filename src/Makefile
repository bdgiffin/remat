CC=g++
CFLAGS=-std=c++17

all: test_Rational test_Fixed test_Dual

.PHONY: clean test_Rational test_Fixed test_Dual

test_Rational: test_Rational.cpp Rational.h types.h
	$(CC) $(CFLAGS) test_Rational.cpp -o test_Rational

test_Fixed: test_Fixed.cpp Fixed.h types.h
	$(CC) $(CFLAGS) test_Fixed.cpp -o test_Fixed

test_Dual: test_Dual.cpp Dual.h Fixed.h types.h
	$(CC) $(CFLAGS) test_Dual.cpp -o test_Dual

clean:
	rm -rf test_Rational test_Fixed test_Dual
