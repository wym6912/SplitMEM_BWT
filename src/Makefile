SDSLLITE = ../include/sdsl
THREADS  = 1

CC       = g++
CC_FLAGS = -O9 -std=c++17 -DNDEBUG -ffast-math -funroll-loops -msse4.2 -Wall -I${SDSLLITE}/include -L${SDSLLITE}/lib -static
CCLIB    = -lsdsl -ldivsufsort -ldivsufsort64 -lpardivsufsort -lprange -Wl,--whole-archive -ldl -Wl,--no-whole-archive -fopenmp
SOURCES  =$(wildcard *.cpp)
EXECS    =$(SOURCES:.cpp=)
TARGETS  = SplitMEM

all: sdsl $(TARGETS)

sdsl:
	make -C ../include/sdsl -j$(THREADS)

SplitMEM: SplitMEM.cpp *.hpp
	$(CC) $(CC_FLAGS) $< -o $@ $(CCLIB)

clean:
	rm -f $(TARGETS)
	make -C ../include/sdsl clean

