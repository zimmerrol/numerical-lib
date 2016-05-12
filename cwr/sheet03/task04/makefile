CC=g++
SOURCE=task04.cpp ../../../lib/diff.cpp
CFLAGS=-Wall
PROGRAM=task04

OBJECTS=$(SOURCE:.cpp=.o)

all: $(PROGRAM)
build:
	$(CC) -o $(PROGRAM) $(SOURCE)
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -f $(PROGRAM) *.o
