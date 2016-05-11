CC=g++
SOURCE=task03.cpp ../../lib/root.cpp
CFLAGS=-Wall -std=c++11
PROGRAM=task03

OBJECTS=$(SOURCE:.cpp=.o)

all: $(PROGRAM)
build:
	$(CC) $(CFLAGS) -o $(PROGRAM) $(SOURCE)
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -f $(PROGRAM) *.o
