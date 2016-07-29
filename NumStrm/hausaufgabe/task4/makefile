CC=g++
SOURCE=task04.cpp ../../../lib/linear_system.cpp ../../../lib/matrix.cpp ../../../lib/square_matrix.cpp 
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
