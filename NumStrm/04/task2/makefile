CC=g++
SOURCE=task02.cpp ../../../lib/ode.cpp ../../../lib/linear_system.cpp ../../../lib/square_matrix.cpp ../../../lib/matrix.cpp
CFLAGS=-Wall
PROGRAM=task02

OBJECTS=$(SOURCE:.cpp=.o)

all: $(PROGRAM)
build:
	$(CC) -o $(PROGRAM) $(SOURCE)
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -f $(PROGRAM) *.o
