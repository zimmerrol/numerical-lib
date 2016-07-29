CC=g++
SOURCE=task03.cpp ../../../lib/ode.cpp
CFLAGS=-Wall
PROGRAM=task03

OBJECTS=$(SOURCE:.cpp=.o)

all: $(PROGRAM)
build:
	$(CC) -o $(PROGRAM) $(SOURCE)
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -f $(PROGRAM) *.o
