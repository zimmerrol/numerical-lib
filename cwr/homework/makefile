CC=g++
SOURCE=task.cpp ../../lib/ode.cpp
CFLAGS=-Wall
PROGRAM=task

SOURCE2=bifurcation.cpp ../../lib/ode.cpp
PROGRAM2=bifurcation

OBJECTS=$(SOURCE:.cpp=.o)

all: $(PROGRAM)
build:
	$(CC) -o $(PROGRAM) $(SOURCE)
	$(CC) -o $(PROGRAM2) $(SOURCE2)
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -f $(PROGRAM) *.o
