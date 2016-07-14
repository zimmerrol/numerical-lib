CC=g++
SOURCE=pendulum.cpp lib/ode.cpp
CFLAGS=-Wall
PROGRAM=bin/pendulum

SOURCE2=bifurcation.cpp lib/ode.cpp
PROGRAM2=bin/bifurcation

OBJECTS=$(SOURCE:.cpp=.o)

all: $(PROGRAM)
build:
	#create directory structure
	mkdir -p results
	mkdir -p results/bifurcation
	mkdir -p results/phase_diagrams

	$(CC) -o $(PROGRAM) $(SOURCE)
	$(CC) -o $(PROGRAM2) $(SOURCE2)
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -f $(PROGRAM) *.o
	rm -f $(PROGRAM2) *.o
