CC=g++
SOURCE=task01.cpp ../../../lib/ode.cpp fft_fn.c
CFLAGS=-Wall
PROGRAM=task01

OBJECTS=$(SOURCE:.cpp=.o)

all: $(PROGRAM)
build:
	$(CC) -o $(PROGRAM) $(SOURCE)
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -f $(PROGRAM) *.o
