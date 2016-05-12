CC=gcc
SOURCE=task02.c
CFLAGS=-Wall
PROGRAM=task02

all: $(PROGRAM)

build:
	$(CC) -o $(PROGRAM) $(SOURCE)
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
