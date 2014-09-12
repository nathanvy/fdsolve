# Makefile pour fdsolve
# Auteur Nathan Van Ymeren
# https://github.com/nathanvy
CC=g++
CFLAGS=-c -Wall

all: fdsolve

fdsolve: main.o pointgrille.o
	$(CC) pointgrille.o main.o -o fdsolve

main.o:
	$(CC) $(CFLAGS) main.cpp

pointgrille.o:
	$(CC) $(CFLAGS) pointgrille.cpp

clean:
	rm -rf *o fdsolve
