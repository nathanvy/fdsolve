# Makefile pour fdsolve
# Auteur Nathan Van Ymeren
# https://github.com/nathanvy
CC=g++
CFLAGS=-c -Wall -std=c++11 -g3

all: fdsolve

fdsolve: main.o pointgrille.o quadrillage.o readwrite.o
	$(CC) pointgrille.o quadrillage.o readwrite.o main.o -o fdsolve

main.o:
	$(CC) $(CFLAGS) main.cpp

pointgrille.o:
	$(CC) $(CFLAGS) pointgrille.cpp

quadrillage.o:
	$(CC) $(CFLAGS) quadrillage.cpp

readwrite.o:
	$(CC) $(CFLAGS) readwrite.cpp

clean:
	rm -rf *o fdsolve
