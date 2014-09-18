# Makefile pour fdsolve
# Auteur Nathan Van Ymeren
# https://github.com/nathanvy
CC=g++
CFLAGS=-c -std=c++11
RELEASE_FLAGS=-O3
DEBUG_FLAGS=-Wall -g3

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

release: CFLAGS += $(RELEASE_FLAGS)
release: all

debug: CFLAGS += $(DEBUG_FLAGS)
debug: all

clean:
	rm -rf *o fdsolve
