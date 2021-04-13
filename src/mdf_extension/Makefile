
# Makefile for the mdf extension of the parser

CC = gcc
# Adder no-maybe-uninitialized, as it doesn't help debugging
CFLAGS = -std=gnu99 -march=native -O3 -Wall -Wno-maybe-uninitialized -D_SHARED #-D_UNROLL

LIBRARIES = -lm -lz

EXECNAME = qLD-parse-2MDF

SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c, build/%.o, $(SRCS))

all: $(EXECNAME)

qLD-parse-2MDF : $(OBJS)
	$(CC) $(CFLAGS) -o $(EXECNAME) $(OBJS) $(LIBRARIES)

build/%.o: src/%.c
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm $(EXECNAME)
	rm $(OBJS)
