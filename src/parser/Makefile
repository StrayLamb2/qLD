
# Makefile for the parser

CC = gcc
# Adder no-maybe-uninitialized, as it doesn't help debugging
# MPAMPIS
CFLAGS = -O3 -Wall -Wno-maybe-uninitialized -D_SHARED -g#-D_UNROLL

LIBRARIES = -lm -lz

EXECNAME = qLD-parse-VCF

SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c, build/%.o, $(SRCS))

all: $(EXECNAME)

qLD-parse-VCF : $(OBJS)
	$(CC) $(CFLAGS) -o $(EXECNAME) $(OBJS) $(LIBRARIES)

build/%.o: src/%.c
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm $(EXECNAME)
	rm $(OBJS)
