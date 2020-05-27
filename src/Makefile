CFLAGS=-03 -std=c99
CC=gcc

default: conv fntt

conv: ntt.c ntt.h
	$(CC) ntt.c -o conv -DCONV_TEST

fntt: ntt.c ntt.h
	$(CC) ntt.c -o fntt -DFNTT_TEST

clean:
	rm conv fntt
