CFLAGS=-03 -std=c99
CC=gcc

default: conv fntt fntt_time

conv: ntt.c ntt.h
	$(CC) ntt.c -o conv -DCONV_TEST

fntt: ntt.c ntt.h
	$(CC) ntt.c -o fntt -DFNTT_TEST

fntt_time: ntt.c ntt.h
	$(CC) ntt.c -o fntt_time -DFNTT_TIME

clean:
	rm conv fntt fntt_time
