CC=gcc
CFLAGS=-std=c11 -lm -Wall -Werror -g

OBJS = fft.o array.o

all: fft

check-bck: fft
	valgrind --leak-check=full \
		--show-leak-kinds=all \
		--track-origins=yes \
		--verbose \
		--log-file=valgrind-out-benchmark.txt \
		./fft FFT 128 benchmark 100

check-run: fft
	valgrind --leak-check=full \
		--show-leak-kinds=all \
		--track-origins=yes \
		--verbose \
		--log-file=valgrind-out-run.txt \
		./fft FFT 128

fft: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o fft

fft.o: fft.c
	$(CC) $(CFLAGS) -c fft.c

array.o: array.c array.h
	$(CC) $(CFLAGS) -c array.c

clean:
	rm -f *~ *.o fft