CC=gcc
CFLAGS=-std=c11 -lm

build:
	$(CC) fft.c -o fft -O3 $(CFLAGS)

fft-debug: fft.c
	$(CC) fft.c -o fft-debug -pg $(CFLAGS)

profile: fft-debug
	./fft-debug FFT 1024
	gprof fft-debug gmon.out > analysis.txt

check-run: fft-debug
	valgrind --leak-check=full \
		--show-leak-kinds=all \
		--track-origins=yes \
		--verbose \
		--log-file=valgrind-out-run.txt \
		./fft-debug FFT 64

check-bck: fft-debug
	valgrind --leak-check=full \
		--show-leak-kinds=all \
		--track-origins=yes \
		--verbose \
		--log-file=valgrind-out-bck.txt \
		./fft-debug FFT 64 benchmark 100

clean:
	rm fft fft-debug gmon.out analysis.txt valgrind-out-run.txt valgrind-out-bck.txt