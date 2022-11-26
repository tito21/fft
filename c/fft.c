#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "utils.c"

# define PI 3.14159265358979323846264338327950288419716939937510

typedef double complex *(complex_func)(double complex *, size_t);

enum MODE { BCK, RUN };

void benchmark(int runs, int num,  complex_func func)
{
    double complex *x = random_array(num);
    double *times = (double *) malloc(sizeof(double)*runs);

    for (int i=0; i < runs; i++) {
        clock_t start = clock();
        double complex *out = func(x, num);
        clock_t stop = clock();
        free(out);
        times[i] = ((double)(stop - start))/CLOCKS_PER_SEC;
    }

    printf("Run took: %f ± %f s\n", mean(times, runs), std(times, runs));

    free(x);
    free(times);
}

double complex *dft(double complex *x, size_t N) {

    double complex *X = (double complex *) malloc(sizeof(double complex)*N);
    double complex w = cexp(-2.0 * I * PI / (double)N);

    for (int j=0; j < N; j++) {
        double complex ss = 0.0 + I*0.0;
        for (int k=0; k < N; k++) {
            ss += x[k] * cpow(w, k*j);
        }
        X[j] = ss;
    }

    return X;
}

double complex *fft(double complex *x, size_t N) {

    if (!((N & (N - 1)) == 0)) {
        printf("Number of elements must be a power of 2");
    }
    else if (N < 16) {
        return dft(x, N);
    }
    else {

        double complex *x_even = get_even(x, N);
        double complex *X_even = fft(x_even, N/2);
        double complex *x_odd = get_odd(x, N);
        double complex *X_odd = fft(x_odd, N/2);

        double complex *factor = (double complex *) malloc(sizeof(double complex)*N);

        for (int n=0; n < N; n++) {
            factor[n] = cexp(-2.0 * PI * I * (double)  n / (double) N);
        }

        double complex *s1 = slice(factor, 0, (N / 2) - 1);
        double complex *sum1 = mult(s1, X_odd, N/2);
        double complex *first = add(X_even, sum1, N/2);

        double complex *s2 = slice(factor, N / 2, N - 1);
        double complex *sum2 = mult(s2, X_odd, N/2);
        double complex *last = add(X_even, sum2, N/2);

        double complex *out = concatenate(first, N/2, last, N/2);

        free(x_even);
        free(x_odd);
        free(s1);
        free(s2);
        free(sum1);
        free(sum2);
        free(last);
        free(first);
        free(factor);
        free(X_even);
        free(X_odd);

        return out;
    }
}

int main(int argc, char** argv) {

    enum MODE mode = RUN;
    int num, reps;
    complex_func *func;

    if (argc > 1) {

        char *arg = argv[1];
        if (!strcmp(arg, "FFT")) {
            func = &fft;
        }
        else if (!strcmp(arg, "DFT")) {
            func = &dft;

        }
        num = atoi(argv[2]);
    }

    if (argc > 3) {
        char *arg = argv[3];
        if (!strcmp(arg, "benchmark")) mode = BCK;
        reps = atoi(argv[4]);
    }

    switch (mode)
    {
    case BCK:
        printf("Benchmarking\n");
        benchmark(reps, num, func);
        break;

    case RUN:
        double complex *x = random_array(num);
        double complex * out = func(x, num);
        free(x);
        free(out);
    }

    // double complex *x = random_array(num);
    // double complex *X_dft = dft(x, num);
    // double complex *X_fft = fft(x, num);

    // printf("Input\n");
    // print_complex_array(x, num);
    // printf("DFT\n");
    // print_complex_array(X_dft, num);
    // printf("FFT\n");
    // print_complex_array(X_fft, num);

    // free(x);
    // free(X_dft);
    // free(X_fft);

    return 0;
}