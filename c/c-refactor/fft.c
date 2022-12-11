#include <complex.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "array.h"

#define PI 3.14159265358979323846264338327950288419716939937510

enum MODE { BCK, RUN };
enum FUNC { FFT, FFT_ITR, DFT };

double complex random_complex(int n);
void benchmark(int runs, int num, enum FUNC func);
CMP_Array dft(CMP_Array x);
CMP_Array fft(CMP_Array x, int stride);
CMP_Array fft_iterative(CMP_Array x);
CMP_Array bit_reverse_array(CMP_Array x);
unsigned int reverse_bit(unsigned int v, unsigned int m);
double mean(double *arr, size_t N);
double std(double *arr, size_t N);

int main(int argc, char** argv)
{

    enum MODE mode = RUN;
    int num, reps;
    enum FUNC func = FFT;

    if (argc > 1) {

        char *arg = argv[1];
        if (!strcmp(arg, "FFT")) {
            func = FFT;
        }
        else if (!strcmp(arg, "DFT")) {
            func = DFT;

        }
        else if (!strcmp(arg, "FFT-ITR")) {
            func = FFT_ITR;

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
        printf("Running\n");
        double complex *buffer = malloc(sizeof(*buffer) * num);
        CMP_Array x = { .N = num, .data = buffer };
        fill_func(x, random_complex);
        print_arr(x);
        CMP_Array out = func == FFT ? fft(x, 1) : func == FFT_ITR ? fft_iterative(x) : dft(x);
        print_arr(out);

        CMP_Array out1 = fft(x, 1);
        print_arr(out1);


        CMP_Array out2 = dft(x);
        print_arr(out2);

        free(x.data);
        free(out.data);
        free(out1.data);
        free(out2.data);
    }

    // int num_elem = 32;
    // double complex *buffer = malloc(sizeof(*buffer) * num_elem);
    // CMP_Array arr = { .N = num_elem, .data = buffer };

    // fill_func(arr, random_complex);
    // print_arr(arr);

    // CMP_Array dft_arr = dft(arr);
    // print_arr(dft_arr);

    // CMP_Array fft_arr = fft(arr, 1);
    // print_arr(fft_arr);

    // free(arr.data);
    // free(dft_arr.data);
    // free(fft_arr.data);

    return 0;

}

double complex random_complex(int n)
{
    double num1 = (double) rand() / RAND_MAX * 2.0 - 1.0;
    double num2 = (double) rand() / RAND_MAX * 2.0 - 1.0;
    return num1 + I*num2;
}

double mean(double *arr, size_t N) {
    double sum = 0;
    for (int i=0; i < N; i++) {
        sum += arr[i];
    }
    return sum / N;
}

double std(double *arr, size_t N) {

    double m = mean(arr, N);
    double sum = 0;
    for (int i=0; i < N; i++) {
        sum += (arr[i] - m)*(arr[i] - m);
    }
    return sqrt(sum / N);
}

void benchmark(int runs, int num, enum FUNC func)
{
    double complex *buffer = malloc(sizeof(*buffer) * num);
    CMP_Array x = { .N = num, .data = buffer };

    fill_func(x, random_complex);

    double *times = (double *) malloc(sizeof(double)*runs);

    for (int i=0; i < runs; i++) {
        clock_t start = clock();
        CMP_Array out = func == FFT ? fft(x, 1) : func == FFT_ITR ? fft_iterative(x) : dft(x);
        clock_t stop = clock();
        free(out.data);
        times[i] = ((double)(stop - start))/CLOCKS_PER_SEC;
    }

    printf("Run took: %f ± %f s\n", mean(times, runs), std(times, runs));

    free(x.data);
    free(times);
}

CMP_Array dft(CMP_Array x)
{
    int num_elem = x.N;
    double complex *buffer = malloc(sizeof(*buffer) * num_elem);
    CMP_Array X = { .N = num_elem, .data = buffer };

    double complex w = cexp(-2.0 * I * PI / (double)num_elem);
    for (int j=0; j < num_elem; j++) {
        double complex ss = 0.0 + I*0.0;
        for (int k=0; k < num_elem; k++) {
            ss += x.data[k] * cpow(w, k*j);
        }
        X.data[j] = ss;
    }

    return X;
}

CMP_Array fft(CMP_Array x, int stride)
{
    int N = x.N;
    double complex *buffer = malloc(sizeof(*buffer) * N);
    CMP_Array X = { .N = N, .data = buffer };

    if (N == 1) {
        X.data[0] = x.data[0];
        return X;
    }
    else if (!((N & (N - 1)) == 0)) {
        printf("Number of elements must be a power of 2");
        free(X.data);
        return x;
    }
    else {
        CMP_Array x_even = { .N = N/2, .data = x.data };
        CMP_Array x_odd = { .N = N/2, .data = x.data+stride };

        CMP_Array X_1 = fft(x_even, 2*stride);
        CMP_Array X_2 = fft(x_odd, 2*stride);

        double complex factor = cexp(-2.0 * PI * I / (double) N);

        for (int k=0; k < N/2; k++) {

            double complex p = X_1.data[k];
            double complex q = cpow(factor, k) * X_2.data[k];

            X.data[k] = p + q;
            X.data[k+N/2] = p - q;
        }

        free(X_1.data);
        free(X_2.data);

        return X;
    }
}

CMP_Array fft_iterative(CMP_Array x)
{
    printf("Starting bit reverse\n");
    CMP_Array X = bit_reverse_array(x);
    printf("Reversed array\n");
    print_arr(X);
    int N = x.N;
    for (int s=1; s < (int)log2(N); s++) {
        int m = pow(2, s);
        // printf("%d, %d\n", s, m);
        double complex w_m = cexp(-2.0 * PI * I / (double) m);
        // printf("%f, %f\n", creal(w_m), cimag(w_m));

        for (int k=0; k < N; k+=m) {
            double complex w = 1.0 + I*0.0;
            for (int j=0; j < m/2; j++) {
                // printf("%f, %f\n", creal(w), cimag(w));
                double complex t = w * X.data[k + j + m/2];
                double complex u = X.data[k + j];
                X.data[k + j] = u + t;
                X.data[k + j + m/2] = u - t;
                w = w * w_m;
            }
        }
    }
    return X;
}

unsigned int reverse_bit(unsigned int v, unsigned int m)
{
    unsigned int r = 0;
    // unsigned int s = sizeof(v)*8 - 1;
    // unsigned int s = m - 1;

    // for (v >>= 1; v; v >>= 1) {
    //     r <<= 1;
    //     r |= v & 1;
    //     s--;
    // }
    // r <<= s;
    // printf("Num: %d, log2: %d\n", m, (int) log2(m));

    // if (v == 0) return m /2;

    unsigned int num_half = m / 2;
    unsigned int carry = 0;
    // printf("input: %u num_half: %u\n", v, num_half);

    for (int i=(int) log2(m) - 1; i > -1; i--) {
        unsigned int bit_input = v & (1 << i);
        unsigned int bit_num_half = num_half & (1 << i);
        unsigned int carry_bit = (carry << i);
        r += (bit_input + bit_num_half + carry_bit) & (1 << i);

        carry = (bit_input && bit_num_half) || (bit_input && carry_bit) || (bit_num_half && carry_bit);
        // printf("%u, %u, %u, %u, %u\n", bit_input, bit_num_half, carry, r, (1 << i));
    }

    // printf("v: %d, r: %d\n", v, r);
    return r;
}

CMP_Array bit_reverse_array(CMP_Array a)
{
    unsigned int num = a.N;
    unsigned int r = 0;
    double complex *buffer = malloc(sizeof(*buffer) * num);
    CMP_Array A = { .N = num, .data = buffer };
    for (unsigned int i=0; i < num; i++) {
        // printf("i: %u, rev(i): %u\n", i, r);
        A.data[r] = a.data[i];
        r = reverse_bit(r, num);
    }

    return A;

}