#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void print_complex_array(double complex *arr, size_t N) {

    printf("[ ");
    for (int i=0; i < N; i++) {
        printf("(%f,%f), ", creal(arr[i]), cimag(arr[i]));
    }
    printf("]\n");

}

double complex *concatenate(double complex *arr1, size_t N1, double complex *arr2, size_t N2) {

    double complex *out = (double complex *) malloc(sizeof(double complex)*(N1 + N2));

    for (int i=0; i < N1; i++) {
        out[i] = arr1[i];
    }

    for (int i=0; i < N2; i++) {
        out[i+N1] = arr2[i];
    }
    return out;
}

double complex *slice(double complex *arr, int m, int n)
{
    int range = n+1 - m;
    double complex *out = (double complex *) malloc(sizeof(double complex)*range);

    for (int i=0; i < range; i++) {
        out[i] = arr[m+i];
    }
    return out;

}

double complex *add(double complex *arr1, double complex *arr2, size_t N) {

    double complex *out = (double complex *) malloc(sizeof(double complex)*N);
    for (int i=0; i < N; i++) {
        out[i] = arr1[i] + arr2[i];
    }
    return out;
}

double complex *mult(double complex *arr1, double complex *arr2, size_t N) {

    double complex *out = (double complex *) malloc(sizeof(double complex)*N);
    for (int i=0; i < N; i++) {
        out[i] = arr1[i] * arr2[i];
    }
    return out;
}

double complex *get_even(double complex *arr, size_t N) {

    double complex *out = (double complex *) malloc(sizeof(double complex)*N/2);
    for (int i=0; i < N; i += 2) {
        out[i/2] = arr[i];
    }
    return out;
}

double complex *get_odd(double complex *arr, size_t N) {

    double complex *out = (double complex *) malloc(sizeof(double complex)*N/2);
    for (int i=1; i < N; i += 2) {
        out[i/2] = arr[i];
    }
    return out;
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

double complex *random_array(size_t N)
{
    double complex *arr = (double complex *) malloc(sizeof(double complex)*N);

    for (int i=0; i < N; i++) {
        double num1 = (double) rand() / RAND_MAX * 2.0 - 1.0;
        double num2 = (double) rand() / RAND_MAX * 2.0 - 1.0;
        arr[i] = num1 + I*num2;
    }

    return arr;
}