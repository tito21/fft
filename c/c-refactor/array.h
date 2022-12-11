#include <stdio.h>
#include <complex.h>

typedef struct CMP_Array {
    size_t N;
    double complex* data;
} CMP_Array;

void fill_value(CMP_Array arr, double complex value);

void fill_func(CMP_Array arr, double complex (*value)(int));

void add(CMP_Array arr1, CMP_Array arr2);

void mult(CMP_Array arr1, CMP_Array arr2);

void print_arr(CMP_Array arr);

CMP_Array slice(CMP_Array arr, int m, int n);



