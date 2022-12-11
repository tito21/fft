#include "array.h"

void fill_value(CMP_Array arr, double complex value)
{
    int N = arr.N;
    for (int i=0; i < N; i++) {
        arr.data[i] = value;
    }
}

void fill_func(CMP_Array arr, double complex (*value)(int))
{
    int N = arr.N;
    for (int i=0; i < N; i++) {
        arr.data[i] = value(i);
    }
}

void add(CMP_Array arr1, CMP_Array arr2)
{
    if (arr1.N != arr2.N) {
        printf("Array length must be the same.\n");
        return;
    }
    int N = arr1.N;
    for (int i=0; i < N; i++) {
        arr1.data[i] += arr2.data[i];
    }
}

void mult(CMP_Array arr1, CMP_Array arr2)
{
    if (arr1.N != arr2.N) {
        printf("Array length must be the same.\n");
        return;
    }
    int N = arr1.N;
    for (int i=0; i < N; i++) {
        arr1.data[i] *= arr2.data[i];
    }
}

void print_arr(CMP_Array arr)
{
    int N = arr.N;
    printf("[ ");
    for (int i=0; i < N; i++) {
        printf("(%f,%f), ", creal(arr.data[i]), cimag(arr.data[i]));
    }
    printf("]\n");

}

CMP_Array slice(CMP_Array arr, int m, int n)
{
    int range = n+1 - m;
    CMP_Array out = {.N = range, .data = &arr.data[m] };
    return out;
}
