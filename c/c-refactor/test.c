#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include "array.h"

double complex number(int n);

int main(void)
{
    int num_elem = 10;
    double complex *buffer1 = malloc(sizeof(*buffer1) * num_elem);
    CMP_Array arr1 = { .N = num_elem, .data = buffer1 };

    double complex *buffer2 = malloc(sizeof(*buffer2) * num_elem);
    CMP_Array arr2 = { .N = num_elem, .data = buffer2 };

    double complex *buffer3 = malloc(sizeof(*buffer3) * num_elem);
    CMP_Array arr3 = { .N = num_elem, .data = buffer3 };


    fill_value(arr1, 1.0 + I*2.0);
    fill_value(arr2, 1.0 + I*2.0);
    print_arr(arr1);
    print_arr(arr2);
    add(arr1, arr2);
    print_arr(arr1);

    mult(arr1, arr2);
    print_arr(arr2);

    fill_func(arr3, number);
    print_arr(arr3);

    CMP_Array arr4 = slice(arr3, 2, 5);
    print_arr(arr4);

    free(arr1.data);
    free(arr2.data);
    free(arr3.data);

    return 0;

}

double complex number(int n) {
    return (double) n - I*n;
}
