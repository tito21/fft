#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>
#include <random>
#include <stdexcept>


std::vector<std::complex<double>> random_array(size_t N)
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dis(-1.0, 1.0);

    std::vector<std::complex<double>> arr(N);
    std::generate(arr.begin(), arr.end(), [&dis, &gen]() { return std::complex<double>(dis(gen), dis(gen)); });

    return arr;
}

template<typename T>
void print_vector(std::vector<T> const vec)
{
    std::cout << "[ ";
    for (auto x : vec) {
        std::cout << x << ", ";
    }
    std::cout << " ]\n";
}

template<typename T>
std::vector<T> get_odd(std::vector<T> const &arr)
{
    int N = arr.size();
    std::vector<T> arr_odd(N/2);
    for (int i=1; i < N; i += 2) {
        arr_odd[i/2] = arr[i];
    }
    return arr_odd;
}

template<typename T>
std::vector<T> get_even(std::vector<T> const &arr)
{
    int N = arr.size();
    std::vector<T> arr_even(N/2);
    for (int i=0; i < N; i += 2) {
        arr_even[i/2] = arr[i];
    }
    return arr_even;
}

template<typename T>
std::vector<T> mult_vector(std::vector<T> const &arr1, std::vector<T> const &arr2)
{
    if (arr1.size() != arr2.size())
    {
        throw std::invalid_argument("Vector sizes must be equal");
    }

    std::vector<T> out(arr1.size());
    std::transform(arr1.begin(), arr1.end(), arr2.begin(), out.begin(), std::multiplies<T>());

    return out;
}

template <typename T>
std::vector<T> add_vector(std::vector<T> const &arr1, std::vector<T> const &arr2)
{
    if (arr1.size() != arr2.size())
        {
            throw std::invalid_argument("Vector sizes must be equal");
        }

    std::vector<T> out(arr1.size());
    std::transform(arr1.begin(), arr1.end(), arr2.begin(), out.begin(), std::plus<T>());

    return out;

}

template <typename T>
T mean_vector(std::vector<T> const &arr)
{
    T sum = 0.0;
    for (auto v : arr) sum += v;
    return sum / arr.size();
}

template <typename T>
T std_vector(std::vector<T> const &arr)
{
    auto mean = mean_vector(arr);
    T sum = 0.0;
    for (auto v : arr) sum += (v - mean)*(v - mean);
    return std::sqrt(sum / arr.size());
}

template <typename T>
std::vector<T> slice_vector(std::vector<T> const &arr, int m, int n)
{
    auto first = arr.begin() + m;
    auto last = arr.begin() + n + 1;

    std::vector<T> out(first, last);
    return out;
}