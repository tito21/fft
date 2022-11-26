#include <chrono>
#include <complex>
#include <iostream>
#include <functional>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include "vector_func.cpp"

enum MODE { BCK, RUN };

void benchmark(int runs, int num, std::function<std::vector<std::complex<double>>(std::vector<std::complex<double>>)> func)
{
    auto x = random_array(num);
    std::vector<float> times(runs);
    using FpSeconds = std::chrono::duration<float, std::chrono::seconds::period>;

    static_assert(std::chrono::treat_as_floating_point<FpSeconds::rep>::value,
                  "Rep required to be floating point");

    for (int i=0; i < runs; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        func(x);
        auto stop = std::chrono::high_resolution_clock::now();
        times[i] = FpSeconds(stop - start).count();
    }

    std::cout << "Run took: " << mean_vector(times) << " ± " << std_vector(times) << " s" << std::endl;

}

std::vector<std::complex<double>> dft(std::vector<std::complex<double>> const &x)
{
    int N = x.size();
    std::vector<std::complex<double>> X(N);
    std::complex<double> I(0.0, 1.0);
    std::complex<double> w = std::exp(-2.0 * I * std::numbers::pi / (double)N);;

    for (int j=0; j < N; j++) {
        std::complex<double> ss(0.0, 0.0);
        for (int k=0; k < N; k++) {
            ss += x[k] * std::pow(w, k*j);
        }
        X[j] = ss;
    }
    return X;
}

std::vector<std::complex<double>> fft(std::vector<std::complex<double>> const &x)
{
    int N = x.size();
    std::complex<double> I(0.0, 1.0);
    if (!((N & (N - 1)) == 0)) {
        throw std::invalid_argument("Number of elements must be a power of 2");
    }
    else if (N < 16) {
        return dft(x);
    }
    else {
        auto X_even = fft(get_even(x));
        auto X_odd = fft(get_odd(x));

        std::vector<std::complex<double>> factor(N);
        for (int n=0; n < N; n++) {
            factor[n] = std::exp(-2.0 * std::numbers::pi * I * (double)  n / (double) N);
        }

        auto sum1 = mult_vector(slice_vector(factor, 0, (N / 2) - 1), X_odd);
        auto first = add_vector(X_even, sum1);

        auto sum2 = mult_vector(slice_vector(factor, N / 2, N - 1), X_odd);
        auto last = add_vector(X_even, sum2);

        first.insert(first.end(), last.begin(), last.end());
        return first;
    }

}

int main(int argc, char** argv)
{

    MODE mode = RUN;
    int num, reps;
    std::function<std::vector<std::complex<double>>(std::vector<std::complex<double>>)> func;

    if (argc > 1) {

        std::string arg(argv[1]);
        if (arg == "FFT") {
            func = fft;
        }
        else if (arg == "DFT") {
            func = dft;
        }
        else {
            throw std::invalid_argument("Must pass a valid algorithm name (FFT, DFT)");
        }
        num = std::stoi(argv[2]);

    }
    if (argc > 3) {
        std::string arg(argv[3]);
        if (arg == "benchmark") mode = BCK;
        reps = std::stoi(argv[4]);

    }

    switch (mode)
    {
    case BCK:
        std::cout << "Benchmarking" << std::endl;
        benchmark(reps, num, func);
        break;

    case RUN:
        auto x = random_array(num);
        func(x);
    }
    return 0;

}