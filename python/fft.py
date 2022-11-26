import sys
import time

import numpy as np

BCK, RUN = 0, 1

def random_array(num):
    return np.random.rand(num)*2.0-1.0 + 1j*np.random.rand(num)*2.0-1.0

def benchmark(runs, num, func):

    x = random_array(num)
    times = np.empty(runs)

    for i in range(runs):
        start = time.perf_counter()
        func(x)
        stop = time.perf_counter()
        times[i] = stop - start

    print(f"Run took: {times.mean()} ± {times.std()} s")

def fft_np(x):
    return np.fft.fft(x)

def dft_vec(x):
    N = x.shape[0]
    n = np.arange(N)
    k = n.reshape((N, 1))
    W = np.exp(-2j * np.pi * k * n / N)
    return W @ x

def fft(x):
    N = x.shape[0]
    if not ((N & (N - 1)) == 0):
        raise ValueError("Number of elements must be a power of 2")
    elif N <= 16:
        return dft_vec(x)
    else:
        X_even = fft(x[::2])
        X_odd = fft(x[1::2])
        n = np.arange(N)
        factor = np.exp(-2j*np.pi*n/N)
        return np.concatenate([X_even + factor[:N // 2] * X_odd,
                               X_even + factor[N // 2:] * X_odd])

def dft(x):
    N = x.shape[0]
    X = np.zeros(N, dtype=complex)

    for j in range(N):
        ss = 0.0 + 1j*0.0
        for k in range(N):
            w = np.exp(-j*k*2j*np.pi/N)
            ss += x[k]*w
        X[j] = ss
    return X

if __name__ == "__main__":

    mode = RUN
    num, reps = 1024, 100
    func = fft

    if len(sys.argv) > 1:

        arg = sys.argv[1]
        if arg == "FFT":
            func = fft
        elif arg == "DFT-VEC":
            func = dft_vec
        elif arg == "DFT":
            func = dft
        elif arg == "FFT-NP":
            func = fft_np
        num = int(sys.argv[2])

    if len(sys.argv) > 3:
        arg = sys.argv[3]
        if arg =="benchmark": mode = BCK
        reps = int(sys.argv[4])

    if mode == BCK:
        print("Benchmarking")
        benchmark(reps, num, func)

    elif mode == RUN:
        x = random_array(num)
        func(x)