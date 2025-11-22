# Sparse Matrix-Vector Multipication with OpenMP

This repository contains the code and experimental setup for a project on sparse matrix–vector multiplication (SpMV) in C with OpenMP.

## 1. Overview

The code:

* reads a sparse matrix **A** from a Matrix Market (`.mtx`) file in COO format,
* converts it to **CSR** (Compressed Sparse Row),
* constructs a dense input vector **x** with deterministic pseudo-random values in ([-1, 1]),
* computes **y = A x** using a CSR SpMV kernel parallelised with **OpenMP** over rows,
* measures performance in terms of:

  * **p90 time** Each reported time is the p90 over 15 repetitions, where every repetition consists of 1000 SpMV calls to reduce timing noise.

The same binary is used for the “sequential” baseline by running with `OMP_NUM_THREADS=1`.


## 2. Prerequisites
You need a reasonably recent GNU/Linux environment with:

- **gcc** with OpenMP support  
  (e.g. `gcc` ≥ 9.0; OpenMP 4.5+ is sufficient),
- a POSIX shell (`bash` or similar),
- standard build tools (`make`).

## 3. Building the Project

The project is built using the provided `Makefile`.

```bash
git clone https://github.com/lillo24/parallel_SpMV
cd parallel_SpMV

make
```

The Makefile uses optimisation and OpenMP flags such as:

```bash
-std=c11 -O3 -march=native -fopenmp -Wall -Wextra
```

After compilation you obtain a single executable:

* `mmread`

You can also build a debug version with symbols by running:

```bash
make DEBUG=1
```

## 4. Repository Contents

Source and build files:

* `main.c` – driver program: reads a `.mtx` file, builds CSR, runs SpMV and reports timing statistics.
* `mmio.c`, `mmio.h` – Matrix Market I/O (adapted from the NIST example).
* `Makefile` – build rules and a convenience `run` target.
* `README.md` – this document.

Matrices used in the experiments (Matrix Market format, from the SuiteSparse collection):

* `bcspwr01.mtx`
* `arc130.mtx`
* `bcspwr04.mtx`
* `ash608.mtx`
* `1138_bus.mtx`

Place these files in the repository root (they are the filenames expected by the examples below).

## 5. Running the Program

The Makefile defines a convenience rule that compiles (if needed) and runs the executable.

```bash
# default matrix: 1138_bus.mtx
make run

# run on another matrix, e.g. bcspwr01.mtx
make run FILE=bcspwr01.mtx
```

This is equivalent to running the executable directly:

```bash
./mmread bcspwr01.mtx
```

The program prints:

* basic CSR information,
* the first few elements of the result vector,
* a checksum,
* timing statistics, e.g.:

```text
SpMV: runs=15  p90=0.006 ms  mean=0.005 ms  (min=0.004 ms, max=0.010 ms)
```


## 6. OpenMP Settings

Before running, OpenMP parameters are set via environment variables.

The experiments in the report use:

```bash
export OMP_SCHEDULE=static
```

and varying `OMP_NUM_THREADS` for the sequential baseline and parallel runs, e.g.:

```bash
# sequential baseline (1 thread)
OMP_NUM_THREADS=1 ./mmread 1138_bus.mtx

# parallel runs (2, 4, 8 threads) with the same binary
OMP_NUM_THREADS=2 ./mmread 1138_bus.mtx
OMP_NUM_THREADS=4 ./mmread 1138_bus.mtx
OMP_NUM_THREADS=8 ./mmread 1138_bus.mtx
```

To compare scheduling strategies on `1138_bus.mtx` with 4 threads:

```bash
OMP_NUM_THREADS=4 OMP_SCHEDULE=static      ./mmread 1138_bus.mtx
OMP_NUM_THREADS=4 OMP_SCHEDULE=dynamic,64  ./mmread 1138_bus.mtx
OMP_NUM_THREADS=4 OMP_SCHEDULE=guided      ./mmread 1138_bus.mtx
```



## 7. Reproducing the Experiments from the Report

To reproduce the main scaling experiments with the five matrices and `schedule(static)`, you can use:

```bash
export OMP_SCHEDULE=static

for M in 1138_bus ash608 bcspwr04 arc130 bcspwr01; do
  echo "=== Matrix: $M.mtx ==="
  for T in 1 2 4 8; do
    echo "--- OMP_NUM_THREADS=$T ---"
    OMP_NUM_THREADS=$T ./mmread "$M.mtx"
  done
done
```

Each run internally performs 15 SpMV calls and reports the 90th-percentile (p90) time, which is the metric used in the report.
