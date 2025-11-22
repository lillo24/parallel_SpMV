# Sparse Matrix-Vector Multipication with OpenMP

This repository contains the code and experimental setup for the project

## 1. Overview

The code:

- reads a sparse matrix **A** from a Matrix Market (`.mtx`) file in COO format,
- converts it to **CSR** (Compressed Sparse Row),
- constructs a dense input vector **x** with deterministic pseudo-random values in \([-1, 1]\),
- computes **y = A x** using:
  - a **sequential** CSR SpMV kernel, and
  - a **parallel** CSR SpMV kernel with **OpenMP**, parallelising over rows,
- measures performance in terms of:
  - **p90 time** (90th-percentile time per SpMV, in ms),
  - **GFLOP/s** (throughput, assuming 2 FLOPs per non-zero),
  - **effective memory bandwidth (GB/s)** from a simple byte-count model.

## 2. Prerequisites

You need a reasonably recent GNU/Linux environment with:

- **gcc** with OpenMP support  
  (e.g. `gcc` ≥ 9.0; OpenMP 4.5+ is sufficient),
- a POSIX shell (`bash` or similar),
- standard build tools (`make`).


## 3. Building the Project

The project is built using a standard `Makefile`.

```bash
git clone https://github.com/lillo24/delivery1_leonardocolli
cd delivery1_leonardocolli

make
```

The Makefile uses the same compiler flags described in the report, e.g.:
```bash
-std=c11 -O3 -march=native -fopenmp -Wall -Wextra
```

After compilation, you should obtain a single executable
