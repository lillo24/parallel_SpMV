#ifdef _OPENMP
  #include <omp.h>
#endif


/* 
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing 
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include <string.h>
#include <time.h>

#include "mmio.h"



//HELPER FUNCTIONS
static void coo_to_csr(int M, int nnz,
                       const int *I, const int *J, const double *V,
                       int **rowptr_out, int **col_out, double **val_out)
{
    int *rowptr = (int*)calloc((size_t)M + 1, sizeof *rowptr);
    if (!rowptr) { fprintf(stderr, "CSR rowptr\n"); exit(1); }

    // Adding +1 for the row that has a value (scanning through I[])
    for (int k = 0; k < nnz; ++k) rowptr[I[k] + 1]++;

    // Enabling the rowptr[i] < rowptr[i+1] loops
    for (int i = 0; i < M; ++i) rowptr[i + 1] += rowptr[i];


    // Re-ordering values in Col_vec and Value_vec
    int *col  = (int*)   malloc((size_t)nnz * sizeof *col);
    double *av = (double*)malloc((size_t)nnz * sizeof *av);
    if (!col || !av) { fprintf(stderr, "CSR col/val\n"); exit(1); }

    int *next = (int*)malloc((size_t)M * sizeof *next);
    if (!next) { fprintf(stderr, "CSR next\n"); exit(1); }
    memcpy(next, rowptr, (size_t)M * sizeof *next);

    for (int k = 0; k < nnz; ++k) {
        int r = I[k];
        int p = next[r]++;
        col[p] = J[k];
        av[p]  = V[k];
    }
    free(next);

    *rowptr_out = rowptr;
    *col_out    = col;
    *val_out    = av;
}








/***** Parallel SpMV ******/

void spmv_csr(int M,
              const int *rowptr, const int *col, const double *value_nnzs,
              const double *x, double *y)
{
    #pragma omp parallel for schedule(runtime)
    /* sum is private per iteration; rowptr, col, value_nnzs, x are read-only, so no data races */
    for (int i = 0; i < M; ++i) {
        double sum = 0;
        for (int p = rowptr[i]; p < rowptr[i+1]; ++p)
            sum += value_nnzs[p] * x[col[p]];
        y[i] = sum;
    }
}





/***** Analysis ******/ 

// Get Seconds
static inline double now_sec(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9*ts.tv_nsec;
#endif
}

static int cmp_dbl(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return  1;
    return 0;
}


//Sorts the array of run-times and return p-th value.
static double percentile_ms(double *ms, int n, double p) {
    qsort(ms, (size_t)n, sizeof(double), cmp_dbl);
    int idx = (int)((p*(n-1)) + 0.5);   // nearest rank
    return ms[idx];
}

//Mean ms
static double mean_ms(const double *ms, int n) {
    double s=0; for (int i=0;i<n;++i) s+=ms[i]; return s/n;
}




int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    bool is_symmetric = false, is_pattern = false;

    FILE *f;
    int M, N, nz;   
    int *I, *J;
    double *val;
    
    if (argc < 2) {
        fprintf(stderr, "Usage: %s [matrix-market-filename]\n", argv[0]);
        return EXIT_FAILURE;
    }
    f = fopen(argv[1], "r");
    if (!f) {
        perror("fopen");
        return EXIT_FAILURE;
    }

    if (mm_read_banner(f, &matcode) != 0) {
        fprintf(stderr, "Could not process Matrix Market banner.\n");
        return EXIT_FAILURE;
    }


    //My version
    if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode)){
        fprintf(stderr, "Only  sparse matrix (coordinate) supported.\n"); exit(1);
    }

    if (mm_is_complex(matcode)) {
        fprintf(stderr, "Complex matrices not supported for this exercise.\n"); exit(1);
    }


    // Check if pattern
    if (mm_is_pattern(matcode)) {
        is_pattern = true;
    }

    // Check if Symmetry
    if (mm_is_symmetric(matcode)) {
        is_symmetric = true;
    }




    //find out size of sparse matrix...

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    //reseve memory (64bytes alignment for cache efficiency) for matrices

    int cap = is_symmetric ? 2 * nz : nz;
    size_t bytesI = (size_t)cap * sizeof *I;
    size_t bytesJ = (size_t)cap * sizeof *J;
    size_t bytesV = (size_t)cap * sizeof *val;
    if (posix_memalign((void**)&I,   64, bytesI) ||
        posix_memalign((void**)&J,   64, bytesJ) ||
        posix_memalign((void**)&val, 64, bytesV)) {
        fprintf(stderr,"posix_memalign failed\n"); exit(1);
    }




    /************************/
    /* my algorithm */
    /************************/

    int nnz_eff = 0; //If it is diagonal we don't know how many elements are on the diagonal, so we need to "check it"
    for (int t = 0; t < nz; ++t) {
        int ii, jj; double dv = 1.0;

        //Check for pattern
        if (is_pattern) {
            if (fscanf(f, "%d %d", &ii, &jj) != 2) { fprintf(stderr,"read error\n"); exit(1); }
        // Check for integer
        } else if (mm_is_integer(matcode)) {
            long long tmp;
            if (fscanf(f, "%d %d %lld", &ii, &jj, &tmp) != 3) { fprintf(stderr,"read error\n"); exit(1); }
            dv = (double)tmp;
        } else {
            if (fscanf(f, "%d %d %lf", &ii, &jj, &dv) != 3) { fprintf(stderr,"read error\n"); exit(1); }
        }

        --ii; --jj; // 1-based -> 0-based
        I[nnz_eff] = ii; 
        J[nnz_eff] = jj;
        val[nnz_eff] = dv; ++nnz_eff;
        
        // Expand for symmetry
        if (is_symmetric && ii != jj) {
            I[nnz_eff] = jj; 
            J[nnz_eff] = ii; 
            val[nnz_eff] = dv;
            ++nnz_eff;
        }
    }

    if (f != stdin) fclose(f);






    /************************/
    /* Coo --> CSR */
    /************************/

    int *rowptr = NULL, *col = NULL;
    double *value_nnzs = NULL;
    coo_to_csr(M, nnz_eff, I, J, val, &rowptr, &col, &value_nnzs);





    /************************/
    /* Make randomly generated vector (+ final result vector) 
    /************************/

    double *random_vector = NULL, *result = NULL;

    if (posix_memalign((void**)&random_vector, 64, (size_t)N * sizeof *random_vector) ||
        posix_memalign((void**)&result, 64, (size_t)M * sizeof *result)) {
        fprintf(stderr, "posix_memalign x/y failed\n"); 
        exit(1);
    }

    srand(12345u);
    for (int j = 0; j < N; ++j)
        random_vector[j] = 2.0 * (rand() / (double)RAND_MAX) - 1.0;




    /************************/
    /* TIMING */
    /************************/

    const int runs = 15;
    const int inner_reps = 1000;
    double time_runs[runs];

    // Warm-up

    for (int w = 0; w < 3; ++w) {
        spmv_csr(M, rowptr, col, value_nnzs, random_vector, result);
    }
    
    //Measured Runs
    for (int r = 0; r < runs; ++r) {
        double t0 = now_sec();

        for (int it = 0; it < inner_reps; ++it) { // To avoid noise
            spmv_csr(M, rowptr, col, value_nnzs, random_vector, result);
        }
        
        double t1 = now_sec();
        time_runs[r] = (t1 - t0) * 1e3 / (double)inner_reps; // milliseconds
    }





    /************************/
    /* Logging */
    /************************/

#if 0

    /* Matrix metadata */
    double density = (double)nnz_eff / ((double)M * (double)N);
    fprintf(stderr,
        "META: M=%d N=%d nnz=%d density=%g\n",
        M, N, nnz_eff, density);

    // Check if pattern
    if (is_pattern) {
        printf("Pattern Matrix.\n");
    }

    // Check if Symmetry
    if (is_symmetric) {
        printf("Symmetric Matrix.\n");
    }


    fprintf(stderr, "CSR: M=%d, nnz=%d, rowptr[M]=%d\n", M, nnz_eff, rowptr[M]);
    for (int i = 0; i < (M < 5 ? M : 5); ++i)
        fprintf(stderr, "row %d: [%d..%d)\n", i, rowptr[i], rowptr[i+1]);

    for (int i = 0; i < (M < 5 ? M : 5); ++i)
        fprintf(stderr, "result[%d] = %.6f\n", i, result[i]);
#endif


    //Check sum
    double chk = 0.0;
    #pragma omp parallel for reduction(+:chk)
    for (int i = 0; i < M; ++i) chk += result[i];
    fprintf(stderr, "checksum=%.17g\n", chk);


    //Stats in ms
    double p90 = percentile_ms(time_runs, runs, 0.90);
    double mean = mean_ms(time_runs, runs);

    fprintf(stderr,
            "SpMV: runs=%d  inner_reps=%d  p90=%.3f ms  mean=%.3f ms  (min=%.3f ms, max=%.3f ms)\n",
            runs, inner_reps, p90, mean, time_runs[0], time_runs[runs-1]);



    // Performance model: GFLOP/s and bandwitdh
    double t_p90_s = p90 / 1e3;
    double flops  = 2.0 * (double)nnz_eff;
    double gflops = flops / (t_p90_s * 1e9); //GFLOP/s


    /*
    8B*nnz for values
    4B*nnz for column index (col[p]) 
    8B*nnz for x[]

    8B*M for y[]
    */  
    double bytes  = 20.0 * (double)nnz_eff + 8.0 * (double)M;
    double bw_GBs = bytes / (t_p90_s * 1e9);

    fprintf(stderr,
            "Perf (p90): ~%.3f GFLOP/s, ~%.3f GB/s (approx)\n",
            gflops, bw_GBs);



    /************************/
    /* Simple summary line for results.txt */
    /************************/

    // Read thread count and schedule from environment
    const char *threads_env = getenv("OMP_NUM_THREADS");
    int threads = threads_env ? atoi(threads_env) : 1;

    const char *sched_env = getenv("OMP_SCHEDULE");
    if (!sched_env) sched_env = "static";

    // Matrix filename (argv[1])
    const char *matrix_name = (argc >= 2) ? argv[1] : "unknown";

    // One-line summary (goes to stdout)
    printf("BENCH matrix=%s threads=%d schedule=%s runs=%d inner_reps=%d p90_ms=%.6f gflops=%.6f bw_GBs=%.6f\n",
       matrix_name, threads, sched_env, runs, inner_reps, p90, gflops, bw_GBs);


           
    //Freeing
    free(rowptr);
    free(col);
    free(value_nnzs);
    free(I); free(J); free(val);
    free(random_vector); 
    free(result);
	return 0;

}
