#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>

/* No OpenMP for this HW */
#ifdef _OPENMP
#include <omp.h>
#endif

/* set this define to 0 to build serial code */
#define USE_MPI 1

#if USE_MPI
#include <mpi.h>
#endif

#define SEND_B 0   // tag for mpi communication of matrix B

static double timer() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);

    /* The code below is for another high resolution timer */
    /* I'm using gettimeofday because it's more portable */

    /*
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return ((double) (tp.tv_sec) + 1e-9 * tp.tv_nsec);
    */
}

int main(int argc, char **argv) {
    int rank, num_tasks;

    /* Initialize MPI */
#if USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // printf("Hello world from rank %3d of %3d\n", rank, num_tasks);
    
    // set id of next processor after the current (proc 0 and last one loop around)
    int nextProc = (rank+1) % num_tasks;
    int prevProc = (rank+num_tasks-1) % num_tasks; //can't remember how % is defined on negatives, so add num_tasks for safety
#else
    rank = 0;
    num_tasks = 1;
#endif

    if (argc != 2) {
        if (rank == 0) {
            fprintf(stderr, "%s <n>\n", argv[0]);
            fprintf(stderr, "Program for parallel dense matrix-matrix multiplication\n");
            fprintf(stderr, "with 1D row partitioning\n");
            fprintf(stderr, "<n>: matrix dimension (an nxn dense matrix is created)\n");
#if USE_MPI
            MPI_Abort(MPI_COMM_WORLD, 1);
#else
            exit(1);
#endif
        }
    }

    int n;

    n = atoi(argv[1]);
    assert(n > 0);
    assert(n < 10000);

    /* ensure that n is a multiple of num_tasks */
    n = (n/num_tasks) * num_tasks;
    
    int n_p = (n/num_tasks);

    /* print new n to let user know n has been modified */
    if (rank == 0) {
        fprintf(stderr, "n: %d, n_p: %d\n", n, n_p);
        fprintf(stderr, "Requires %3.6lf MB of memory per task\n", ((3*4.0*n_p)*n/1e6));
    }

    float *A, *B, *C;
    
    A = (float *) malloc(n_p * n * sizeof(float));
    assert(A != 0);

    B = (float *) malloc(n_p * n * sizeof(float));
    assert(B != 0);
    
    C = (float *) malloc(n_p * n * sizeof(float));
    assert(C != 0);

    /* linearized matrices in row-major storage */
    /* A[i][j] would be A[i*n+j] */

    int i, j;

    /* static initalization, so that we can verify output */
    /* using very simple initialization right now */
    /* this isn't a good check for parallel debugging */
#ifdef _OPENMP
#pragma omp parallel for private(i,j)
#endif
    for (i=0; i<n_p; i++) {
        for (j=0; j<n; j++) {
            A[i*n+j] = (rank+1);
            B[i*n+j] = 1;
            C[i*n+j] = 0;
        }
    }

#if USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    double elt = 0.0;
    if (rank == 0) 
        elt = timer();

#if USE_MPI
    /* Parallel matmul code goes here, see lecture slides for idea */
    /* The matrix C should be updated correctly */
    
    float*     tempB;   // the B received from the previous process; used so B
                        //   isn't overwritten in the receive communication
    MPI_Status status;  // status after an mpi communication
    
    // allocate tempB
    tempB = (float *) malloc(n_p * n * sizeof(float));
    assert(tempB != 0);
    
    for (int i = 0; i < num_tasks; ++i) {
        // C = C + AB
        int aOffset = (i + rank) % num_tasks;
        for(int j = 0; j < n_p; ++j) {
        // for each row of C this proc deals with
            for(int k = 0; k < n; ++k) {
            // for each element in the row
                int cIndex = j * n + k;
                for(int m = 0; m < n_p; ++m) {
                // for each row of B
                    int aIndex = j * n + aOffset + m;
                    int bIndex = m * n + k;
                    C[cIndex] += A[aIndex] * B[bIndex];
                }
            }
        }
        
        // send B to the next process, and will be receiving from the previous
        MPI_Sendrecv(B,     n_p * n, MPI_FLOAT, prevProc, SEND_B,
                     tempB, n_p * n, MPI_FLOAT, nextProc, SEND_B,
                                          MPI_COMM_WORLD, &status);
                                          
        MPI_Barrier(MPI_COMM_WORLD);
        
        // copy tempB back to B
        for(int j = 0; j < n_p * n; j++) {
            B[j] = tempB[j];
        }
    }
    
    free(tempB);

#else
    int k;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
    for (i=0; i<n_p; i++) {
        for (j=0; j<n; j++) {
            float c_ij = 0;
            for (k=0; k<n; k++) {
                c_ij += A[i*n+k]*B[k*n+j];
            }
            C[i*n+j] = c_ij;
        }
    }
#endif

    if (rank == 0) 
        elt = timer() - elt;

    /* Verify */
    int verify_failed = 0;
    for (i=0; i<n_p; i++) {
        for (j=0; j<n; j++) {
            if (C[i*n+j] != ((rank+1)*n))
                verify_failed = 1;
        }
    }

    if (verify_failed) {
        fprintf(stderr, "ERROR: rank %d, verification failed, exiting!\n", rank);
#if USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 2);
#else
        exit(2);
#endif
    }

    if (rank == 0) {
        fprintf(stderr, "Time taken: %3.3lf s.\n", elt);
        fprintf(stderr, "Performance: %3.3lf GFlop/s\n", (2.0*n*n)*n/(elt*1e9));
    }

    /* free memory */
    free(A); free(B); free(C);

    /* Shut down MPI */
#if USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
