
// mpi_matvec_1d.c
// Parallel matrix-vector multiplication with row-wise 1-D block partitioning
// Each rank holds a contiguous block of rows.
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

static void die(const char *msg, int code){ fprintf(stderr, "%s\n", msg); MPI_Abort(MPI_COMM_WORLD, code); exit(code); }

static void compute_row_partition(int nrows, int p, int *rows_per_rank, int *displs_rows){
    // Distribute 'nrows' as evenly as possible across 'p'
    int base = nrows / p;
    int rem  = nrows % p;
    int acc = 0;
    for(int r=0;r<p;r++){
        rows_per_rank[r] = base + (r < rem ? 1 : 0);
        displs_rows[r] = acc;
        acc += rows_per_rank[r];
    }
}

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(argc < 2){
        if(rank==0) fprintf(stderr, "Usage: mpirun -np <p> ./mpi_mv1d <n>\n");
        MPI_Finalize();
        return 1;
    }
    int n = atoi(argv[1]);
    if(n <= 0){
        if(rank==0) fprintf(stderr, "n must be positive\n");
        MPI_Finalize();
        return 1;
    }

    // Determine rows per rank
    int *rows_per_rank = (int*)malloc(sizeof(int)*size);
    int *displs_rows   = (int*)malloc(sizeof(int)*size);
    compute_row_partition(n, size, rows_per_rank, displs_rows);
    int my_rows = rows_per_rank[rank];

    // Allocate local storage
    double *A_local = (double*)malloc(sizeof(double)* (long)my_rows * n);
    double *x = (double*)malloc(sizeof(double)*n);
    double *y_local = (double*)malloc(sizeof(double)*my_rows);
    if(!A_local || !x || !y_local) die("Allocation failed", 2);

    // Root allocates and initializes global A and x
    double *A = NULL;
    double *y = NULL;
    if(rank == 0){
        A = (double*)malloc(sizeof(double)*(long)n*n);
        y = (double*)malloc(sizeof(double)*n);
        if(!A || !y) die("Root allocation failed", 3);
        // Initialize A and x (deterministic)
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                A[(long)i*n + j] = (double)(i + j + 1);
            }
            x[i] = 1.0;
        }
    }

    // Broadcast x to all ranks
    if(rank != 0){
        // ensure x allocated before bcast
        for(int i=0;i<n;i++) x[i]=0.0;
    }
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Build counts/displs for Scatterv (in elements, not bytes)
    int *sendcountsA = NULL, *displsA = NULL;
    if(rank == 0){
        sendcountsA = (int*)malloc(sizeof(int)*size);
        displsA     = (int*)malloc(sizeof(int)*size);
        for(int r=0;r<size;r++){
            sendcountsA[r] = rows_per_rank[r] * n;
            displsA[r]     = displs_rows[r] * n;
        }
    }

    // Scatter rows of A
    MPI_Scatterv(
        A, sendcountsA, displsA, MPI_DOUBLE,
        A_local, my_rows * n, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    // Compute local y = A_local * x
    for(int i=0;i<my_rows;i++){
        double sum = 0.0;
        for(int j=0;j<n;j++){
            sum += A_local[(long)i*n + j] * x[j];
        }
        y_local[i] = sum;
    }

    // Gather results
    int *recvcountsY = NULL, *displsY = NULL;
    if(rank == 0){
        recvcountsY = (int*)malloc(sizeof(int)*size);
        displsY     = (int*)malloc(sizeof(int)*size);
        for(int r=0;r<size;r++){
            recvcountsY[r] = rows_per_rank[r];
            displsY[r]     = displs_rows[r];
        }
    }

    MPI_Gatherv(
        y_local, my_rows, MPI_DOUBLE,
        y, recvcountsY, displsY, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    if(rank == 0){
        int m = n < 10 ? n : 10;
        for(int i=0;i<m;i++){
            printf("y[%d] = %.4f\n", i, y[i]);
        }
    }

    // Cleanup
    if(rank == 0){
        free(A); free(y);
        free(sendcountsA); free(displsA);
        free(recvcountsY); free(displsY);
    }
    free(rows_per_rank); free(displs_rows);
    free(A_local); free(x); free(y_local);

    MPI_Finalize();
    return 0;
}
