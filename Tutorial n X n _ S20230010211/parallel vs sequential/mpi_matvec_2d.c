
// mpi_matvec_2d.c
// Parallel matrix-vector multiplication with 2-D block partitioning (q x q process grid)
// Assumes p = q^2 and n divisible by q.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

static void die(const char *msg, int code){ fprintf(stderr, "%s\n", msg); MPI_Abort(MPI_COMM_WORLD, code); exit(code); }

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);
    int rank, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if(argc < 2){
        if(rank==0) fprintf(stderr, "Usage: mpirun -np <q*q> ./mpi_mv2d <n>\n");
        MPI_Finalize();
        return 1;
    }
    int n = atoi(argv[1]);
    if(n <= 0){
        if(rank==0) fprintf(stderr, "n must be positive\n");
        MPI_Finalize();
        return 1;
    }

    // Determine q such that p = q*q
    int q = (int)(sqrt((double)p) + 0.5);
    if(q*q != p){
        if(rank==0) fprintf(stderr, "Number of processes p=%d is not a perfect square.\n", p);
        MPI_Finalize();
        return 2;
    }
    if(n % q != 0){
        if(rank==0) fprintf(stderr, "For simplicity this demo requires n divisible by q=%d.\n", q);
        MPI_Finalize();
        return 3;
    }
    int b = n / q; // block size

    // Create a 2D Cartesian grid
    int dims[2] = {q, q};
    int periods[2] = {0, 0};
    MPI_Comm grid;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid);

    int coords[2];
    MPI_Cart_coords(grid, rank, 2, coords);
    int myrow = coords[0];
    int mycol = coords[1];

    // Row and column communicators
    MPI_Comm row_comm, col_comm;
    int row_color = myrow;
    int col_color = mycol;
    MPI_Comm_split(MPI_COMM_WORLD, row_color, mycol, &row_comm);
    MPI_Comm_split(MPI_COMM_WORLD, col_color, myrow, &col_comm);

    // Allocate local blocks
    double *A_local = (double*)malloc(sizeof(double) * (long)b * b);
    double *x_local = (double*)malloc(sizeof(double) * b);
    double *y_partial = (double*)malloc(sizeof(double) * b);
    double *y_block = (mycol==0) ? (double*)malloc(sizeof(double)*b) : NULL;
    if(!A_local || !x_local || !y_partial || (mycol==0 && !y_block)){
        die("Allocation failed", 4);
    }

    // Root holds global A and x for initialization and distribution
    double *A = NULL;
    double *x = NULL;
    double *y = NULL;
    if(rank == 0){
        A = (double*)malloc(sizeof(double)*(long)n*n);
        x = (double*)malloc(sizeof(double)*n);
        y = (double*)malloc(sizeof(double)*n);
        if(!A || !x || !y) die("Root allocation failed", 5);
        for(int i=0;i<n;i++){
            x[i] = 1.0;
            for(int j=0;j<n;j++){
                A[(long)i*n + j] = (double)(i + j + 1);
            }
        }
    }

    // Distribute A blocks using subarray datatypes with point-to-point sends
    // Each destination defines its own subarray view into A
    if(rank == 0){
        for(int r=0;r<p;r++){
            int ccoords[2];
            int rcoords[2];
            MPI_Cart_coords(grid, r, 2, ccoords); // reuse variable name
            int rr = ccoords[0];
            int cc = ccoords[1];
            int starts[2]   = { rr*b, cc*b };
            int sizes[2]    = { n, n };
            int subsizes[2] = { b, b };
            MPI_Datatype block_t;
            MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &block_t);
            MPI_Type_commit(&block_t);
            if(r == 0){
                // Copy my own block into A_local using MPI_Pack with same datatype into a contiguous buffer is complex;
                // we'll just copy via loops for clarity.
                for(int ii=0; ii<b; ++ii){
                    for(int jj=0; jj<b; ++jj){
                        A_local[ii*b + jj] = A[(long)(rr*b + ii)*n + (cc*b + jj)];
                    }
                }
            }else{
                MPI_Send(A, 1, block_t, r, 777, MPI_COMM_WORLD);
            }
            MPI_Type_free(&block_t);
        }
    }else{
        MPI_Recv(A_local, b*b, MPI_DOUBLE, 0, 777, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Distribute x blocks: root scatters to row 0, then broadcast down each column
    if(myrow == 0){
        if(rank == 0){
            // scatter contiguous chunks of size b to q processes along the top row
            // Build a communicator for the top row (myrow==0)
        }
    }
    // Build communicator for top row
    MPI_Comm top_row_comm;
    int top_row_color = (myrow == 0) ? 0 : MPI_UNDEFINED;
    MPI_Comm_split(MPI_COMM_WORLD, top_row_color, mycol, &top_row_comm);

    if(myrow == 0){
        if(rank == 0){
            // Prepare send counts and displs for scatter to the top row ranks in order of columns
            double *x_top = x; // contiguous
            // Gather ranks of top row in column order
            int *top_ranks = (int*)malloc(sizeof(int)*q);
            for(int c=0;c<q;c++){
                int crd[2] = {0, c};
                int rr;
                MPI_Cart_rank(grid, crd, &rr);
                top_ranks[c] = rr;
            }
            // If communicator is ordered by mycol, we can MPI_Scatter directly on top_row_comm
            // Create a temporary buffer for the top row ranks to receive in rank order of top_row_comm
            // But MPI_Scatter works within the communicator with ranks 0..q-1. We need to ensure rank 0 in top_row_comm is the process with (0,0).
            // We split with key=mycol, so in top_row_comm, rank equals mycol. Thus rank 0 is (0,0) indeed.
            MPI_Scatter(x_top, b, MPI_DOUBLE, x_local, b, MPI_DOUBLE, 0, top_row_comm);
            free(top_ranks);
        }else{
            // Top row but not global root
            MPI_Scatter(NULL, b, MPI_DOUBLE, x_local, b, MPI_DOUBLE, 0, top_row_comm);
        }
    }

    // Broadcast x_local down the column
    MPI_Bcast(x_local, b, MPI_DOUBLE, 0, col_comm);

    // Compute y_partial = A_local * x_local
    for(int i=0;i<b;i++){
        double sum = 0.0;
        for(int j=0;j<b;j++){
            sum += A_local[i*b + j] * x_local[j];
        }
        y_partial[i] = sum;
    }

    // Reduce partial sums across columns to column 0 (one per row)
    MPI_Reduce(y_partial, y_block, b, MPI_DOUBLE, MPI_SUM, 0, row_comm);

    // Gather y blocks at global root from column-0 processes
    if(mycol == 0){
        if(rank == 0){
            // Gather in row order
            double *y_blocks = y; // final vector
            MPI_Gather(y_block, b, MPI_DOUBLE, y_blocks, b, MPI_DOUBLE, 0, col_comm); // gather along column 0 communicator
        }else{
            MPI_Gather(y_block, b, MPI_DOUBLE, NULL, b, MPI_DOUBLE, 0, col_comm);
        }
    }

    if(rank == 0){
        int m = n < 10 ? n : 10;
        for(int i=0;i<m;i++){
            printf("y[%d] = %.4f\n", i, y[i]);
        }
    }

    if(top_row_comm != MPI_COMM_NULL) MPI_Comm_free(&top_row_comm);
    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&grid);
    if(rank == 0){ free(A); free(x); free(y); }
    free(A_local); free(x_local); free(y_partial); if(y_block) free(y_block);

    MPI_Finalize();
    return 0;
}
