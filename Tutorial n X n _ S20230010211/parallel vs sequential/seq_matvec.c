
// seq_matvec.c
// Sequential dense n x n matrix times n x 1 vector
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static void die(const char *msg){ fprintf(stderr, "%s\n", msg); exit(1); }

int main(int argc, char **argv){
    if(argc < 2) die("Usage: ./seq_mv <n>");
    long n = atol(argv[1]);
    if(n <= 0) die("n must be positive");
    double *A = (double*)malloc(sizeof(double)*n*n);
    double *x = (double*)malloc(sizeof(double)*n);
    double *y = (double*)malloc(sizeof(double)*n);
    if(!A || !x || !y) die("Allocation failed");

    // Initialize with simple values (deterministic)
    for(long i=0;i<n;i++){
        x[i] = 1.0; // all ones
        for(long j=0;j<n;j++){
            A[i*n + j] = (double)(i + j + 1); // any dense numbers
        }
    }

    // y = A * x
    for(long i=0;i<n;i++){
        double sum = 0.0;
        for(long j=0;j<n;j++){
            sum += A[i*n + j] * x[j];
        }
        y[i] = sum;
    }

    // Print first few entries
    long m = n < 10 ? n : 10;
    for(long i=0;i<m;i++){
        printf("y[%ld] = %.4f\n", i, y[i]);
    }

    free(A); free(x); free(y);
    return 0;
}
