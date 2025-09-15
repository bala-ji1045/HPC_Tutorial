// pthread_mergesort.c
// Compile: gcc -pthread -O2 pthread_mergesort.c -o pmsort
// Run:     ./pmsort 1000000 8

#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

typedef struct {
    int *A, *B;
    int left, right;
    int depth, max_depth;
} Job;

static inline void insertion_sort(int *A, int l, int r){
    for (int i=l+1;i<r;i++){
        int x=A[i], j=i-1;
        while (j>=l && A[j]>x){ A[j+1]=A[j]; j--; }
        A[j+1]=x;
    }
}

static void merge(int *A, int *B, int l, int m, int r){
    int i=l, j=m, k=l;
    while (i<m && j<r) B[k++] = (A[i]<=A[j])? A[i++]:A[j++];
    while (i<m) B[k++] = A[i++];
    while (j<r) B[k++] = A[j++];
    for (i=l;i<r;i++) A[i]=B[i];
}

static void *merge_sort_job(void *arg);

static void spawn_or_run(Job *L, Job *R){
    if (L->depth < L->max_depth) {
        pthread_t t;
        pthread_create(&t, NULL, merge_sort_job, R);
        merge_sort_job(L);
        pthread_join(t, NULL);
    } else {
        merge_sort_job(L);
        merge_sort_job(R);
    }
}

static void *merge_sort_job(void *arg){
    Job *J = (Job*)arg;
    int n = J->right - J->left;
    if (n <= 32) { insertion_sort(J->A, J->left, J->right); return NULL; }
    int mid = J->left + n/2;

    Job JL = *J; JL.right = mid; JL.depth = J->depth + 1;
    Job JR = *J; JR.left  = mid; JR.depth = J->depth + 1;

    spawn_or_run(&JL, &JR);
    merge(J->A, J->B, J->left, mid, J->right);
    return NULL;
}

int main(int argc, char **argv){
    int N = (argc > 1) ? atoi(argv[1]) : 100000;
    int threads = (argc > 2) ? atoi(argv[2]) : 8;

    int depth = 0, t=threads; while (t>1) { depth++; t>>=1; }

    int *A = malloc(N*sizeof(int));
    int *B = malloc(N*sizeof(int));
    srand(42);
    for (int i=0;i<N;i++) A[i] = rand()%1000000;

    Job root = { .A=A, .B=B, .left=0, .right=N, .depth=0, .max_depth=depth };
    merge_sort_job(&root);

    // verify
    int ok=1; for (int i=1;i<N;i++) if (A[i-1]>A[i]){ ok=0; break; }
    printf("Sorted? %s\n", ok?"YES":"NO");

    free(A); free(B);
    return 0;
}
