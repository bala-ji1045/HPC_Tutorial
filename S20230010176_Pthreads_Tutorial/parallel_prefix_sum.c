// parallel_prefix_sum.c
// Compile: gcc -pthread -O2 parallel_prefix_sum.c -o pscan
// Run:     ./pscan 1000000 8   (N=1e6, threads=8)

#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>

typedef struct {
    int *arr, *prefix;
    long long *block_sum;
    size_t start, end;
    int block_id;
} Task;

static void *phase1_local_prefix(void *arg){
    Task *t = (Task*)arg;
    if (t->start >= t->end) { t->block_sum[t->block_id] = 0; return NULL; }

    t->prefix[t->start] = t->arr[t->start];
    for (size_t i = t->start + 1; i < t->end; i++)
        t->prefix[i] = t->prefix[i-1] + t->arr[i];

    t->block_sum[t->block_id] = t->prefix[t->end - 1];
    return NULL;
}

static void *phase3_add_offset(void *arg){
    Task *t = (Task*)arg;
    if (t->start >= t->end) return NULL;
    long long offset = t->block_sum[t->block_id];
    for (size_t i = t->start; i < t->end; i++)
        t->prefix[i] += offset;
    return NULL;
}

int main(int argc, char **argv){
    size_t N = (argc > 1) ? strtoull(argv[1], NULL, 10) : 32;
    int T = (argc > 2) ? atoi(argv[2]) : 4;
    if (T < 1) T = 1;

    int *arr = malloc(N*sizeof(int));
    int *prefix = malloc(N*sizeof(int));
    for (size_t i=0;i<N;i++) arr[i] = 1; // fill with 1's

    pthread_t *th = malloc(T*sizeof(pthread_t));
    Task *tasks = malloc(T*sizeof(Task));
    long long *block_sum = calloc(T, sizeof(long long));
    long long *block_offset = calloc(T, sizeof(long long));

    // Split into T blocks
    size_t chunk = (N + T - 1)/T;
    int blocks = T;

    // Phase 1: local prefix
    for (int b=0;b<blocks;b++){
        size_t s = b * chunk;
        size_t e = (s + chunk > N) ? N : s + chunk;
        tasks[b] = (Task){ .arr=arr, .prefix=prefix, .block_sum=block_sum,
                           .start=s, .end=e, .block_id=b };
        pthread_create(&th[b], NULL, phase1_local_prefix, &tasks[b]);
    }
    for (int b=0;b<blocks;b++) pthread_join(th[b], NULL);

    // Phase 2: compute block offsets
    block_offset[0] = 0;
    for (int b=1;b<blocks;b++)
        block_offset[b] = block_offset[b-1] + block_sum[b-1];

    for (int b=0;b<blocks;b++) block_sum[b] = block_offset[b];

    // Phase 3: add offsets
    for (int b=1;b<blocks;b++)
        pthread_create(&th[b], NULL, phase3_add_offset, &tasks[b]);
    for (int b=1;b<blocks;b++) pthread_join(th[b], NULL);

    printf("N=%zu, T=%d, final prefix sum=%d\n",
           N, T, prefix[N-1]);

    free(arr); free(prefix); free(th); free(tasks);
    free(block_sum); free(block_offset);
    return 0;
}
