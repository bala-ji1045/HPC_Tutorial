// parallel_sum.c
// Compile: gcc -pthread -O2 parallel_sum.c -o psum
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#define MAX_THREADS 16
#define N 100000000

long long *arr;
long long sum[MAX_THREADS]={0};

typedef struct {
    int start, end, tid;
} ThreadData;

void* sum_array(void* arg){
    ThreadData *d=(ThreadData*)arg;
    long long local=0;
    for (int i=d->start;i<d->end;i++) local+=arr[i];
    sum[d->tid]=local;
    return NULL;
}

double Time_T(int threads){
    pthread_t th[threads];
    ThreadData data[threads];
    clock_t st=clock();
    int chunk=N/threads;

    for(int i=0;i<threads;i++){
        data[i].start=i*chunk;
        data[i].end=(i==threads-1)?N:(i+1)*chunk;
        data[i].tid=i;
        pthread_create(&th[i],NULL,sum_array,&data[i]);
    }
    long long total=0;
    for(int i=0;i<threads;i++){ pthread_join(th[i],NULL); total+=sum[i]; }
    clock_t et=clock();
    return (double)(et-st)/CLOCKS_PER_SEC;
}

int main(){
    arr=malloc(N*sizeof(long long));
    for(int i=0;i<N;i++) arr[i]=1;
    for(int t=1;t<=MAX_THREADS;t++){
        double tm=Time_T(t);
        printf("Threads=%d, Time=%lf sec\n",t,tm);
    }
    free(arr);
    return 0;
}
