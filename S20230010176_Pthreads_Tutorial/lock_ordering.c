// lock_ordering.c
// Compile: gcc -pthread -O2 lock_ordering.c -o locks
#include <stdio.h>
#include <pthread.h>

pthread_mutex_t A=PTHREAD_MUTEX_INITIALIZER, B=PTHREAD_MUTEX_INITIALIZER;

void* t1(void* a){
    pthread_mutex_lock(&A);
    pthread_mutex_lock(&B);
    printf("Thread 1 executing\n");
    pthread_mutex_unlock(&B);
    pthread_mutex_unlock(&A);
    return NULL;
}
void* t2(void* a){
    pthread_mutex_lock(&A);   // same order as t1
    pthread_mutex_lock(&B);
    printf("Thread 2 executing\n");
    pthread_mutex_unlock(&B);
    pthread_mutex_unlock(&A);
    return NULL;
}

int main(){
    pthread_t x,y;
    pthread_create(&x,NULL,t1,NULL);
    pthread_create(&y,NULL,t2,NULL);
    pthread_join(x,NULL);
    pthread_join(y,NULL);
    return 0;
}
