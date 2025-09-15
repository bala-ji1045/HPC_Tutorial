// dining_philosophers.c
// Compile: gcc -pthread -O2 dining_philosophers.c -o dining
// Run:     ./dining 5  (default N=5 if omitted)

#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>

typedef struct {
    int N;
    int id;
} Philosopher;

static pthread_mutex_t *chop;
static pthread_mutex_t waiter_lock = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  waiter_cv   = PTHREAD_COND_INITIALIZER;
static int permits; // == N-1

// simple counting semaphore
static void waiter_acquire(void){
    pthread_mutex_lock(&waiter_lock);
    while (permits == 0) pthread_cond_wait(&waiter_cv, &waiter_lock);
    permits--;
    pthread_mutex_unlock(&waiter_lock);
}
static void waiter_release(void){
    pthread_mutex_lock(&waiter_lock);
    permits++;
    pthread_cond_signal(&waiter_cv);
    pthread_mutex_unlock(&waiter_lock);
}

void think(int id){ usleep(1000 * (rand()%30 + 5)); }
void eat(int id){   usleep(1000 * (rand()%30 + 5)); }

void *philosopher(void *arg){
    Philosopher *p = (Philosopher*)arg;
    int left  = p->id;
    int right = (p->id + 1) % p->N;

    for (int r = 0; r < 5; r++) { // each eats 5 times
        think(p->id);

        waiter_acquire(); // permission

        int first  = left < right ? left : right;
        int second = left ^ right ^ first;

        pthread_mutex_lock(&chop[first]);
        pthread_mutex_lock(&chop[second]);

        printf("Philosopher %d is eating (round %d)\n", p->id, r+1);
        eat(p->id);

        pthread_mutex_unlock(&chop[second]);
        pthread_mutex_unlock(&chop[first]);

        waiter_release();
    }
    return NULL;
}

int main(int argc, char **argv){
    int N = (argc > 1) ? atoi(argv[1]) : 5;
    if (N < 2) N = 5;

    chop = calloc(N, sizeof(pthread_mutex_t));
    for (int i = 0; i < N; i++) pthread_mutex_init(&chop[i], NULL);

    pthread_t *ths = malloc(N*sizeof(pthread_t));
    Philosopher *ps = malloc(N*sizeof(Philosopher));

    permits = N - 1;

    for (int i = 0; i < N; i++){
        ps[i] = (Philosopher){ .N = N, .id = i };
        pthread_create(&ths[i], NULL, philosopher, &ps[i]);
    }
    for (int i = 0; i < N; i++) pthread_join(ths[i], NULL);

    for (int i = 0; i < N; i++) pthread_mutex_destroy(&chop[i]);
    free(chop); free(ths); free(ps);
    return 0;
}
