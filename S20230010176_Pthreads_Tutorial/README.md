# Pthreads Tutorial Solutions

This folder contains working C implementations for the tutorial problems.

---

## 1. Dining Philosophers (deadlock-free)

File: `dining_philosophers.c`

- Uses a waiter semaphore (`N-1` permits).
- Ensures no deadlock while philosophers eat & think.

```bash
gcc -pthread -O2 dining_philosophers.c -o dining
./dining 5

## 2. prefix sum

gcc -pthread -O2 parallel_prefix_sum.c -o pscan
./pscan 1000000 8

## 3. parallel merge sort

gcc -pthread -O2 pthread_mergesort.c -o pmsort
./pmsort 1000000 8

## 4. shared counter with mutex

gcc -pthread -O2 counter_mutex.c -o counter
./counter

## 5. dead lock fixing with lock ordering

gcc -pthread -O2 lock_ordering.c -o locks
./locks

## 6. multi threaded array sum

gcc -pthread -O2 parallel_sum.c -o psum
./psum