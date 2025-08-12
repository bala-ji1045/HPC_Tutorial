#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 1000000
#define INSERTION_SORT_THRESHOLD 32

void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void insertionSort(int arr[], int low, int high) {
    for (int i = low + 1; i <= high; i++) {
        int key = arr[i];
        int j = i - 1;
        while (j >= low && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
}

int partition(int arr[], int low, int high) {
    int pivot = arr[high];
    int i = (low - 1);
    for (int j = low; j < high; j++) {
        if (arr[j] < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

void quickSort(int arr[], int low, int high) {
    while (low < high) {
        if (high - low + 1 < INSERTION_SORT_THRESHOLD) {
            insertionSort(arr, low, high);
            break;
        } else {
            int pi = partition(arr, low, high);
            if (pi - low < high - pi) {
                quickSort(arr, low, pi - 1);
                low = pi + 1;
            } else {
                quickSort(arr, pi + 1, high);
                high = pi - 1;
            }
        }
    }
}

int main() {
    int *arr = (int*)malloc(SIZE * sizeof(int));
    if (!arr) {
        printf("Memory allocation failed\n");
        return 1;
    }
    srand(time(NULL));
    for (int i = 0; i < SIZE; i++) {
        arr[i] = rand();
    }
    clock_t start, end;
    start = clock();
    quickSort(arr, 0, SIZE - 1);
    end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Hybrid QuickSort+InsertionSort of %d elements took %.2f seconds.\n", SIZE, time_taken);
    free(arr);
    return 0;
}
