#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <chrono>
#include <future>
#include <omp.h>

#define LIMIT 1000 // max value of array elements
#define DOT_LENGTH 35 // Length for displaying dots
#define CORES_NUMBER 10

#define MERGE_SORT 0
#define PARALLEL_MERGE_SORT 1
#define PARALLEL_MERGE_SORT_OMP 2
#define QUICK_SORT_SINGLE_THREAD 3
#define QUICK_SORT_PARALLEL_WITH_FUTURE 4
#define PARALLEL_MERGE_SORT_THREAD 5
const int FUNCTIONS_NUMBER = 6;

long long arrayMinSizeForThread = 100000;

void merge(int* arr, long long first, long long last) {
    long long middle = first + (last - first) / 2;
    long long size = last - first + 1;
    int* temp = new int[size];
    long long left = first;
    long long right = middle + 1;

    for (long long i = 0; i < size; i++) {
        if (left > middle)
            temp[i] = arr[right++];
        else if (right > last)
            temp[i] = arr[left++];
        else
            temp[i] = arr[left] < arr[right] ? arr[left++] : arr[right++];
    }

    for (long long i = 0, j = first; i < size; i++, j++)
        arr[j] = temp[i];

    delete[] temp;
}

void mergeSort(int* arr, long long first, long long last) {
    if (first < last) {
        int middle = first + (last - first) / 2;
        mergeSort(arr, first, middle);
        mergeSort(arr, middle + 1, last);
        merge(arr, first, last);
    }
}

void parallelMergeSort(int* arr, long long first, long long last) {
    if (first < last) {
        long long middle = first + (last - first) / 2;

        if (last - first >= arrayMinSizeForThread) { // Use min fixed size of array to handle with new thread
            std::future<void> f1 = std::async(std::launch::async, [&] {
                parallelMergeSort(arr, first, middle);
                });
            parallelMergeSort(arr, middle + 1, last);
            //f1.wait();
        }
        else {
            mergeSort(arr, first, middle);
            mergeSort(arr, middle + 1, last);
        }

        merge(arr, first, last);
    }
}

void parallelMergeSortWithThread(int* arr, long long first, long long last) {
    if (first < last) {
        long long middle = first + (last - first) / 2;

        if (last - first >= arrayMinSizeForThread ) {
            std::thread t1([&] {
                parallelMergeSort(arr, first, middle);
                });
            parallelMergeSort(arr, middle + 1, last);
            t1.join();
        }
        else {
            mergeSort(arr, first, middle);
            mergeSort(arr, middle + 1, last);
        }

        merge(arr, first, last);
    }
}


void mergeSortOMP(int* arr, long long first, long long last) {
    if (first < last) {
        int middle = first + (last - first) / 2;

        // Parallelize the two recursive calls
#pragma omp parallel sections
        {
#pragma omp section
            mergeSort(arr, first, middle);
#pragma omp section
            mergeSort(arr, middle + 1, last);
        }

        merge(arr, first, last);
    }
}

void quickSort(int* arr, int left, int right) {
    int i = left, j = right;
    int pivot = arr[(left + right) / 2];

    while (i <= j) {
        while (arr[i] < pivot)
            i++;
        while (arr[j] > pivot)
            j--;
        if (i <= j) {
            std::swap(arr[i], arr[j]);
            i++;
            j--;
        }
    }

    if (left < j)
        quickSort(arr, left, j);
    if (i < right)
        quickSort(arr, i, right);
}

long long partition(int* arr, long long low, long long high) {
    long long pivot = arr[high];
    long long i = (low - 1);

    for (long long j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[high]);
    return (i + 1);
}


void quickSortSingleThread(int* arr, long long first, long long last) {
    if (first < last) {
        long long pivotIndex = partition(arr, first, last);
        quickSortSingleThread(arr, first, pivotIndex - 1);
        quickSortSingleThread(arr, pivotIndex + 1, last);
    }
}

void quickSortParallelWithLimitedThreads(int* arr, long long left, long long right) {

        if (left >= right) return;

        long long leftBound = left;
        long long rightBound = right;
        long long middle = arr[(leftBound + rightBound) / 2];

        do {
            while (arr[leftBound] < middle) {
                leftBound++;
            }
            while (arr[rightBound] > middle) {
                rightBound--;
            }
            if (leftBound <= rightBound) {
                std::swap(arr[leftBound], arr[rightBound]);
                leftBound++;
                rightBound--;
            }
        } while (leftBound <= rightBound);

        if ((rightBound - left) > 10000) {
            std::future<void> f = std::async(std::launch::async, [&] {
                quickSort(arr, left, rightBound);
                });
            quickSort(arr, leftBound, right);
        }
        else {
            quickSort(arr, left, rightBound);
            quickSort(arr, leftBound, right);
        }
    }


bool sortedArray(int* arr, long long size) {
    for (long long i = 0; i < size - 1; i++) {
        if (arr[i] > arr[i + 1]) {
            return false;
        }
    }
    return true;
}

double durationMeasure(void (*sortFunction)(int*, long long, long long), long long arraySize) {
    int* arr = new int[arraySize];
    srand(0);

    for (long long i = 0; i < arraySize; i++)
        arr[i] = rand() % (LIMIT);

    double start = omp_get_wtime();
    sortFunction(arr, 0, arraySize - 1);
    double end = omp_get_wtime();

    if (!sortedArray(arr, arraySize))
        std::cout << "!!! Sorting failed !!!"; // Sort check

    delete[] arr;
    return end - start;
}

void presentDuration(std::string name, double basis, double duration) {
    std::cout << std::setw(40) << name << std::fixed << std::setprecision(4) << std::setw(10) << duration << " seconds  " << std::string(DOT_LENGTH * duration / basis, '.') << std::endl;
}

void runSortingTests() {
    std::vector<long long> arraySizes = { 100000, 400000, 1600000, 3200000};
    int numTests = arraySizes.size(); // Number of tests for each sorting method
    double d[FUNCTIONS_NUMBER];

    for (long long arraySize : arraySizes) {
        std::cout << "\nArray size: " << arraySize << "\n";

        d[MERGE_SORT] = durationMeasure(mergeSort, arraySize);
        d[PARALLEL_MERGE_SORT] = durationMeasure(parallelMergeSort, arraySize);
        d[PARALLEL_MERGE_SORT_OMP] = durationMeasure(mergeSortOMP, arraySize);
        d[PARALLEL_MERGE_SORT_THREAD] = durationMeasure(parallelMergeSortWithThread, arraySize);
        d[QUICK_SORT_SINGLE_THREAD] = durationMeasure(quickSortSingleThread, arraySize);
        d[QUICK_SORT_PARALLEL_WITH_FUTURE] = durationMeasure(quickSortParallelWithLimitedThreads, arraySize);

        double basis = 0;
        for (auto e : d) {
            basis = basis < e ? e : basis;
        }

        presentDuration("Single thread merge sort", basis, d[MERGE_SORT]);
        presentDuration("Single thread quick sort", basis, d[QUICK_SORT_SINGLE_THREAD]);
        presentDuration("Parallel merge sort with future", basis, d[PARALLEL_MERGE_SORT]);
        presentDuration("Parallel merge sort with thread", basis, d[PARALLEL_MERGE_SORT_THREAD]);
        presentDuration("Parallel merge sort with OMP", basis, d[PARALLEL_MERGE_SORT_OMP]);
        presentDuration("Parallel quick sort with future", basis, d[QUICK_SORT_PARALLEL_WITH_FUTURE]);
    }
}

int main() {
    std::cout.imbue(std::locale("en_US"));
    std::cout << "Sorting by different ways\n";

    //int coreNumber = omp_get_max_threads();
    //std::cout << coreNumber << " cores available\n" << std::endl;

    runSortingTests();

    return 0;
}
