// #############################################################################
// Provides quicksortMap(...), which uses quicksort to return a permutation map 
// corresponding to an ordered input array
// #############################################################################
#pragma once

namespace sorts{

    // Anonymous namespace creates function only accessible to the sorts namespace
    namespace {
        int partitionMap(double *arr, int *map, int lo, int hi) {
            int tmp;
            double pivot = arr[map[hi]];

            int i = lo - 1;

            for (int j=lo; j<hi; j++) {
                if (arr[map[j]] <= pivot) {
                    i++;
                    tmp = map[i];
                    map[i] = map[j];
                    map[j] = tmp;
                }
            }

            i++;
            tmp = map[i];
            map[i] = map[hi];
            map[hi] = tmp;

            return i;
        }
    }

    void quicksortMap(double *arr, int *map, int lo, int hi) 
    {
        if ((lo >= hi) || (lo < 0) ) return;

        int p = partitionMap(arr, map, lo, hi);

        quicksortMap(arr, map, lo, p-1);
        quicksortMap(arr, map, p+1, hi);
    }

}