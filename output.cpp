#include <iostream>

void output(double* array, int n, int r, int l) {
    int h = std::min(l, r);
    int w = std::min(n, r);
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            printf("%10.3e ", array[n*i + j]);
        }
        printf("\n");
    }
}
