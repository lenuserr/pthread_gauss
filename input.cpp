#include <fstream>
#include <cmath>
#include "inc.h"

double f(int s, int n, int i, int j) { 
    switch(s) {
        case 1:
            return n - std::max(i, j);
        case 2:
            return std::max(i, j) + 1;
        case 3:
            return std::fabs(i-j);
        case 4:
            return 1. / (i + j + 1);
        default:
            return 1;
    }
}

void init_arrays(double* matrix, double* b, int n, int m, int s, int p, int k) {
    for (int j = k*m; j < n; j += p*m) {
        int h = (j + m < n) ? m : n - j;
        for (int i = 0; i < n; ++i) {
            for (int v = 0; v < h; ++v) {
                matrix[i*n + j + v] = f(s, n, i, j + v);
            }   
        }
    }

    reduce_sum<int>(p); 
    for (int i = k*m; i < n; i += p*m) { 
        int h = (i + m < n) ? m : n - i;
        for (int t = i; t < i + h; ++t) {
            double sum = 0;
            for (int k = 0; k <= (n-1) / 2; ++k) {
                sum += matrix[n * t + 2*k];
            }
            b[t] = sum; 
        }
    }   

    reduce_sum<int>(p);
}

void init_b(double* matrix, double* b, int n, int m, int p, int k) {
    for (int i = k*m; i < n; i += p*m) { 
        int h = (i + m < n) ? m : n - i;
        for (int t = i; t < i + h; ++t) {
            double sum = 0;
            for (int k = 0; k <= (n-1) / 2; ++k) {
                sum += matrix[n * t + 2*k];
            }
            b[t] = sum; 
        }
    }

    reduce_sum<int>(p);    
}

void simple_input_matrix(int s, int n, double* matrix) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[n * i + j] = f(s, n, i, j);
        }
    }       
}

void simple_input_b(int n, double* matrix, double* b) {
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int k = 0; k <= (n-1) / 2; ++k) {
            sum += matrix[n * i + 2*k];
        }
        b[i] = sum; 
    } 
}

int read_array(double* array, int n, const std::string& name_file) {
    std::ifstream fin(name_file);
    
    if (!fin.is_open()) {
        return -1;
    }
    
    double x;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(fin >> x)) {
                return -2;  
            }
            
            array[n * i + j] = x;
        }
    }

    return 0;
}
