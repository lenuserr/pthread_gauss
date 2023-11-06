#include <fstream>
#include <cmath>

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

void init_matrix(double* matrix, int n, int m, int s, int p, int k) {
    for (int i = k*m; i < n; i += p*m) { 
        int h = (i + m < n) ? m : n - i;
        for (int t = i; t < i + h; ++t) {
            for (int j = 0; j < n; ++j) {
                matrix[n * t + j] = f(s, n, t, j);
            }
        }
    }       
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
