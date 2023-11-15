#define EPS 1e-14
#include <cmath>
#include "inc.h"

void swap_rows(double* matrix, int n, int i, int j) {
    for (int k = 0; k < n; ++k) {
        std::swap(matrix[n*i + k], matrix[n*j + k]);
    }
}

void copy_block_col(double* a, double* main_col, int n, int m, int t) {
    for (int i = t*m; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            main_col[m*i + j] = a[n*i + t*m + j];
        }
    }      
}

void swap_block_rows(double* a, double* b, int row_max_block, int t, int n, int m, int f, int l, int h, int p, int k) {
    for (int j = t + k; j < h; j += p) {
        int v_max = (j < f) ? m : l;
        for (int u = 0; u < m; ++u) {
            for (int v = 0; v < v_max; ++v) {
                std::swap(a[n * (m * row_max_block + u) + m * j + v], a[n * (m * t + u) + m * j + v]);
            }
        }
    }

    // переставляю две блочные строки у вектора b.
    int k1 = t % p;
    int k2 = row_max_block % p;
    if (k == k1) {
        for (int u = 0; u < m / 2; u++) {
            std::swap(b[m*row_max_block + u], b[m*t + u]);    
        }
    }

    if (k == k2) {
        for (int u = m / 2; u < m; u++) {
            std::swap(b[m*row_max_block + u], b[m*t + u]);    
        }
    }

    reduce_sum<int>(p); // вторая точка синхронизации в алгоритме (по отчёту).
}

void get_block(int i, int j, int n, int m, int f, int l, double* matrix, double* block1) {
    int h = i < f ? m : l;
    int w = j < f ? m : l;
    
    int ind = 0;

    for (int p = 0; p < h; ++p) {
        for (int q = 0; q < w; ++q) {
            block1[ind] = matrix[n * (m * i + p) + m * j + q];
            ind++;
        }
    }
}

void print_block(double* block, int m) {
    for (int u = 0; u < m; ++u) {
        for (int v = 0; v < m; ++v) {
            std::cout << block[m*u + v] << " ";
        }
        std::cout << "\n";
    }
}

void put_block(int i, int j, int n, int m, int k, int l, double* block, double* matrix) {
    int h = i < k ? m : l;
    int w = j < k ? m : l;

    int ind = 0;
    for (int p = 0; p < h; ++p) {
        for (int q = 0; q < w; ++q) {
            matrix[n * (m * i + p) + m * j + q] = block[ind];
            ind++;
        }
    }
}

void put_vector(int i, int m, int k, int l, double* b_i, double* b) {
    int length = i < k ? m : l;
    for (int p = 0; p < length; ++p) {
        b[m*i + p] = b_i[p];
    }
}

double matrix_norm(int n, int m, double* matrix) {
    double norm = -1;
    for (int j = 0; j < m; ++j) {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += std::fabs(matrix[m*i + j]);           
        }

        norm = std::max(norm, sum);
    }

    return norm;
}

void subtract_matrix_inplace(int n, int m, double* a, double* b) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            a[m*i + j] -= b[m*i + j];           
        }
    }
}

bool inverse_matrix(int m, double* matrix, double* identity, double a_norm) {

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            identity[i * m + j] = (i != j) ? 0 : 1;
        }
    }

    for (int i = 0; i < m; ++i) {
        double max_elem = matrix[i * m + i];
        int row_max_elem = i;
        for (int j = i + 1; j < m; ++j) {
            if (std::fabs(matrix[j * m + i]) > std::fabs(max_elem)) {
                max_elem = matrix[j * m + i];
                row_max_elem = j;
            }
        }

        swap_rows(matrix, m, i, row_max_elem);
        swap_rows(identity, m, i, row_max_elem);

        if (std::fabs(max_elem) < EPS * a_norm) {
            return false;    
        }

        double factor = 1 / max_elem;
        for (int s = 0; s < i; ++s) {
            identity[i * m + s] *= factor;
        }

        for (int s = i; s < m; ++s) {
            matrix[i * m + s] *= factor;
            identity[i * m + s] *= factor;
        }

        for (int k = i + 1; k < m; ++k) {
            double multiplier = -matrix[k * m + i];
            for (int p = 0; p < i + 1; ++p) {
                identity[k * m + p] += identity[i * m + p] * multiplier;
            }

            for (int p = i + 1; p < m; ++p) { 
                matrix[k * m + p] += matrix[i * m + p] * multiplier;
                identity[k * m + p] += identity[i * m + p] * multiplier;
            }
        }
    }

    for (int i = m - 1; i > 0; --i) {
        for (int k = i - 1; k >= 0; --k) {
            double multiplier = -matrix[k * m + i];
            for (int p = 0; p < m; ++p) { 
                identity[k * m + p] += identity[i * m + p] * multiplier;
            }
        }
    }

    return true;
}

bool is_inv(int m, double* matrix, double a_norm) {
    for (int i = 0; i < m; ++i) {
        double max_elem = matrix[i * m + i];
        int row_max_elem = i;
        for (int j = i + 1; j < m; ++j) {
            if (std::fabs(matrix[j * m + i]) - EPS * a_norm > std::fabs(max_elem)) {
                max_elem = matrix[j * m + i];
                row_max_elem = j;
            }
        }

        swap_rows(matrix, m, i, row_max_elem);

        if (std::fabs(max_elem) < EPS * a_norm) {
            return false;    
        }

        double factor = 1 / max_elem;
        for (int s = i; s < m; ++s) {
            matrix[i * m + s] *= factor;
        }

        for (int k = i + 1; k < m; ++k) {
            double multiplier = -matrix[k * m + i];
            for (int p = i + 1; p < m; ++p) { 
                matrix[k * m + p] += matrix[i * m + p] * multiplier;
            }
        }
    }

    return true;
}

void matrix_product(int n, int m, int k, double* a, double* b, double* c) {
    
    for (int i = 0; i < n*k; ++i) {
        c[i] = 0;
    }
    
    double sum00, sum01, sum02, sum10, sum11, sum12, sum20, sum21, sum22;
    
    int v3 = n%3;
    int h3 = k%3;
    
    for (int i = 0; i < v3; ++i) {
        for (int j = 0; j < h3; ++j) {
            sum00 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*i + p] * b[k*p + j];
            }
            
            c[k*i + j] = sum00;
        }
        for (int j = h3; j < k; j+=3) {
            sum00 = 0; sum01 = 0; sum02 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*i + p] * b[k*p + j];
                sum01 += a[m*i + p] * b[k*p + j + 1];
                sum02 += a[m*i + p] * b[k*p + j + 2];
            }
            c[k*i + j] = sum00;
            c[k*i + j + 1] = sum01;
            c[k*i + j + 2] = sum02;
        }
    }
    
    for (int i = v3; i < n; i+=3) {   
        for (int j = 0; j < h3; ++j) {
            sum00 = 0; sum01 = 0; sum02 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*i + p] * b[k*p + j];
                sum01 += a[m*(i + 1) + p] * b[k*p + j];
                sum02 += a[m*(i + 2) + p] * b[k*p + j];
            }
            c[k*i + j] = sum00;
            c[k*(i + 1) + j] = sum01;
            c[k*(i + 2) + j] = sum02;
        }
        
        for (int j = h3; j < k; j+=3) {
            sum00 = 0; sum01 = 0; sum02 = 0; sum10 = 0; sum11 = 0; sum12 = 0;
            sum20 = 0; sum21 = 0; sum22 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*i + p] * b[k*p + j];
                sum01 += a[m*i + p] * b[k*p + j + 1];
                sum02 += a[m*i + p] * b[k*p + j + 2];
                sum10 += a[m*(i + 1) + p] * b[k*p + j];
                sum11 += a[m*(i + 1) + p] * b[k*p + j + 1];
                sum12 += a[m*(i + 1) + p] * b[k*p + j + 2];
                sum20 += a[m*(i + 2) + p] * b[k*p + j];
                sum21 += a[m*(i + 2) + p] * b[k*p + j + 1];
                sum22 += a[m*(i + 2) + p] * b[k*p + j + 2];
            }
            c[k*i + j] = sum00;
            c[k*i + j + 1] = sum01; 
            c[k*i + j + 2] = sum02;
            c[k*(i + 1) + j] = sum10;
            c[k*(i + 1)  + j + 1] = sum11;
            c[k*(i + 1) + j + 2] = sum12;
            c[k*(i + 2) + j] = sum20;
            c[k*(i + 2) + j + 1] = sum21;
            c[k*(i + 2) + j + 2] = sum22;
        }
    }    
}
