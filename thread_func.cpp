#include <sys/sysinfo.h>
#include <string.h>
#include <sched.h>
#include "inc.h"

void* thread_func(void* ptr) {
    Args* ap = (Args*)ptr;
    double* a = ap->a;
    double* b = ap->b;
    double* x = ap->x;
    
    int n = ap->n;
    int m = ap->m;
    int p = ap->p;
    int k = ap->k;
    int s = ap->s;
    int r = ap->r;
    std::string name = ap->filename;
    
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus - 1 - (k % n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();

    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);

    for (int j = k*m; j < n; j += p*m) {
        int h = (j + m < n) ? m : n - j;
        for (int i = 0; i < n; ++i) {
            memset(a + i*n + j, 0, h*sizeof(double));   
        }

        memset(b + j, 0, h*sizeof(double));
        memset(x + j, 0, h*sizeof(double));
    }

    reduce_sum<int>(p);
    double* block1 = new double[m*m];
    double* block2 = new double[m*m];
    double* block3 = new double[m*m];
    double* main_col = new double[n*m]; 
    double* vec = new double[2];

    if (!name.empty()) {
        int res = 0;
        if (k == 0) { 
            res = read_array(a, n, name);
        }

        reduce_sum(p, &res, 1, &sum);
        if (res < 0) { 
            ap->res = res; 
            delete[] block1;
            delete[] block2;
            delete[] block3;
            delete[] main_col;
            delete[] vec;
            return nullptr;
        }

        init_b(a, b, n, m, p, k);

    } else {
        init_arrays(a, b, n, m, s, p, k);
    }

    double a_norm = -1;
    if (k == 0) {
        std::cout << "A:\n";
        output(a, n, r, n);
        std::cout << "\n";
        a_norm = matrix_norm(n, n, a);
    }
    reduce_sum(p, &a_norm, 1, &maximum);

    ap->t = get_cpu_time();

    int f = n / m;
    int l = n - f * m;
    int h = l ? f + 1 : f;

    for (int t = 0; t < f; ++t) {
        copy_block_col(a, main_col, n, m, t);
        double max_norm_block = -1;
        int row_max_block = -1;
        for (int i = t + k; i < f; i += p) { 
            get_block(i, 0, m, m, f, l, main_col, block1);
            if (is_inv(m, block1, a_norm)) {
                double block_norm = matrix_norm(m, m, block1);
                if (block_norm - EPS * a_norm > max_norm_block) {
                    max_norm_block = block_norm;
                    row_max_block = i;
                }
            }
        }

        vec[0] = max_norm_block;
        vec[1] = row_max_block;
        reduce_sum(p, vec, 2, &max); 
        max_norm_block = vec[0];
        row_max_block = vec[1];
        if (row_max_block == -1) {
            ap->method_not_applicable = true; 
            delete[] block1;
            delete[] block2;
            delete[] block3;
            delete[] main_col;
            delete[] vec;
            return nullptr;
        }

        for (int u = 0; u < m; ++u) {
            for (int v = 0; v < m; ++v) {
                std::swap(main_col[m*(m*row_max_block + u) + v], main_col[m*(m*t + u) + v]);
            }
        }

        swap_block_rows(a, b, row_max_block, t, n, m, f, l, h, p, k); 
        get_block(t, 0, m, m, f, l, main_col, block1);     
        inverse_matrix(m, block1, block2, a_norm); 

        int main_thread = (t + 1) % p;
        int q = (k >= main_thread) ? t + 1 + k - main_thread : t + 1 + k - main_thread + p;
        for (; q < f; q += p) { 
            get_block(t, q, n, m, f, l, a, block1); 
            matrix_product(m, m, m, block2, block1, block3); 
            put_block(t, q, n, m, f, l, block3, a);
        }

        if (l && k == f % p) {
            get_block(t, f, n, m, f, l, a, block1);
            matrix_product(m, m, l, block2, block1, block3);
            put_block(t, f, n, m, f, l, block3, a);
        }

        if (k == t % p) {
            matrix_product(m, m, 1, block2, b + m * t, block3); 
            put_vector(t, m, f, l, block3, b);
        }
        reduce_sum<int>(p); 

        for (q = t + 1; q < h; ++q) { 
            get_block(q, t, n, m, f, l, a, block1); 
            int multiplier_rows = q < f ? m : l;
            main_thread = (t + 1) % p;
            int v = (k >= main_thread) ? t + 1 + k - main_thread : t + 1 + k - main_thread + p;
            for (; v < h; v += p) { 
                get_block(t, v, n, m, f, l, a, block2); 
                int block_cols = v < f ? m : l;
                matrix_product(multiplier_rows, m, block_cols, block1, block2, block3);
                get_block(q, v, n, m, f, l, a, block2);
                subtract_matrix_inplace(multiplier_rows, block_cols, block2, block3);
                put_block(q, v, n, m, f, l, block2, a);
            }

            if (k == t % p) {
                matrix_product(multiplier_rows, m, 1, block1, b + m*t, block3);
                subtract_matrix_inplace(1, multiplier_rows, b + m*q, block3);
            }

            reduce_sum<int>(p);
        }
    }
    
    if (l) {
        int res = 0;
        if (k == f % p) {
            get_block(f, f, n, m, f, l, a, block1);  
            if (inverse_matrix(l, block1, block2, a_norm)) {
                matrix_product(l, l, 1, block2, b + m*f, block3);
                put_vector(f, m, f, l, block3, x);
            } else {
                res = -1;
            }
        }

        reduce_sum(p, &res, 1, &sum);
        if (res < 0) { 
            ap->method_not_applicable = true; 
            delete[] block1;
            delete[] block2;
            delete[] block3;
            delete[] main_col;
            delete[] vec;
            return nullptr;
        }
    }

    for (int i = f - 1; i >= 0; --i) {
        for (int vv = 0; vv < m; ++vv) { 
            block2[vv] = 0;
        }

        int main_thread = (i + 1) % p;
        int j = (k >= main_thread) ? i + 1 + k - main_thread : i + 1 + k - main_thread + p;
        for (; j < h; j += p) { 
            get_block(i, j, n, m, f, l, a, block1);
            int cols = j < f ? m : l;
            matrix_product(m, cols, 1, block1, x + m*j, block3);

            for (int p = 0; p < m; ++p) {
                block2[p] += block3[p];
            }
        }

        reduce_sum(p, block2, m, &sum);

        if (k == i % p) {
            for (int p = 0; p < m; ++p) {
                b[m*i + p] -= block2[p]; 
            }

            put_vector(i, m, f, l, b + m*i, x);
        }
        
        reduce_sum<int>(p);
    }
    

    delete[] block1;
    delete[] block2;
    delete[] block3;
    delete[] main_col;
    delete[] vec;

    reduce_sum<int>(p);
    return nullptr;
}
