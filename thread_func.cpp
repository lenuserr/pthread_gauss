#include <iostream>
#include <cmath>
#include <sys/sysinfo.h>
#include <string.h>
#include <sched.h>
#include "inc.h"

void* thread_func(void* ptr) {
    static pthread_mutex_t M = PTHREAD_MUTEX_INITIALIZER;
    Args* ap = (Args*)ptr;
    double* a = ap->a;
    double* b = ap->b;
    double* c = ap->c;

    int n = ap->n;
    int m = ap->m;
    int p = ap->p;
    int k = ap->k;
    int s = ap->s;
    //int r = ap->r;
    std::string name = ap->filename;

    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus - 1 - (k % n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();

    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
    
    int i;
    for (i = k*m; i < n; i += p*m) { // работаем со "своими" строками
        int h = (i + m < n) ? m : n - i; 
        memset(a + i*n, 0, h*n*sizeof(double));
        memset(b + i, 0, h*sizeof(double));
        memset(c + i, 0, h*sizeof(double));
    }

    reduce_sum<int>(p);

    // здесь выделение локальных буферных массивов вспомогательных.
    //double* block1 = new double[n]; // для копирования "главной" строки во все потоки.
    double* block = new double[m*m];
    int* rows = new int[m];
    double* vec = new double[2];

    if (!name.empty()) {
        int res = 0;
        if (k == 0) { // главный поток
            res = read_array(a, n, name);
        }

        reduce_sum(p, &res, 1);

        if (res < 0) { // такой res во всех потоках будет, т.к. сделали reduce sum.
            // освободить память и всё прочее в main.
            ap->res = res; // по этому флагу поймём в main всё норм отработало или нет.
            // этот res проверим после pthread_join всех потоков в мэйне.
            return nullptr;
        }

        init_b(a, b, n, m, p, k);

    } else {
        init_matrix(a, n, m, s, p, k);
        init_b(a, b, n, m, p, k);
    }

    if (k == 0) {
        //output(a, n, r, n);      
        //output(b, 1, r, n);  
    }
    reduce_sum<int>(p);
    //t = get_cpu_time();
    // solve(a, b, n, p, m, k);
    // t = get_cpu_time() - t;
    /*
    pthread_mutex_lock(&M);
    std::cout << "РАБОТАЕТ ПОТОК ПОД НОМЕРОМ " << k << "\n\n\n";
    int f = n / m;
    for (int t = 0; t < f; ++t) {
        std::cout << "ШАГ ПОД НОМЕРОМ " << t << "\n";
        int q = t % p;
        i = (k - q) >= 0 ? (t + k - q) * m : (t + k + p - q) * m;
        for (; i < n; i += p*m) {
            int h = (i + m < n) ? m : n - i;
            for (int u = i; u < i + h; ++u) {
                for (int j = t*m; j < n; ++j) {
                    std::cout << a[n*u + j] << " ";                    
                }
                std::cout << "\n";
            }                        
        }
    } 
    std::cout << "ПОТОК ПОД НОМЕРОМ " << k << " ЗАКОНЧИЛ РАБОТУ" << "\n\n\n";
    pthread_mutex_unlock(&M);
    */

    int f = n / m;
    int l = n - f * m;
    for (int t = 0; t < f; ++t) {
        int q = t % p;
        i = (k - q) >= 0 ? (t + k - q) : (t + k + p - q);
        double max_norm_block = -1;
        int row_max_block = -1;
        for (; i < f; i += p) {
            get_block(i, t, n, m, f, l, a, block);
            double block_norm = matrix_norm(m, m, block);
            // block_norm - EPS * a_norm надо. потом посчитаю a_norm, пока не надо и лень.
            if (block_norm - EPS > max_norm_block && is_inv(m, block, 1, rows)) {
                max_norm_block = block_norm;
                row_max_block = i;
            }
            /*
            for (int v = 0; v < m; ++v) {
                for (int z = 0; z < m; ++z) {
                    std::cout << block[m*v + z] << " ";  
                }
                std::cout << "\n";
            }
            */
        }

        pthread_mutex_lock(&M);
        std::cout << "ШАГ " << t << ". РАБОТАЕТ ПОТОК ПОД НОМЕРОМ " << k << "\n";
        std::cout << "Local max norm block: " << max_norm_block << " Row: " << row_max_block << "\n\n\n";
        pthread_mutex_unlock(&M);
        vec[0] = max_norm_block;
        vec[1] = row_max_block;
        reduce_sum(p, vec, 2, &max); //чтобы выбрать главный элемент по столбцу.
        max_norm_block = vec[0];
        row_max_block = vec[1];
        pthread_mutex_lock(&M);
        std::cout << "ШАГ " << t << ". РАБОТАЕТ ПОТОК ПОД НОМЕРОМ " << k << "\n";
        std::cout << "Global max norm block: " << max_norm_block << " Row: " << row_max_block << "\n\n\n";
        pthread_mutex_unlock(&M);
    } 

    delete[] block;
    delete[] rows;
    delete[] vec;
    // delete[] block1;

    return nullptr;
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

bool is_inv(int m, double* matrix, double a_norm, int* rows) {
    for (int k = 0; k < m; ++k) {
        rows[k] = k;
    }

    for (int i = 0; i < m; ++i) {
        double max_elem = matrix[rows[i] * m + i];
        int row_max_elem = i;
        for (int j = i + 1; j < m; ++j) {
            if (std::fabs(matrix[rows[j] * m + i]) - EPS * a_norm > std::fabs(max_elem)) {
                max_elem = matrix[rows[j] * m + i];
                row_max_elem = j;
            }
        }

        std::swap(rows[i], rows[row_max_elem]);

        if (std::fabs(max_elem) < EPS * a_norm) {
            return false;    
        }

        double factor = 1 / max_elem;
        for (int s = i; s < m; ++s) {
            matrix[rows[i] * m + s] *= factor;
        }

        for (int k = i + 1; k < m; ++k) {
            double multiplier = -matrix[rows[k] * m + i];
            for (int p = i + 1; p < m; ++p) { 
                matrix[rows[k] * m + p] += matrix[rows[i] * m + p] * multiplier;
            }
        }
    }

    return true;
}
