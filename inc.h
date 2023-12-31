#define EPS 1e-14
#pragma once
#include <iostream>
#include <string>
#include <pthread.h>

void init_b(double* matrix, double* b, int n, int m, int p, int k) ;
void init_arrays(double* matrix, double* b, int n, int m, int s, int p, int k);
int read_array(double* array, int n, const std::string& name_file);
void output(double* array, int n, int r, int l);
void print_block(double* block, int m);
void get_block(int i, int j, int n, int m, int f, int l, double* matrix, double* block1);
void put_block(int i, int j, int n, int m, int k, int l, double* block, double* matrix);
void put_vector(int i, int m, int k, int l, double* b_i, double* b);
double matrix_norm(int n, int m, double* matrix);
void swap_rows(double* matrix, int n, int i, int j);
void swap_block_rows(double* a, double* b, int row_max_block, int t, int n, int m, int f, int l, int h, int p, int k);
bool is_inv(int m, double* matrix, double a_norm);
bool inverse_matrix(int m, double* matrix, double* identity, double a_norm);
void copy_block_col(double* a, double* main_col, int n, int m, int t);
void matrix_product(int n, int m, int k, double* a, double* b, double* c);
void subtract_matrix_inplace(int n, int m, double* a, double* b);
double r1_eval(int n, double* matrix, double* x, double* b, double* c, double* d);
double r2_eval(int n, double* x);
void subtract_vector(int n, double* a, double* b, double* c);
double vector_norm(int n, double* vec);
void simple_input_matrix(int s, int n, double* matrix);
void simple_input_b(int n, double* matrix, double* b);
double get_cpu_time();
void* thread_func(void* ptr);

class Args {
public:
    double* a = nullptr;    
    double* b = nullptr;    
    double* x = nullptr;      

    double t = 0;
    int n = 0;
    int m = 0;
    int p = 0; 
    int k = 0; 
    int r = 0;
    int s = 0; 
    int res = 0; 
    bool method_not_applicable = false;
    std::string filename;

    Args() = default;

    Args(const Args& x) = delete;
    Args(Args&& x);
    Args& operator=(const Args& x) = delete;
    Args& operator=(Args&& x);

    ~Args();
};

template<class T>
void sum(T* r, T* a, int n) {
    for (int i = 0; i < n; ++i) {
        r[i] += a[i];
    }    
}

template<class T>
void min(T* r, T* a, int n) {
    for (int i = 0; i < n; ++i) {
        r[i] = std::min(r[i], a[i]);
    }
}

template<class T>
void maximum(T* r, T* a, int n) {
    for (int i = 0; i < n; ++i) {
        r[i] = std::max(r[i], a[i]);
    }
}

template<class T>
void max(T* r, T* a, int n) {
    if (n != 2) {
        std::cout << "Ты не ту функцию передал в reduce sum" << "\n";
        return;
    }
    
    if (a[0] > r[0] + EPS) {
        r[0] = a[0];
        r[1] = a[1]; 
    }
}

template<class T>
void reduce_sum(int p, T* a = nullptr, int n = 0, void (*func)(T*, T*, int) = &sum) {
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static T* r = nullptr;
    int i;

    if (p <= 1) {
        return;
    }

    pthread_mutex_lock(&m);

    if (r == nullptr) {
        r = a;
    } else {
        func(r, a, n);
    }

    t_in++;
    if (t_in >= p) {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
        
    } else {
        while (t_in < p) {
            pthread_cond_wait(&c_in, &m);
        }
    }

    if (r != a) {
        for (i = 0; i < n; ++i) {
            a[i] = r[i];
        }
    } 

    t_out++;
    if (t_out >= p) {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);        
    } else {
        while (t_out < p) {
            pthread_cond_wait(&c_out, &m);
        }
    }

    pthread_mutex_unlock(&m);
}
