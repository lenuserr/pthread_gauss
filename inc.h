#define EPS 1e-14
#pragma once
#include <string>
#include <pthread.h>

void init_matrix(double* matrix, int n, int m, int s, int p, int k);
void init_b(double* matrix, double* b, int n, int m, int p, int k);
int read_array(double* array, int n, const std::string& name_file);
void output(double* array, int n, int r, int l);
void get_block(int i, int j, int n, int m, int k, int l, double* matrix, double* block1);
double matrix_norm(int n, int m, double* matrix);
bool is_inv(int m, double* matrix, double a_norm, int* rows);
void* thread_func(void* ptr);

class Args {
public:
    double* a = nullptr;    
    double* b = nullptr;    
    double* c = nullptr;      

    int n = 0;
    int m = 0;
    int p = 0; // общее число потоков
    int k = 0; // номер потока
    int r = 0;
    int s = 0; // номер формулы для инициализации
    int res = 0; // результат чтения из файла
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
void max(T* r, T* a, int n) {
    if (n != 2) {
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
