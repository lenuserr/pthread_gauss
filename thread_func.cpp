#include <iostream>
#include <sys/sysinfo.h>
#include <string.h>
#include <sched.h>
#include "inc.h"
// вспоминаю как коммитить изменения на гитхаб
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
    double* block1 = new double[n]; // для копирования "главной" строки во все потоки.

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

    delete[] block1;

    return nullptr;
}
