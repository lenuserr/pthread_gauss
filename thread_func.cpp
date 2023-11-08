#include <sys/sysinfo.h>
#include <string.h>
#include <sched.h>
#include "inc.h"

// 15:15. 08.11.23. let's go code pthread gauss.
// 18:50. напиши хоть как-то крч для начала. хотя я и так уже это 3 часа делаю)
// потом надо на reduce sum проверить: недостающие поставить, лишние убрать, если такие будут.

void* thread_func(void* ptr) {
    //static pthread_mutex_t M = PTHREAD_MUTEX_INITIALIZER;
    Args* ap = (Args*)ptr;
    double* a = ap->a;
    double* b = ap->b;
    double* x = ap->x;

    int n = ap->n;
    int m = ap->m;
    int p = ap->p;
    int k = ap->k;
    int s = ap->s;
    [[maybe_unused]] int r = ap->r;
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
        memset(x + i, 0, h*sizeof(double));
    }

    reduce_sum<int>(p);

    // здесь выделение локальных буферных массивов вспомогательных.
    //double* block1 = new double[n]; // для копирования "главной" строки во все потоки.
    double* block1 = new double[m*m];
    double* block2 = new double[m*m];
    double* block3 = new double[m*m];
    double* main_row = new double[m*n];
    double* main_b = new double[m];
    double* vec = new double[2];

    if (!name.empty()) {
        int res = 0;
        if (k == 0) { // главный поток
            res = read_array(a, n, name);
        }

        reduce_sum(p, &res, 1, &sum);
        if (res < 0) { 
            ap->res = res; 
            delete[] block1;
            delete[] block2;
            delete[] block3;
            delete[] main_row;
            delete[] main_b;
            delete[] vec;
            return nullptr;
        }

        init_b(a, b, n, m, p, k);

    } else {
        init_arrays(a, b, n, m, s, p, k);
    }

    
    double a_norm = -1;
    if (k == 0) {
        //output(a, n, r, n);      
        a_norm = matrix_norm(n, n, a);
    }
    reduce_sum(p, &a_norm, 1, &maximum);

    /*
    t = get_cpu_time();
    solve(a, b, n, p, m, k);
    t = get_cpu_time() - t;
    */

    // вот она ниже моя solve функция пишется по сути. пока тут будет прост.
    int f = n / m;
    int l = n - f * m;
    int h = l ? f + 1 : f;

    for (int t = 0; t < f; ++t) {
        int q = t % p; // главная строка принадлежит q-ому потоку.
        i = (k - q) >= 0 ? (t + k - q) : (t + k + p - q);
        double max_norm_block = -1;
        int row_max_block = -1;
        for (; i < f; i += p) {
            get_block(i, t, n, m, f, l, a, block1);
            double block_norm = matrix_norm(m, m, block1);
            if (block_norm - EPS * a_norm > max_norm_block && is_inv(m, block1, a_norm)) {
                max_norm_block = block_norm;
                row_max_block = i;
            }
        }

        vec[0] = max_norm_block;
        vec[1] = row_max_block;
        reduce_sum(p, vec, 2, &max); // чтобы выбрать главный элемент по столбцу.
        max_norm_block = vec[0];
        row_max_block = vec[1];
        if (row_max_block == -1) {
            ap->method_not_applicable = true; // не применим, т.к. ни один поток не нашел обратный.
            delete[] block1;
            delete[] block2;
            delete[] block3;
            delete[] main_row;
            delete[] main_b;
            delete[] vec;
            return nullptr;
        }
        
        // у меня тупо один поток меняет строки и умножает на обратную, а все остальные потоки ждут...
        // пока хз как это переписывать, напишу остальную логику и потом вернусь сюда.
        if (k == q) {
            // t-ую блочную строку надо поменять с row_max_block-ой.
            swap_block_rows(a, t, row_max_block, 0, n, n, m);
            
            get_block(t, t, n, m, f, l, a, block1);
            inverse_matrix(m, block1, block2, a_norm); 
            
            for (int j = t + 1; j < f; ++j) { // j = t + 1 тут надо, чтобы более эффективно было!
                get_block(t, j, n, m, f, l, a, block1);
                matrix_product(m, m, m, block2, block1, block3); 
                put_block(t, j, n, m, f, l, block3, a);
            }       

            if (l) {
                get_block(t, f, n, m, f, l, a, block1);
                matrix_product(m, m, l, block2, block1, block3);
                put_block(t, f, n, m, f, l, block3, a); 
            }

            //std::cout << row_max_block << "\n";
            //output(a, n, n, n);
            //std::cout << "\n\n";     

            matrix_product(m, m, 1, block2, b + m * t, block3); 
            put_vector(t, m, f, l, block3, b);
        } 

        /*
        if (k == row_max_block % p) {
            swap_block_rows(a, t, row_max_block, n / 2, n, n, m);            
        } 
        */

        reduce_sum<int>(p);
        // копируем в каждый поток главную блочную строку и главную часть вектора b не забываем.
        copy_block_row(a, main_row, n, m, t);
        for (int v = 0; v < m; ++v) {
            main_b[v] = b[m*t + v];
        }

        reduce_sum<int>(p);

        i = (k - q) >= 0 ? (t + k - q) : (t + k + p - q);
        i = (i > t) ? i : i + p; // нижележащие ведь. на себя бы самого не попасться бы типа.
        for (; i < h; i += p) {
            get_block(i, t, n, m, f, l, a, block1);
            int multiplier_rows = i < f ? m : l;
            for (int r = t + 1; r < h; ++r) {
                get_block(0, r, n, m, f, l, main_row, block2); 
                int block_cols = r < f ? m : l;
                matrix_product(multiplier_rows, m, block_cols, block1, block2, block3);
                get_block(i, r, n, m, f, l, a, block2);
                subtract_matrix_inplace(multiplier_rows, block_cols, block2, block3);
                put_block(i, r, n, m, f, l, block2, a);
            }

            matrix_product(multiplier_rows, m, 1, block1, main_b, block3);
            subtract_matrix_inplace(1, multiplier_rows, b + m*i, block3);
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
            delete[] main_row;
            delete[] main_b;
            delete[] vec;
            return nullptr;
        }
    }

    // прямой ход метода Гаусса завершен.
    // остался обратный. Подумай, почитай идеи в отчете других(Влада) ес чо.

    delete[] block1;
    delete[] block2;
    delete[] block3;
    delete[] main_row;
    delete[] main_b;
    delete[] vec;

    reduce_sum<int>(p);
    return nullptr;
}
