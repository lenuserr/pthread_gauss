#include <sys/sysinfo.h>
#include <string.h>
#include <sched.h>
#include "inc.h"

// 15.11.23. 13:00.
// Переписываю код по столбцам (ровно так как у меня в отчёте в общем)
// 14:00. Уже полным ходом пишу выбор главного элемента, затем прямой ход. Work, work.
// Сегодня-завтра должен всё написать (даже с учетом просмотра лекции по град.бустингу today).
// Просто работай и всё получится.
// 15:15. Выбор главного элемента сделал. Переставляем строки как написано в отчете.
// 15:50. Переставляю строки. Ещё шуть-шуть. Прямой метод Гаусса должен до 19 00 кончить с едой.
// 16:00. Строки-то я поменял. Теперь надо подумать как менять у вектора b строки? Ну пока походу в тупую. Двумя потоками как бгч говорил крч.
// 17:00. Какой-то пиздец.
// 17:15. Бля, я даун, всё правильно работало. Идём дальше.

void* thread_func(void* ptr) {
    [[maybe_unused]] static pthread_mutex_t M = PTHREAD_MUTEX_INITIALIZER;
    
    Args* ap = (Args*)ptr;
    double* a = ap->a;
    double* b = ap->b;
    double* x = ap->x;
    
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

    for (int j = k*m; j < n; j += p*m) {
        int h = (j + m < n) ? m : n - j;
        for (int i = 0; i < n; ++i) {
            memset(a + i*n + j, 0, h*sizeof(double));   
        }

        memset(b + j, 0, h*sizeof(double));
        memset(x + j, 0, h*sizeof(double));
    }

    reduce_sum<int>(p);
    // здесь выделение локальных буферных массивов вспомогательных.
    double* block1 = new double[m*m];
    double* main_col = new double[n*m]; // копируем t-ый столбец сюда в каждом потоке типа.
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
        //for (int i = 0; i < n; ++i) { std::cout << b[i] << " "; }
        //std::cout << "\n\n";
        /*
        std::cout << "A:\n";
        output(a, n, r, n);
        std::cout << "\n";
        std::cout << "b:\n";
        output(b, n, r, 1);
        std::cout << "\n";
        */
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
    [[maybe_unused]] int l = n - f * m;
    [[maybe_unused]] int h = l ? f + 1 : f;

    for (int t = 0; t < f; ++t) {
        // скопировать t-ый столбец в каждый поток. (в main_col)
        copy_block_col(a, main_col, n, m, t);
        double max_norm_block = -1;
        int row_max_block = -1;
        for (int i = t + k; i < f; i += p) {
            get_block(i, 0, m, m, f, l, main_col, block1);
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

        // в каждом потоке поменяем main_col на тот, который теперь будет после перестановки строк. 
        for (int u = 0; u < m; ++u) {
            for (int v = 0; v < m; ++v) {
                std::swap(main_col[m*(m*row_max_block + u) + v], main_col[m*(m*t + u) + v]);
            }
        }

        // row_max_block нужно переставить с t-ой строкой.
        swap_block_rows(a, b, row_max_block, t, n, m, f, l, h, p, k);   
        // переходим к пункту 3.2 из отчёта.     
    }

    delete[] block1;
    delete[] main_col;
    delete[] vec;

    return nullptr;
}
