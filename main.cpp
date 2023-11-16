#include <iostream>
#include "inc.h"

int main(int argc, char* argv[]) {
    const int task = 9;

    int n = std::stoi(argv[1]);
    int m = std::stoi(argv[2]);
    int p = std::stoi(argv[3]);
    int r = std::stoi(argv[4]);
    int s = std::stoi(argv[5]);
    int k = 0; 
    std::string filename;
    if (!s) {
        filename = argv[argc - 1];
    }   

    double* a = new double[n*n]; 
    double* b = new double[n];
    double* x = new double[n];

    Args* ap = nullptr;
    pthread_t* tid = nullptr;

    ap = new Args[p];
    tid = new pthread_t[p];
    
    for (k = 0; k < p; ++k) {
        ap[k].a = a; ap[k].b = b; ap[k].x = x;
        ap[k].n = n; ap[k].m = m; ap[k].p = p; ap[k].k = k;
        ap[k].r = r; ap[k].s = s; ap[k].filename = filename;
    }
    
    for (k = 1; k < p; ++k) {
        pthread_create(tid + k, 0, &thread_func, ap + k);
    }

    thread_func(ap + 0);

    for (k = 1; k < p; ++k) {
        pthread_join(tid[k], 0);
    }

    double t1 = get_cpu_time() - ap[0].t;

    if (ap[0].res < 0 || ap[0].method_not_applicable) {
        printf (
        "%s : Task = %d Res1 = %d Res2 = %d T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
        argv[0], task, -1, -1, t1, 0., s, n, m, p);
        delete[] a; delete[] b; delete[] x;
        delete[] ap; delete[] tid;   
        return -1;  
    }        

    std::cout << "x:\n";
    output(x, n, r, 1);
    std::cout << "\n";

    if (s) {
        simple_input_matrix(s, n, a);
    } else {
        std::string name_file = argv[argc - 1];
        if (read_array(a, n, name_file)) {
            std::cout << "Problems reading from a file " << filename << "\n";

            delete[] a; delete[] b; delete[] x;
            delete[] ap; delete[] tid;   
            return -2;
        }
    }

    simple_input_b(n, a, b);   

    double* c = new double[n];
    double* d = new double[n];
    double t2 = get_cpu_time();
    double r1 = r1_eval(n, a, x, b, c, d);
    double r2 = r2_eval(n, x);
    t2 = get_cpu_time() - t2;

    printf (
    "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
    argv[0], task, r1, r2, t1, t2, s, n, m, p);

    delete[] a; delete[] b; delete[] x;
    delete[] ap; delete[] tid;    
    delete[] c; delete[] d;
    
    return 0;
}
