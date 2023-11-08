#include <iostream>
#include "inc.h"

int main(int argc, char* argv[]) {
    // ./a.out n m p r s filename
    int n = std::stoi(argv[1]);
    int m = std::stoi(argv[2]);
    int p = std::stoi(argv[3]);
    int r = std::stoi(argv[4]);
    int s = std::stoi(argv[5]);
    int k = 0; // номер потока
    std::string filename;
    if (!s) {
        filename = argv[argc - 1];
    }   

    double* a = new double[n*n]; 
    double* b = new double[n];
    double* c = new double[n];

    Args* ap = nullptr;
    pthread_t* tid = nullptr;

    ap = new Args[p];
    tid = new pthread_t[p];
    
    for (k = 0; k < p; ++k) {
        ap[k].a = a; ap[k].b = b; ap[k].c = c;
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

    if (ap[0].res < 0) {
        std::cout << "Problems reading from a file " << filename << "\n";
        delete[] a; delete[] b; delete[] c;
        delete[] ap; delete[] tid;   
        return -1;  
    }        

    if (ap[0].method_not_applicable) {
        std::cout << "Method is not applicable" << "\n";
        delete[] a; delete[] b; delete[] c;
        delete[] ap; delete[] tid;   
        return -2;
    }

    delete[] a; delete[] b; delete[] c;
    delete[] ap; delete[] tid;    

    return 0;
}
