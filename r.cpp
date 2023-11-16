#include <cmath>
#include "inc.h"

double r1_eval(int n, double* matrix, double* x, double* b, double* c, double* d) {
    
    matrix_product(n, n, 1, matrix, x, c);
    subtract_vector(n, c, b, d);

    return vector_norm(n, d) / vector_norm(n, b);
}

double r2_eval(int n, double* x) {
    
    double r2 = 0;
    for (int i = 1; i <= n; ++i) {
        r2 += std::fabs(x[i - 1] - (i % 2));
    }

    return r2;
}
