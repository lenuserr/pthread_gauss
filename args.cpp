#include "inc.h"

Args::Args(Args&& x) {
    a = x.a; x.a = nullptr;
    b = x.b; x.b = nullptr;
    c = x.c; x.c = nullptr;
    n = x.n;
    m = x.m;
    p = x.p;
    k = x.k;
    s = x.s;
    filename = x.filename;
}

Args& Args::operator=(Args&& x) {
    a = x.a; x.a = nullptr;
    b = x.b; x.b = nullptr;
    c = x.c; x.c = nullptr;
    n = x.n;
    m = x.m;
    p = x.p;
    k = x.k;
    s = x.s;
    filename = x.filename;
    return *this;
}

Args::~Args() {
    a = nullptr;
    b = nullptr;
    c = nullptr;        
}  
