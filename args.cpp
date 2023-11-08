#include "inc.h"

Args::Args(Args&& v) {
    a = v.a; v.a = nullptr;
    b = v.b; v.b = nullptr;
    x = v.x; v.x = nullptr;
    n = v.n;
    m = v.m;
    p = v.p;
    k = v.k;
    s = v.s;
    filename = v.filename;
}

Args& Args::operator=(Args&& v) {
    a = v.a; v.a = nullptr;
    b = v.b; v.b = nullptr;
    x = v.x; v.x = nullptr;
    n = v.n;
    m = v.m;
    p = v.p;
    k = v.k;
    s = v.s;
    filename = v.filename;
    return *this;
}

Args::~Args() {
    a = nullptr;
    b = nullptr;
    x = nullptr;        
}  
