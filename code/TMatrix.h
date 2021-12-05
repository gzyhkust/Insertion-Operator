#ifndef TDGT_TMATRIX_H
#define TDGT_TMATRIX_H

#include <memory>
#include <vector>
//#include "../Entity/PLF.h"
#include "PLF.h"

struct TMatrix {
    unsigned long n;  //TMatrix size is n*n
    PLF *a;//TMatrix
    TMatrix() : n(0), a(nullptr) {}

    ~TMatrix() { delete[] a; }

    //diagonal is not used
    void init(unsigned long _n) {
        this->n = _n;
        a = new PLF[_n * _n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i != j)
                    a[i * n + j] = {Segment()};

    }

    PLF *operator[](int x) { return a + x * n; }
};

#endif //TDGT_TMATRIX_H
