#ifndef TDGT_UTILS_H
#define TDGT_UTILS_H

#include <iostream>
#include <stack>
#include <vector>
#include "Constants.h"


#define PRINT_INFO(info) std::cout<<__TIME__ << info<<__FILE__ <<__LINE__ <<std::endl;
#define PRINT_BUILD(info)  std::cout<< info << std::endl;

inline bool le(const double &x, const double &y) { return x <= y + EPSILON; }

inline bool lt(const double &x, const double &y) { return x + EPSILON < y; }

inline bool eq(const double &x, const double &y) { return fabs(x - y) <= EPSILON; }

inline bool neq(const double &x, const double &y) { return !eq(x, y); }

inline bool gt(const double &x, const double &y) { return lt(y, x); }

inline bool ge(const double &x, const double &y) { return le(y, x); }




#endif //TDGT_UTILS_H
