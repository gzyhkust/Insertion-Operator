#ifndef TDGT_POINT_H
#define TDGT_POINT_H

#include <climits>
#include <iostream>
#include "misc.h"

/**
 * left point of a linear piece,
 * slop determined by following point,
 * or zero for last point.
 */
struct Segment {
    double t;   //departure time
    double w;   //weight, i.e., travel time
    int intv;   // intermediate vertex of current segment or connenction info of intv belong to(negative)
    Segment() : t(0), w(INT_MAX), intv(DE_INTV) {}

    Segment(double _t, double _w) noexcept : t(_t), w(_w), intv(DE_INTV) {}

    Segment(double _t, double _w, int _intv) noexcept : t(_t), w(_w), intv(_intv) {}

    friend std::ostream &operator<<(std::ostream &os, const Segment &seg) {
        os << "(" << seg.t << ", " << seg.w << ", " << seg.intv << ")";
        return os;
    }

    inline bool operator==(const Segment &rhs) {
        return eq(t, rhs.t) && eq(w, rhs.w) && intv == rhs.intv;
    }

    inline bool operator!=(const Segment &rhs) {
        return !(*this == rhs);
    }

};

#endif //TDGT_POINT_H
