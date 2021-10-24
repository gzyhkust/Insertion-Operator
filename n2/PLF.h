#ifndef TDGT_PL_H
#define TDGT_PL_H

#include <vector>
#include <utility>
#include <cmath>
#include <memory>
#include <initializer_list>
#include "Segment.h"
#include "Constants.h"

using namespace std;

/**
 * Piecewise linear function for edge weight, travel time function of paths
 */
struct PLF {
    typedef shared_ptr<vector<Segment>> FPtr;
    typedef vector<Segment>::size_type size_type;

    FPtr f;

    PLF() : f(make_shared<vector<Segment>>()) {}

    PLF(initializer_list<Segment> segList) : f(make_shared<vector<Segment>>(segList)) {}

//copy constructors
    PLF(const PLF &other) {
        if (other.f)
            f = make_shared<vector<Segment>>(*other.f);
    }

    PLF(PLF &&rhs) noexcept : f(std::move(rhs.f)) {}

//copy assignments
    PLF &operator=(const PLF &other) {
        if (other.f)
            f = make_shared<vector<Segment>>(*other.f);
        return *this;
    }

    PLF &operator=(PLF &&rhs) noexcept {
        f = rhs.f;
        return *this;
    }

    Segment &operator[](size_type id) {
        return (*f)[id];
    }

    friend ostream &operator<<(ostream &out, const PLF &plf) {
        for (auto &e:*plf.f)
            out << e << ",  ";
        return out;
    }

    //find segment for time t
    vector<Segment>::iterator dpt2seg(double t) {
        auto beg = (*f).begin(), end = (*f).end();
        auto mid = beg + (end - beg) / 2;
        while (mid != end) {// find the first t* > t, or end
            if (ge(t, mid->t)) // t* in [mid+1,end]
                beg = mid + 1;
            else    //t* in [begin, mid]
                end = mid;
            mid = beg + (end - beg) / 2;
        }
        return mid - 1;
    }

    inline double dpt2wgt(double t, vector<Segment>::iterator s1) {  //weight func
        auto s2 = s1 + 1;
        if (s2 != f->end()) {// adopt s2 to get slop
            //assure s1->t≤t≤s2->t (right equal happen for minimize with same next point)
            assert(lt(s1->t, s2->t) && ge(t, s1->t) && le(t, s2->t));
            double k = (s2->w - s1->w) / (s2->t - s1->t);
            assert(gt(k, -1)); // assure FIFO property
            return k * (t - s1->t) + s1->w;
        } else { //last point
            assert(ge(t, s1->t));//&& le(t, TMAX)
            return s1->w;
        }
    }

    inline double dpt2wgt(double t) {
        auto pseg = dpt2seg(t);
        return dpt2wgt(t, pseg);
    }

    // get dpt time for an arrival time
    inline double arr2dpt(double arr, std::vector<Segment>::iterator s1) {
        auto s2 = s1 + 1;
        if (s2 != f->end()) {
            assert(lt(s1->t, s2->t) && ge(arr, s1->t + s1->w) && le(arr, s2->t + s2->w));  //assure s2 is not redundant
            double k = (s2->w - s1->w) / (s2->t - s1->t);
            assert(gt(k, -1)); // assure FIFO property
            double arr1 = s1->w + s1->t;
            return s1->t + (arr - arr1) / (k + 1);
        } else {
            assert(ge(arr, s1->t + s1->w));
            return arr - s1->w;
        }
    }

    // get arrival time from dpt time
    inline double dpt2arr(double t, std::vector<Segment>::iterator s1) {
        auto s2 = s1 + 1;
        if (s2 != f->end()) {
            assert(lt(s1->t, s2->t) && ge(t, s1->t) && le(t, s2->t));  //assure s2 is not redundant
            double k = (s2->w - s1->w) / (s2->t - s1->t);
            assert(gt(k, -1)); // assure FIFO property
            return k * (t - s1->t) + s1->w + t;
        } else {    //last point
            assert(ge(t, s1->t));//&& le(t, TMAX)
            return s1->w + t;
        }
    }

    inline double dpt2arr(double td) {
        auto pseg = dpt2seg(td);
        return dpt2arr(td, pseg);
    }

    // push a new segment point, initial  value: (ts,INT_MAX), intv=INTV_CNTED (directly connected),
    inline void push(double &pre_k, double dpt, double wgt, int intv) {
        auto &segs = *f;
        assert(le(dpt, TMAX) && le(wgt, INT_MAX) && gt(wgt, 0));
        if (segs.empty())
            segs.emplace_back(dpt, wgt, intv);
        else if (neq(segs.back().t, dpt)) {
            double new_k = (wgt - segs.back().w) / (dpt - segs.back().t);
            assert(gt(dpt, segs.back().t) && gt(new_k, -1)); // assure FIFO property
            // the next point wo pw0 may have different intv
            if (eq(pre_k, new_k)) {
                segs.back().t = dpt;
                segs.back().w = wgt;
                segs.back().intv = intv;
            } else {
                segs.emplace_back(dpt, wgt, intv);
                pre_k = new_k;
            }
        }
    }

    //compute new weight func f(g(t)),t∈[ts,te]  regard g(t) as arr func(only used for time-interval query)
    void compound(double ts, double te, PLF &PLF0, PLF &PLFc, int intv) {
        double pre_k = INT_MAX;
        assert(eq(PLF0.f->front().t, ts));
        auto pw0 = PLF0.f->begin();
        double dpt = pw0->t;
        double wgt0 = pw0->w;
        double arr0 = dpt + wgt0;
        auto pw1 = dpt2seg(arr0);// for edge (v_s, v_e), find dpt@v_s for dpt t at start v
        double wgt1 = dpt2wgt(arr0, pw1);
        PLFc.push(pre_k, dpt, wgt0 + wgt1, intv);
        while (true) {
            bool valid0 = pw0 + 1 < PLF0.f->end();
            bool valid1 = pw1 + 1 < f->end();
            int turn; // next point from base or from cur plf(1)
            if (valid0 && valid1) {
                double next_arr0 = (pw0 + 1)->t + (pw0 + 1)->w;
                if (lt((pw1 + 1)->t, next_arr0)) turn = 1;
                else turn = 0;
            } else if (valid1) turn = 1;
            else if (valid0) turn = 0;
            else break;
            if (turn) {
                pw1++;
                dpt = PLF0.arr2dpt(pw1->t, pw0);
                if (valid0) assert(lt(dpt, (pw0 + 1)->t));
                wgt0 = pw1->t - dpt;
                wgt1 = pw1->w;
            } else {
                pw0++;
                dpt = pw0->t;
                wgt0 = pw0->w;
                arr0 = dpt + wgt0;
                wgt1 = dpt2wgt(arr0, pw1);
            }
            PLFc.push(pre_k, dpt, wgt0 + wgt1, intv);
            if (ge(dpt, te))
                break;
        }
    }

    /**
     * bug: slop>-1, and increase after compound
     * vs to intv, then intv to ve,compound to get travel time func from vs to ve detouring by intv
     * compound for f_intv2ve(f_vs2intv), i.e., PLc = *this->f(PL0)
     * @param PLF0 travel time func from vs to intv
     * @param PLFc compounded travel time function from vs to intv and finally to ve
     * @param intv the detour vertex to compound
     * @param nidx current tree node
     */
    void compound(PLF &PLF0, PLF &PLFc, int intv) {    //suppose wgts belongs to (vs,ve),arr0 starts at v0
        double pre_k = INT_MAX;
        auto pw0 = PLF0.f->begin();
        double dpt = pw0->t; //dpt from start vertex v0
        double wgt0 = pw0->w;
        double arr0 = dpt + wgt0;
        auto pw1 = dpt2seg(arr0);// for edge (v_s, v_e), find dpt@v_s for dpt t at start v
        double wgt1 = dpt2wgt(arr0, pw1);
        PLFc.push(pre_k, dpt, wgt0 + wgt1, intv);
        while (true) {
            bool valid0 = pw0 + 1 < PLF0.f->end();
            bool valid1 = pw1 + 1 < f->end();
            int turn; // next point from base or from cur plf(1)
            if (valid0 && valid1) {
                double next_arr0 = (pw0 + 1)->t + (pw0 + 1)->w;
                if (lt((pw1 + 1)->t, next_arr0)) turn = 1;
                else turn = 0;
            } else if (valid1) turn = 1;
            else if (valid0) turn = 0;
            else break;
            if (turn) {
                pw1++;
                dpt = PLF0.arr2dpt(pw1->t, pw0);
                if (valid0) assert(lt(dpt, (pw0 + 1)->t));
                wgt0 = pw1->t - dpt;
                wgt1 = pw1->w;
            } else {
                pw0++;
                dpt = pw0->t;
                wgt0 = pw0->w;
                arr0 = dpt + wgt0;
                wgt1 = dpt2wgt(arr0, pw1);
            }
            PLFc.push(pre_k, dpt, wgt0 + wgt1, intv);
        }
    }

    //get min from two linear segments, which may cross.
    inline void minseg(double &pre_k, double t1, double t2, double f1, double f2,
                       double g1, double g2, int intvf, int intvg) {
        assert(ge(t2, t1)); //t2≥t1
        if (eq(t1, t2)) return; //not interval
        if (eq(f1, g1)) {// f1=g1
            if (le(f2, g2)) push(pre_k, t1, f1, intvf);// f1=g1, f2≤g2
            else push(pre_k, t1, g1, intvg);// f1=g1, f2>g2
        } else { // f1!=g1
            double delta1 = g1 - f1;
            double delta2 = f2 - g2;
            int flag_cross = 0;
            if (gt(delta1, 0)) { // f1<g1
                push(pre_k, t1, f1, intvf);
                if (gt(delta2, 0)) flag_cross = -1; //f2>g2, crossed
            } else { //f1>g1
                push(pre_k, t1, g1, intvg);
                if (lt(delta2, 0)) flag_cross = 1; //f2<g2, crossed
            }
            if (flag_cross) {
                double denominator = delta1 + delta2;
                double inter_t = (delta1 * t2 + delta2 * t1) / denominator;
                double inter_w = (delta1 * g2 + delta2 * g1) / denominator;
                if (flag_cross == 1) push(pre_k, inter_t, inter_w, intvf);//f2<g2
                else push(pre_k, inter_t, inter_w, intvg);//f2>g2
            }
        }
    }

    //min{f,PLFc->f}
    void minimize(PLF &PLFc) {
        double pre_k = INT_MAX;
        assert(!PLFc.f->empty());
        PLF PL_new;
        auto pw0 = f->begin(), pwc = PLFc.f->begin();
        assert(pw0->t == pwc->t); // same start time
        double t1 = pw0->t, t2, f1 = pw0->w, f2, g1 = pwc->w, g2;
        while (true) {
            int turn; //get next ptr to expand
            bool valid0 = pw0 + 1 < f->end();
            bool valid1 = pwc + 1 < PLFc.f->end();
            if (valid0 && valid1) {
                if (lt((pw0 + 1)->t, (pwc + 1)->t)) turn = 0;
                else turn = 1;
            } else if (valid0) turn = 0;
            else if (valid1) turn = 1;
            else break;
            if (turn) { // pwc turn
                t2 = (pwc + 1)->t;
                f2 = dpt2wgt(t2, pw0);
                g2 = (pwc + 1)->w;
            } else { //pw0 turn
                t2 = (pw0 + 1)->t;
                f2 = (pw0 + 1)->w;
                g2 = PLFc.dpt2wgt(t2, pwc);
            }
            PL_new.minseg(pre_k, t1, t2, f1, f2, g1, g2,
                          pw0->intv, pwc->intv);
            assert(lt(t2, TMAX));
            t1 = t2;
            f1 = f2;
            g1 = g2;
            if (turn) pwc++; else pw0++;
        }
        PL_new.minseg(pre_k, t1, TMAX, f1, pw0->w, g1, pwc->w, // last point ,k = 0
                      pw0->intv, pwc->intv);
        f = PL_new.f;
    }

    bool min4sync(PLF &PLFf, int inter_node_id) {
        bool changed = false;
        double pre_k = INT_MAX;
        assert(!PLFf.f->empty());
        PLF PL_new;
        auto pw0 = f->begin(), pwc = PLFf.f->begin();
        assert(pw0->t == pwc->t); // same start time
        double t1 = pw0->t, t2, f1 = pw0->w, f2, g1 = pwc->w, g2;
        while (true) {
            int turn; //get next ptr to expand
            bool valid0 = pw0 + 1 < f->end();
            bool valid1 = pwc + 1 < PLFf.f->end();
            if (valid0 && valid1) {
                if (lt((pw0 + 1)->t, (pwc + 1)->t)) turn = 0;
                else turn = 1;
            } else if (valid0) turn = 0;
            else if (valid1) turn = 1;
            else break;
            if (turn) { // pwc turn
                t2 = (pwc + 1)->t;
                f2 = dpt2wgt(t2, pw0);
                g2 = (pwc + 1)->w;
            } else { //pw0 turn
                t2 = (pw0 + 1)->t;
                f2 = (pw0 + 1)->w;
                g2 = PLFf.dpt2wgt(t2, pwc);
            }
            if (gt(f1, g1))
                changed = true;
            PL_new.minseg(pre_k, t1, t2, f1, f2, g1, g2,
                          pw0->intv, inter_node_id);
            assert(lt(t2, TMAX));
            t1 = t2;
            f1 = f2;
            g1 = g2;
            if (turn) pwc++; else pw0++;
        }
        if (gt(f1, g1))
            changed = true;
        PL_new.minseg(pre_k, t1, TMAX, f1, pw0->w, g1, pwc->w, // last point ,k = 0
                      pw0->intv, inter_node_id);
        f = PL_new.f;
        return changed;
    }

    std::shared_ptr<std::vector<Segment>> getSlice(double ts, double te) {
        auto ps = dpt2seg(ts), pe = dpt2seg(te);
        auto ans = std::make_shared<std::vector<Segment>>();
        ans->emplace_back(ts, dpt2wgt(ts, ps));
        if (ps == pe)  //only concern one segment
            return ans;
        for (auto p = ps + 1; p <= pe; p++) {//get all in range
            ans->push_back(*p);
        }
        return ans;
    }

    void set_pre_vexs(int marker) {
        for (auto &seg:*f)
            seg.intv = marker;
    }

    // connected in node_id, synced segments.
    void set_pre_nodes(int nega_node_id) {
        for (auto &seg:*f)
            if (seg.intv > 0)
                seg.intv = nega_node_id;
    }
};


#endif //TDGT_PL_H
