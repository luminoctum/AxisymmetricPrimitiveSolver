#ifndef BOUNDARY
#define BOUNDARY
#include "Pooma/Arrays.h"
#include "configure.h"

/**
 * The best way to force boundary condition
 * is to use binary system. 1,2,4,8 represent the four boundaries
 * and their combinations
 */
struct Periodic{
    template<class A> void fix(A& a, const Interval<1>& ci, int edge = 0){
        Interval<1> cx = a.domain();
        Interval<1> lower(ci.first(), cx.last() - ci.length());
        Interval<1> upper(ci.length() + cx.first(), ci.last());
        Interval<1> lowerfix(cx.first(), ci.first() - 1);
        Interval<1> upperfix(ci.last() + 1, cx.last());
        switch (edge) {
            case 0:
                a(lowerfix) = a(upper);
                a(upperfix) = a(lower);
            case 1: 
                a(lowerfix) = a(upper);
                break;
            case 2:
                a(upperfix) = a(lower);
                break;
        }
    }
    template<class A> void fix(A a, const Interval<2>& cij, int edge = 0){
        Interval<2> cxy = a.domain();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        switch (edge) {
            case 0:
                for (int i = xfirst; i < ifirst; i++) a(i, cij[1]) = a(i + ilast - ifirst + 1, cij[1]);
                for (int i = xlast; i > ilast; i--) a(i, cij[1]) = a(i - ilast + ifirst - 1, cij[1]);
                for (int j = yfirst; j < jfirst; j++) a(cij[0], j) = a(cij[0], j + jlast - jfirst + 1);
                for (int j = ylast; j > jlast; j--) a(cij[0], j) = a(cij[0], j - jlast + jfirst - 1);
            case 1:
                for (int i = xfirst; i < ifirst; i++) a(i, cij[1]) = a(i + ilast - ifirst + 1, cij[1]);
                break;
            case 2:
                for (int j = yfirst; j < jfirst; j++) a(cij[0], j) = a(cij[0], j + jlast - jfirst + 1);
                break;
            case 3:
                for (int i = xlast; i > ilast; i--) a(i, cij[1]) = a(i - ilast + ifirst - 1, cij[1]);
                break;
            case 4:
                for (int j = ylast; j > jlast; j--) a(cij[0], j) = a(cij[0], j - jlast + jfirst - 1);
        }
    }
    friend std::ostream& operator<< (std::ostream &os, const Periodic &other){
        os << "Periodic"; return os;
    }
};

struct LinearExtrap{
    template<class A> void fix(A& a, const Interval<1>& ci, int edge = 0){
        Interval<1> cx = a.domain();
        int ifirst = ci.first(), ilast = ci.last();
        int xfirst = cx.first(), xlast = cx.last();
        switch (edge) {
            case 0:
                for (int i = ifirst - 1; i >= xfirst; i--) a(i) = a(i + 1) - a(ifirst + 1) + a(ifirst);
                for (int i = ilast + 1; i <= xlast; i++) a(i) = a(i - 1) + a(ilast) - a(ilast - 1);
                break;
            case 1:
                for (int i = ifirst - 1; i >= xfirst; i--) a(i) = a(i + 1) - a(ifirst + 1) + a(ifirst);
                break;
            case 2:
                for (int i = ilast + 1; i <= xlast; i++) a(i) = a(i - 1) + a(ilast) - a(ilast - 1);
        }
    }
    template<class A> void fix(A& a, const Interval<2>& cij, int edge = 0){
        Interval<2> cxy = a.domain();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        switch (edge) {
            case 0:
                for (int i = ifirst - 1; i >= xfirst; i--) 
                    a(i, cij[1]) = a(i + 1, cij[1]) - a(ifirst + 1, cij[1]) + a(ifirst, cij[1]);
                for (int i = ilast + 1; i <= xlast; i++)
                    a(i, cij[1]) = a(i - 1, cij[1]) + a(ilast, cij[1]) - a(ilast - 1, cij[1]);
                for (int j = jfirst - 1; j >= yfirst; j--) 
                    a(cij[0], j) = a(cij[0], j + 1) - a(cij[0], jfirst + 1) + a(cij[0], jfirst);
                for (int j = jlast + 1; j <= ylast; j++)
                    a(cij[0], j) = a(cij[0], j - 1) + a(cij[0], jlast) - a(cij[0], jlast - 1);
                break;
            case 1:
                for (int i = ifirst - 1; i >= xfirst; i--) 
                    a(i, cij[1]) = a(i + 1, cij[1]) - a(ifirst + 1, cij[1]) + a(ifirst, cij[1]);
                break;
            case 2:
                for (int j = jfirst - 1; j >= yfirst; j--) 
                    a(cij[0], j) = a(cij[0], j + 1) - a(cij[0], jfirst + 1) + a(cij[0], jfirst);
                break;
            case 3:
                for (int i = ilast + 1; i <= xlast; i++)
                    a(i, cij[1]) = a(i - 1, cij[1]) + a(ilast, cij[1]) - a(ilast - 1, cij[1]);
                break;
            case 4:
                for (int j = jlast + 1; j <= ylast; j++)
                    a(cij[0], j) = a(cij[0], j - 1) + a(cij[0], jlast) - a(cij[0], jlast - 1);
        }
    }
    friend std::ostream& operator<< (std::ostream &os, const LinearExtrap &other){
        os << "LinearExtrap"; return os;
    }
};

struct ConstExtrap{
    template<class A> void fix(A& a, const Interval<1>& ci, int edge = 0){
        Interval<1> cx = a.domain();
        int ifirst = ci.first(), ilast = ci.last();
        int xfirst = cx.first(), xlast = cx.last();
        switch (edge) {
            case 0:
                for (int i = ifirst - 1; i >= xfirst; i--) a(i) = a(ifirst);
                for (int i = ilast + 1; i <= xlast; i++) a(i) = a(ilast);
                break;
            case 1:
                for (int i = ifirst - 1; i >= xfirst; i--) a(i) = a(ifirst);
                break;
            case 2:
                for (int i = ilast + 1; i <= xlast; i++) a(i) = a(ilast);
        }
    }
    template<class A> void fix(A& a, const Interval<2>& cij, int edge = 0){
        Interval<2> cxy = a.domain();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        switch (edge) {
            case 0:
                for (int i = xfirst; i < ifirst; i++) a(i, cij[1]) = a(ifirst, cij[1]);
                for (int i = xlast; i > ilast; i--) a(i, cij[1]) = a(ilast, cij[1]);
                for (int j = yfirst; j < jfirst; j++) a(cij[0], j) = a(cij[0], jfirst);
                for (int j = ylast; j > jlast; j--) a(cij[0], j) = a(cij[0], jlast);
                break;
            case 1:
                for (int i = xfirst; i < ifirst; i++) a(i, cij[1]) = a(ifirst, cij[1]);
                break;
            case 2:
                for (int j = yfirst; j < jfirst; j++) a(cij[0], j) = a(cij[0], jfirst);
                break;
            case 3:
                for (int i = xlast; i > ilast; i--) a(i, cij[1]) = a(ilast, cij[1]);
                break;
            case 4:
                for (int j = ylast; j > jlast; j--) a(cij[0], j) = a(cij[0], jlast);
        }
    }
    friend std::ostream& operator<< (std::ostream &os, const ConstExtrap &other){
        os << "ConstExtrap"; return os;
    }
};

template<class T = double>
struct Fixed{
    typedef T Element_t;
    bool unset;
    Element_t lower, upper;
    Array<1, Element_t> left, right, bottom, top;
    Fixed(){ unset = true; }
    template<class A> void fix(A& a, const Interval<1>& ci, int edge = 0) {
        Interval<1> cx = a.domain();
        int ifirst = ci.first(), ilast = ci.last();
        int xfirst = cx.first(), xlast = cx.last();
        if (unset){ 
            lower = a(ifirst); 
            upper = a(ilast); 
            unset = false;
        }
        switch (edge) {
            case 0:
                for (int i = ifirst; i >= xfirst; i--) a(i) = lower;
                for (int i = ilast; i <= xlast; i++) a(i) = upper;
                break;
            case 1:
                for (int i = ifirst; i >= xfirst; i--) a(i) = lower;
                break;
            case 2:
                for (int i = ilast; i <= xlast; i++) a(i) = upper;
        }
    }
    template<class A> void fix(A& a, const Interval<2>& cij, int edge = 0){
        Interval<2> cxy = a.domain();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        Interval<1> ci = cij[0], cj = cij[1];
        if (unset) { 
            left.initialize(cj); right.initialize(cj);
            left = a(ifirst, cj); right = a(ilast, cj);
            bottom.initialize(ci); top.initialize(ci);
            bottom = a(ci, jfirst); top = a(ci, jlast);
            unset = false;
        }
        switch (edge) {
            case 0:
                for (int i = xfirst; i <= ifirst; i++) a(i, cj) =  left;
                for (int i = xlast; i >= ilast; i--) a(i, cj) =  right;
                for (int j = yfirst; j <= jfirst; j++) a(ci, j) = bottom;
                for (int j = ylast; j >= jlast; j--) a(ci, j) = top;
                break;
            case 1:
                for (int i = xfirst; i <= ifirst; i++) a(i, cj) =  left;
                break;
            case 2:
                for (int j = yfirst; j <= jfirst; j++) a(ci, j) = bottom;
                break;
            case 3:
                for (int i = xlast; i >= ilast; i--) a(i, cj) =  right;
                break;
            case 4:
                for (int j = ylast; j >= jlast; j--) a(ci, j) = top;
        }
    }
    friend std::ostream& operator<< (std::ostream &os, const Fixed<T> &other){
        os << "Fixed"; return os;
    }
};

template<class T = double, int width = _AbsorbWidth_, int strength = _AbsorbTimeScale_>
struct Absorbing{
    typedef T Element_t;
    bool unset; 
    Element_t lower, upper;
    Array<1, Element_t> left, right, bottom, top;
    Absorbing() : unset(true) {}
    template<class A> void fix(A& a, const Interval<1>& ci, int edge = 0) {
        Interval<1> cx = a.domain();
        int ifirst = ci.first(), ilast = ci.last();
        int xfirst = cx.first(), xlast = cx.last();
        if (unset){ 
            lower = a(ifirst); 
            upper = a(ilast); 
            unset = false;
        }
        switch (edge) {
            case 0:
                for (int i = ifirst; i >= xfirst; i--) a(i) = lower;
                for (int i = ilast; i <= xlast; i++) a(i) = upper;
                break;
            case 1:
                for (int i = ifirst; i >= xfirst; i--) a(i) = lower;
                break;
            case 2:
                for (int i = ilast; i <= xlast; i++) a(i) = upper;
        }
    }
    template<class A> void fix(A& a, const Interval<2>& cij, int edge = 0){
        Interval<2> cxy = a.domain();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        Interval<1> ci = cij[0], cj = cij[1];
        if (unset){ 
            left.initialize(cj); right.initialize(cj);
            left = a(ifirst, cj); right = a(ilast, cj);
            for (int i = xfirst; i < ifirst; i++) a(i, cj) =  left;
            for (int i = xlast; i > ilast; i--) a(i, cj) =  right;
            bottom.initialize(ci); top.initialize(ci);
            bottom = a(ci, jfirst); top = a(ci, jlast);
            for (int j = yfirst; j < jfirst; j++) a(ci, j) = bottom;
            for (int j = ylast; j > jlast; j--) a(ci, j) = top;
            unset = false;
        } 
        switch (edge) {
            case 0:
                for (int i = xfirst; i < ifirst; i++) a(i, cj) =  left;
                for (int i = ifirst; i < ifirst + width; i++) 
                    a(i, cj) -=  1. / strength * exp(ifirst - i) * (a(i, cj) - left);
                for (int i = xlast; i > ilast; i--) a(i, cj) =  right;
                for (int i = ilast; i > ilast - width; i--)
                    a(i, cj) -=  1. / strength * exp(i - ilast) * (a(i, cj) - right);
                for (int j = yfirst; j < jfirst; j++) a(ci, j) = bottom;
                for (int j = jfirst; j < jfirst + width; j++) 
                    a(ci, j) -=  1. / strength * exp(jfirst - j) * (a(ci, j) - bottom);
                for (int j = ylast; j > jlast; j--) a(ci, j) = top;
                for (int j = jlast; j > jlast - width; j--) 
                    a(ci, j) -=  1. / strength * exp(j - jlast) * (a(ci, j) - top);
                break;
            case 1:
                for (int i = xfirst; i < ifirst; i++) a(i, cj) =  left;
                for (int i = ifirst; i < ifirst + width; i++) 
                    a(i, cj) -=  1. / strength * exp(ifirst - i) * (a(i, cj) - left);
                break;
            case 2:
                for (int j = yfirst; j < jfirst; j++) a(ci, j) = bottom;
                for (int j = jfirst; j < jfirst + width; j++) 
                    a(ci, j) -=  1. / strength * exp(jfirst - j) * (a(ci, j) - bottom);
                break;
            case 3:
                for (int i = xlast; i > ilast; i--) a(i, cj) =  right;
                for (int i = ilast; i > ilast - width; i--)
                    a(i, cj) -=  1. / strength * exp(i - ilast) * (a(i, cj) - right);
                break;
            case 4:
                for (int j = ylast; j > jlast; j--) a(ci, j) = top;
                for (int j = jlast; j > jlast - width; j--) 
                    a(ci, j) -=  1. / strength * exp(j - jlast) * (a(ci, j) - top);
        }
    }
    friend std::ostream& operator<< (std::ostream &os, const Absorbing<T, width, strength> &other){
        os << "Absorbing"; return os;
    }
};

struct Reflective{
    template<class A> void fix(A& a, const Interval<2>& cij, int edge = 0){
        Interval<2> cxy = a.domain();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        Interval<1> ci = cij[0], cj = cij[1];
        switch (edge) {
            case 0:
                for (int i = xfirst; i < ifirst; i++) a(i, cj) = - a(2 * ifirst - i - 1, cj);
                for (int i = xlast; i > ilast; i--) a(i, cj) = - a(2 * ilast - i + 1, cj);
                for (int j = yfirst; j < jfirst; j++) a(ci, j) = - a(ci, 2 * jfirst - j - 1);
                for (int j = ylast; j > jlast; j--) a(ci, j) = - a(ci, 2 * jlast - j + 1);
                break;
            case 1:
                for (int i = xfirst; i < ifirst; i++) a(i, cj) = - a(2 * ifirst - i - 1, cj);
                break;
            case 2:
                for (int j = yfirst; j < jfirst; j++) a(ci, j) = - a(ci, 2 * jfirst - j - 1);
                break;
            case 3:
                for (int i = xlast; i > ilast; i--) a(i, cj) = - a(2 * ilast - i + 1, cj);
                break;
            case 4:
                for (int j = ylast; j > jlast; j--) a(ci, j) = - a(ci, 2 * jlast - j + 1);
        }
    }
    friend std::ostream& operator<< (std::ostream &os, const Reflective &other){
        os << "Reflective"; return os;
    }
};

struct Mirror{
    template<class A> void fix(A& a, const Interval<2>& cij, int edge = 0){
        Interval<2> cxy = a.domain();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        Interval<1> ci = cij[0], cj = cij[1];
        switch (edge) {
            case 0:
                for (int i = xfirst; i < ifirst; i++) a(i, cj) = a(2 * ifirst - i - 1, cj);
                for (int i = xlast; i > ilast; i--) a(i, cj) = a(2 * ilast - i + 1, cj);
                for (int j = yfirst; j < jfirst; j++) a(ci, j) = a(ci, 2 * jfirst - j - 1);
                for (int j = ylast; j > jlast; j--) a(ci, j) = a(ci, 2 * jlast - j + 1);
                break;
            case 1:
                for (int i = xfirst; i < ifirst; i++) a(i, cj) = a(2 * ifirst - i - 1, cj);
                break;
            case 2:
                for (int j = yfirst; j < jfirst; j++) a(ci, j) = a(ci, 2 * jfirst - j - 1);
                break;
            case 3:
                for (int i = xlast; i > ilast; i--) a(i, cj) = a(2 * ilast - i + 1, cj);
                break;
            case 4:
                for (int j = ylast; j > jlast; j--) a(ci, j) = a(ci, 2 * jlast - j + 1);
        }
    }
    friend std::ostream& operator<< (std::ostream &os, const Mirror &other){
        os << "Mirror"; return os;
    }
};

struct Dependent{
    template<class A> void fix(A& a, const Interval<1>& ci, int edge = 0){}
    template<class A> void fix(A& a, const Interval<2>& cij, int edge = 0){}
    friend std::ostream& operator<< (std::ostream &os, const Dependent &other){
        os << "Dependent"; return os;
    }
};

template<class L = Dependent, class B = Dependent, class R = L, class T = B>
class Boundary{
public:
    L left; B bottom; R right; T top;
    template<class A> void fix(A& a, const Interval<2>& cij){
        left.fix(a, cij, 1);
        bottom.fix(a, cij, 2);
        right.fix(a, cij, 3);
        top.fix(a, cij, 4);
    }
    friend std::ostream& operator<< (std::ostream &os, const Boundary<L,B,R,T> &other){
        os  << "-Left- " << other.left << ", "
            << "-Bottom- " << other.bottom << ", "
            << "-Right- " << other.right << ", "
            << "-Top- " << other.top; 
        return os;
    }
};
#endif
