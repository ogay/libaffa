/*
 * aa_interval.cpp -- Implementation of the interval class
 * Copyright (c) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
 * Copyright (c) 2004 LIRIS (University Claude Bernard Lyon 1)
 *
 * This file is part of libaffa.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with libaffa; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include "aa_interval.h"
#include "aa_rounding.h"
#include <cstdio>
#include <iostream>
#include <cmath>

// interval(), modlo(), modhi(), modlohi() are not
// inline functions as they arent used by the lib
// but they are useful for applications

// Create a null interval
// Useful if we only declare a variable

interval:: interval():
    lo(0), hi(0)
{
}


//  Affectation operator

interval & interval::operator = (const interval & I)
{
    if (&I!=this)
    {
        lo = I.lo;
        hi = I.hi;
    }

    return *this;
}


// Modify the lower bound of an interval

void interval::mod_lo(const double low) {
    lo = low;
}


// Modify the higher bound of an interval

void interval::mod_hi(const double high) {
    hi = high;
}


// Modify the two bounds of an interval

void interval::mod_lo_hi(const double low, const double high) {
    lo = low;
    hi = high;
}


// Calculate the midpoint of an interval
// i.e (lo+hi)/2

double interval::mid() const {
    double t0,t1;

    aa_rnd_t mode = aa_fegetround();

    aa_fesetround(AA_DOWNWARD);
    t0 = lo*0.5;

    aa_fesetround(AA_UPWARD);
    t1 = hi*0.5;

    aa_fesetround(mode);

    return t0 + t1;
}


// Calculate the radius of an interval
// i.e (m-lo >=hi-m ? m-lo : hi-m)

double interval::radius() const {
    double m = mid();

    double t0,t1;

    aa_rnd_t mode = aa_fegetround();

    aa_fesetround(AA_DOWNWARD);
    t0 = m-lo;

    aa_fesetround(AA_UPWARD);
    t1 = hi-m;

    aa_fesetround(mode);

    return (t0 >= t1 ? t0 : t1);
}


// Istream input of an interval

std::istream & operator >> (std::istream & s, interval &I) {
    // Accepts the forms
    // x or [x] or [x,y]

    double lo = 0, hi = 0;
    char c = 0;

    s >> c;
    if (c == '[') {
        s >> lo >> c;
        if(c == ',')
            s >> hi >> c;
        //if (c != ']') s.clear(ios_base::badbit);
    } else {
        //s.putback(c);
        s >> lo;
    }

    if(s)
        I = interval(lo, hi);

    return s;

}


// Print an interval to stdout

void interval::int_vprint() const {
    printf("[%f,%f]\n", lo, hi);
}


// Ostream output of an interval

std::ostream & operator << (std::ostream & s, const interval &I) {

    // s.setf(0, ios_base::floatfield);
    // cause we don't want to display in scientific format

    s << "[" << I.left() << "," << I.right() << "]";
    return s;
}


// Caculate the minimal 2PI periodic interval
// of an interval, e.g. : [5*PI, 6*PI] -> [PI, 2PI]
// In fact it is the closest interval to 0

interval min_trigo( const interval &I) {
    // This function is no more needed
    // for our new algorithm of the sine of an AAF

    double a = I.left();
    double b = I.right();
    double t1, t2 ;

    t1 = floor(a/(2*M_PI));
    t2 = b - a;

    a = a - (t1*2*M_PI);
    b = a + t2;

    return interval(a,b);
}

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/


// vim: filetype=c++:expandtab:shiftwidth=4:tabstop=8:softtabstop=4 :
