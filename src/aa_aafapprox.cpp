/*
 * aa_aafapprox.cpp -- Standart non-affine operations
 * Copyright (c) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
 * Copyright (c) 2004 LIRIS (University Claude Bernard Lyon 1)
 * Copyright (c) 2005 Nathan Hurst
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


#include "aa.h"
#include <algorithm>
#include <cmath>

#include "aa_util.h"


// Operator  *

AAF AAF::operator * (const AAF & P) const {
    unsigned l1 = length;
    unsigned l2 = P.length;

    unsigned * id1 = indexes;
    unsigned * id2 = P.indexes;

    double * va1 = coefficients;
    double * va2 = P.coefficients;

    unsigned * pu1 = id1;
    unsigned * pu2 = id2;

    AAF Temp(cvalue*P.cvalue);  // Create our resulting AAF
    
    Temp.indexes = new unsigned [l1+l2+1];
    unsigned * idtemp=Temp.indexes;


    // Fill the indexes array

    unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
    unsigned ltemp=fin-idtemp;

    Temp.coefficients = new double [ltemp+1];
    double * vatempg=Temp.coefficients;

    Temp.length = ltemp+1;


    // Fill the coefficients array

    for (unsigned i = 0; i < ltemp; i++)
    {
        unsigned a = pu1-id1;
        unsigned b = pu2-id2;

        if (a==l1 || id1[a]!=idtemp[i])
	{
            vatempg[i] = cvalue*va2[b];  // cvalue*va2[b]+(P.cvalue)*0
            pu2++;
            continue;
	}

        if (b==l2 || id2[b]!=idtemp[i])
	{
            vatempg[i] = (P.cvalue)*va1[a];  // cvalue*0+(P.cvalue)*va1[a]
            pu1++;
            continue;
	}

        vatempg[i] = cvalue*va2[b] + (P.cvalue)*va1[a];
        pu1++;
        pu2++;
    }


    // Compute the error
    // in a new noise symbol

    Temp.indexes[ltemp]=inclast();
    Temp.coefficients[ltemp]=rad()*(P.rad());

    Temp.special = binary_special(special, P.special);
    
    return Temp;

}


// Operator  /
// It's a non affine-operation
// We use the identity x/y = x * (1/y)

AAF AAF::operator / (const AAF & P) const {
    return (*this)*inv(P);
}


// Square root operator
// It's a non affine-operation
// We use the Chebyshev approximation

AAF sqrt(const AAF & P) {
    handle_infinity(P);
    // sqrt(x) is approximated by f(x)=alpha*x+dzeta
    // delta is the maximum absolute error

    const double a = P.convert().left(); // [a,b] is our interval
    const double b = P.convert().right();
    AAF_TYPE type;
    if(a >= 0)
        type = AAF_TYPE_AFFINE;
    else if(b < 0) 
        type = AAF_TYPE_NAN;
    else if(a < 0) // undefined, can we do better?
        type = (AAF_TYPE)(AAF_TYPE_AFFINE | AAF_TYPE_NAN);
    //type = (AAF_TYPE)(type | P.special);
    
    const double t = (sqrt(a)+sqrt(b));


    const double alpha = 1/t; // alpha is the slope of the line r(x) that
    // interpolate (a, sqrt(a)) and (b, f(b))

    // dzeta calculation:
    const double dzeta = (t/8)+0.5*(sqrt(a*b))/t;

    // Calculation of the error
    const double rdelta = (sqrt(b)-sqrt(a));
    const double delta = rdelta*rdelta/(8*t);

    return AAF(P, alpha, dzeta, delta, type);
}


// Inverse (1/x) operator
// It's a non-affine operation
// We use mini-range approximation
// because undershoot can be high with Chebyshev here

AAF inv(const AAF & P) {
    handle_infinity(P);
    double a = P.convert().left();
    double b = P.convert().right();
    if(P.is_infinite() || (a <= 0) && (b >= 0)) {
        return AAF(interval(-INFINITY, INFINITY));
    }

    // a := min(abs(a), abs(b))
    // b := max(abs(a), abs(b))

    const double t1 = fabs(a);
    const double t2 = fabs(b);

    a= t1 <? t2;  // min(t1,t2)
    b= t1 >? t2;  // max(t1,t2)

    // Derivative of 1/x is -1/x*x

    const double alpha=-1/(b*b);

    interval i((1/a)-alpha*a, 2/b);//-alpha*b);
    double dzeta = i.mid();

    if ((P.convert().left()) < 0) dzeta = -dzeta;

    return AAF(P, alpha, dzeta, i.radius(), P.special);
}

AAF abs(const AAF & P) {
    if(P.strictly_neg())
        return -P;
    if(P.straddles_zero()) {
        AAF Temp(P);
        Temp.cvalue = fabs(Temp.cvalue)/2;
        Temp.special = P.special;
        
        for (unsigned i=0; i<P.length; i++)
            Temp.coefficients[i]=(Temp.coefficients[i])/2;
        return Temp;
    }
    return P;
}

AAF sqr(const AAF & P) {
    return P*P;
}

// Power function
// only for integer exponents

AAF pow(const AAF & P, int exp) {
    handle_infinity(P);
    if (exp == 0) {
        return 1;
    } else if (exp > 0) {
        if(exp & 1)
            return sqr(pow(P, exp>>1))*P;
        else
            return sqr(pow(P, exp>>1));
    } else {
        return inv(pow(P, -exp));
    }
}

AAF pow(const AAF & P, double xp) {
    return exp(xp*log(P));
}

// Exponential operator
// It's a non affine-operation

AAF exp(const AAF & P) {
    handle_infinity(P); // infinity maps to [0, infty) here
    // exp(x) is approximated by f(x)=alpha*x+dzeta
    // delta is the maximum absolute error

    const double a = P.convert().left(); // [a,b] is our interval
    const double b = P.convert().right();
    
    const double ea = exp(a);
    const double eb = exp(b);
    if((ea == INFINITY) || (eb == INFINITY)) {
        // Printing from a numeric library is generally frowned apon,
        // but what to do instead?  This corresponds to an essential 
        // singularity.
        //printf("essential infinity at %g, %g -> (%g, %g)\n", a, b, ea, eb);
        return AAF(interval(-INFINITY, INFINITY));
    }
    
    const double alpha = (eb-ea)/(b-a);
    // alpha is the slope of the line r(x) that
    // interpolate (a, exp(a)) and (b, exp(b))
    const double xs = log(alpha);// the x of the maximum error
    const double maxdelta = alpha*(xs - 1 - a)+ea;

    // dzeta calculation:
    const double dzeta = alpha*(1 - xs);

    // Calculation of the error
    const double delta = maxdelta/2;

    return AAF(P, alpha, dzeta, delta, P.special);
}

// Logarithm operator
// It's a non affine-operation

AAF log(const AAF & P) {
    handle_infinity(P); // infinity needs to map to NaN here

    const double a = P.convert().left(); // [a,b] is our interval
    const double b = P.convert().right();
    
    AAF_TYPE type;
    if(a > 0)
        type = AAF_TYPE_AFFINE;
    else if(b < 0) { // no point in continuing
        type = AAF_TYPE_NAN;
        return AAF(type);
    }
    else if(a <= 0) { // undefined, can we do better?
        type = (AAF_TYPE)(AAF_TYPE_AFFINE | AAF_TYPE_NAN);
        return AAF(type);
        // perhaps we should make a = 0+eps and try to continue?
    }
    
    const double la = log(a);
    const double lb = log(b);
    
    const double alpha = (lb-la)/(b-a);
    // alpha is the slope of the line r(x) that
    // interpolate (a, exp(a)) and (b, exp(b))
    const double xs = 1/(alpha);// the x of the maximum error
    const double ys = (alpha*(xs - a)+la);
    const double maxdelta = log(xs) - ys;

    // dzeta calculation:
    const double dzeta = alpha*(-xs)+(log(xs)+ys)/2;

    // Calculation of the error
    const double delta = maxdelta/2;

    return AAF(P, alpha, dzeta, delta, type);
}



/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0))
  indent-tabs-mode:nil
  fill-column:80
  End:
*/


// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4 :
