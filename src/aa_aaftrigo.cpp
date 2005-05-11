/*
 * aa_aaftrigo.cpp -- Trigonometric operations (all non-affine)
 * Copyright (c) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
 * Copyright (c) 2004 LIRIS (University Claude Bernard Lyon 1)
 * Copyright (c) 2005 Nathan Hurst
 *
 * This file is part of libaa.
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
 * License along with libaa; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */



#include "aa.h"
#include <algorithm>
#include <cmath>
#include "aa_util.h"

#define PI (4*atan(1.0))
#define NPTS 8 // number of points for the linear
               // regression approximation


// Sine function
// sine isn't montonic and the second derivative change its sign
// inside the defined interval
// so we use the approximation by least squares


AAF sin(const AAF & P)
{
    if(P.is_infinite())
        return AAF(interval(-1,1));
    interval i = P.convert();

    double w = i.width();

    const double a = i.left();
    const double b = i.right();


    // y' = alpha*x+dzeta , the regression line
    // approximate y = sin(x)

    double alpha, dzeta, delta;
    
    if (w >= 2*PI ) {
// the trivial case, the interval is larger than 2*PI
        // y' = 0 , delta = 1 cause -1 <= sin(x) <= +1
        alpha = 0.0;
        dzeta = 0.0;
        delta = 1.0;
    } else {
// case of the least squares
        double x[NPTS];
        double y[NPTS];
        double r[NPTS]; // residues, r[i] = y[i]-y'[i]
        
        x[0] = a;
        y[0] = sin(a);
        x[NPTS-1] = b;
        y[NPTS-1] = sin(b);

        double pas = w/(NPTS-1);

        for (unsigned i=1; i< NPTS-1; i++) {
            x[i] = x[i-1]+pas;
            y[i] = sin(x[i]);
	}


        // Calculation of xm and ym , averages of x and y

        double xm = 0;
        double ym = 0;
        
        for (unsigned i=0; i<NPTS; i++) {
            xm = xm + x[i];
            ym = ym + y[i];
	}

        xm = xm/NPTS;
        ym = ym/NPTS;

        // Calculation of alpha and dzeta

        double temp2 = 0;
        alpha = 0;

        for (unsigned i = 0; i < NPTS; i++) {
            const double temp1 = x[i] - xm;
            alpha += y[i]*temp1;
            temp2 += temp1*temp1;
	}

        alpha = alpha/temp2;  // final alpha
        dzeta = ym - alpha*xm; // final dzeta


        // Calculation of the residues
        // We use the absolute value of the residues!

        for (unsigned i = 0; i < NPTS; i++)	{
            r[i] = fabs(y[i] - (dzeta+alpha*x[i]));
	}


        // The error delta is the maximum
        // of the residues (in absolute values)

        delta = *std::max_element(r, r+NPTS);
    }


    return AAF(P, alpha, dzeta, delta, P.special);
}


// Cosine function
// we use the identity cos(x)=sin(x+PI/2)

AAF cos(const AAF & P) {
    return sin(P+PI/2);
}


// Tangent function
// we use the identity tan(x)=sin(x)/cos(x)
// Due to the nature of the tan fct remember that
// we can have infinite value with small intervals

AAF tan(const AAF & P)
{
    return sin(P)/cos(P);
}


// Cotangent function
// we use the identity cotan(x)=cos(x)/sin(x)
// Due to the nature of the cotan fct remember that
// we can have infinite value with small intervals

AAF cotan(const AAF & P){
    return cos(P)/sin(P);
}


// Hyperbolics

AAF cosh(const AAF & P) {
    return (exp(P) + exp(-P))/2;
}

AAF sinh(const AAF & P) {
    return (exp(P) - exp(-P))/2;
}

AAF tanh(const AAF & P) {
    AAF ep = exp(P);
    AAF pe = exp(-P);
    return (ep-pe)/(ep+pe);
}

AAF acosh(const AAF & P) {
    return log(P + sqrt(sqr(P) - 1));
}

AAF asinh(const AAF & P) {
    if(P.getcenter() < 0)
        return -asinh(-P);
    return(log(P + sqrt(sqr(P) + 1)));
}

AAF atanh(const AAF & P) {

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
