/*
 * aa_aafcommon.cpp -- Common functions used to manipulate AAF
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include "aa.h"
#include <cstdio>
#include <iostream>


unsigned AAF::last = 0; // at beginnnig


// Create an AAF from an array of doubles
// ! For debug purposes

AAF:: AAF(double v0, const double * t1, const unsigned * t2, unsigned T) 
    : special(AAF_TYPE_AFFINE)
{

    length=T;
    cvalue=v0;
    coefficients = new double [length];
    indexes = new unsigned [length];

    for (unsigned i = 0; i < length; i++)
    {
        coefficients[i]=t1[i];
        indexes[i]=t2[i];
    }

    if (indexes[length-1] > last) set_default(indexes[length-1]);

}


// Copy constructor

AAF:: AAF(const AAF &P) 
    : special(P.special)
{
    unsigned plength = P.get_length();
    coefficients = new double [plength];
    indexes = new unsigned [plength];

    cvalue = P.cvalue;
    length = plength;

    for (unsigned i = 0; i<plength; i++)
    {
        coefficients[i] = P.coefficients[i];
        indexes[i] = P.indexes[i];
    }

}


// Create an AAF from an interval

AAF:: AAF(interval iv) {
    unsigned en = inclast();
    coefficients = new double [1];
    indexes = new unsigned [1];
  
    if(iv.width() == INFINITY) {
        cvalue = 0;
        length = 1;
        coefficients[0] = INFINITY;
        indexes[0] = en;
        special = AAF_TYPE_INFINITE;
    } else {
        cvalue=(iv.right()+iv.left())/2;
        length = 1;
    
        coefficients[0]=(iv.right()-iv.left())/2;
        indexes[0]=en;
        special = AAF_TYPE_AFFINE;
    }
}


// AAF destructor

AAF::~AAF()
{
    if (length)
    {
        delete [] coefficients;
        delete [] indexes;
    }
}


//  Affectation operator

AAF & AAF::operator = (const AAF & P)
{
    special = P.special;
    unsigned plength = P.get_length();
    
    if (&P!=this)
    {
        if (length != plength)
        {

            if (length)
            {
                delete [] coefficients;
                delete [] indexes;
            }

            coefficients = new double [plength];
            indexes = new unsigned [plength];
        }

        cvalue = P.cvalue;
        length=plength;
        for (unsigned i = 0; i<plength; i++)
        {
            coefficients[i]=P.coefficients[i];
            indexes[i]=P.indexes[i];
        }
    }

    return *this;
}


// Ostream output of an AAF

std::ostream & operator << (std::ostream & s, const AAF &P)
{

    // s.setf(0, ios_base::floatfield);
    s << "-------------\n";
    s << "Length = " << P.length << "\n";
    s << "v0 = " << P.cvalue << "\n";

    for (unsigned i=0; i < P.length ; i++)
        s << "e" << P.indexes[i] << " = " << P.coefficients[i] << "\n";

    return s;

}


// Print length and coefficients of an AAF to stdout

void AAF::aafprint() const
{
    std::cout << "-------------\n";

    printf("Size = %d\n", length);
    printf("v0 = %f\n", cvalue);

    for (unsigned i=0; i < length ; i++)
        printf("e%d = %f\n", indexes[i], coefficients[i]);
}


// Get the central value x0 of an AAF
// not inline as it isn't used by the lib
// but is useful for applications


double AAF::get_center() const
{
    return cvalue;
}


// Convert an AAF to an interval representation

interval AAF::convert() const
{

    // lower bound == central value of the AAF - the total deviation
    // upper bound == central value of the AAF + the total deviation
    if(is_indeterminate())
        return interval(-INFINITY, INFINITY);
    return interval(cvalue-rad(), cvalue+rad());
}


// Get the total deviation of an AAF
// i.e. the sum of all noise symbols (their abs value)

double AAF::rad() const
{
    double sum=0;

    for (unsigned i=0; i< length; i++) {
        if (coefficients[i] >= 0.0)
            sum+=coefficients[i];
        else
            sum+=-coefficients[i];
    }


    return sum;

}

AAF half_plane(const AAF & P) {
    const double a = P.convert().left(); // [a,b] is our interval
    const double b = P.convert().right();
    
    AAF_TYPE type;
    if(P.special == AAF_TYPE_NAN)
        return P;
    
    if(a > 0)
        return P;
    else if(b < 0) {
        type = AAF_TYPE_NAN;
        return AAF(type);
    }
    else if(b == 0) {
        AAF result(0.);
        result.special = (AAF_TYPE)(AAF_TYPE_AFFINE | AAF_TYPE_NAN);
        return result;
    }
    else if(a <= 0) {
        AAF result;
        result.special = (AAF_TYPE)(AAF_TYPE_AFFINE | AAF_TYPE_NAN);
        unsigned plength = P.get_length();
        result.coefficients = new double [plength];
        result.indexes = new unsigned [plength];
        
        result.cvalue = b/2;;
        result.length = plength;
        
        double rescale = (b - a)/b;
        for (unsigned i = 0; i<plength; i++)
        {
            result.coefficients[i] = (P.coefficients[i]*b);
            result.indexes[i] = P.indexes[i];
        }
        return result;
    }
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
