/*
 * aa_interval.cpp -- Implementation of the interval class
 * Copyright (c) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
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


#include "aa_interval.h"
#include "aa_rounding.h"
#include <cstdio>
#include <iostream>
#include <cmath>

#define PI (4*atan(1.0))


// Interval(), modlo(), modhi(), modlohi() are not
// inline functions as they arent used by the lib
// but they are useful for applications

// Create a null interval
// Useful if we only declare a variable

Interval:: Interval():
  lo(0), hi(0)
{
}


//  Affectation operator

Interval & Interval::operator = (const Interval & I)
{
  if (&I!=this)
    {
      lo = I.lo;
      hi = I.hi;
    }

  return *this;
}


// Modify the lower bound of an Interval

void Interval::modlo(const double low)
{
  lo=low;
}


// Modify the higher bound of an Interval

void Interval::modhi(const double high)
{
  hi=high;
}


// Modify the two bounds of an Interval

void Interval::modlohi(const double low, const double high)
{
  lo=low;
  hi=high;
}


// Calculate the midpoint of an Interval
// i.e (lo+hi)/2

double Interval::mid() const
{

  double t0,t1;

  aa_rnd_t mode = aa_fegetround();

  aa_fesetround(AA_DOWNWARD);
  t0 = lo*0.5;

  aa_fesetround(AA_UPWARD);
  t1 = hi*0.5;

  aa_fesetround(mode);

  return t0+t1;
}


// Calculate the radius of an Interval
// i.e (m-lo >=hi-m ? m-lo : hi-m)

double Interval::radius() const
{
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


// Istream input of an Interval

istream & operator >> (istream & s, Interval &I)
{

  // Accepts the forms
  // x or [x] or [x,y]

  double lo = 0, hi = 0;
  char c = 0;

  s >> c;
  if (c == '[')
    {
      s >> lo >> c;
      if (c == ',') s >> hi >> c;
      //if (c != ']') s.clear(ios_base::badbit);
    }
  else
  {
    //s.putback(c);
    s >> lo;
  }

  if (s) I = Interval(lo, hi);

  return s;

}


// Print an Interval to stdout

void Interval::intvprint() const
{

  printf("[%f,%f]\n", lo, hi);

}


// Ostream output of an Interval

ostream & operator << (ostream & s, const Interval &I)
{

  // s.setf(0, ios_base::floatfield);
  // cause we don't want to display in scientific format

  s << "[" << I.getlo() << "," << I.gethi() << "]\n";
  return s;

}


// Caculate the minimal 2PI periodic interval
// of an interval, e.g. : [5*PI, 6*PI] -> [PI, 2PI]
// In fact it is the closer interval to 0

Interval mintrigo( const Interval &I)
{

  // This function is no more needed
  // for our new algorithm of the sine of an AAF

  double a = I.getlo();
  double b = I.gethi();
  double t1, t2 ;

  t1 = floor(a/(2*PI));
  t2 = b-a;

  a=a-(t1*2*PI);
  b=a+t2;

  return Interval(a,b);
}
