/*
 * aa_aaf.h -- Affine Arithmetic class
 * Copyright (C) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
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


#ifndef AA_AAF_H
#define AA_AAF_H

#include "aa_interval.h"
#include <iostream>

// Affine Arithmetic Form

class AAF
{

 private:

  static unsigned last;  // highest noise symbol in use

  double cvalue;       // central value vo
  unsigned length;         // lenght of indexes

  // At creation we don't store null coefficients

  double * coefficients; // values of noise sym
  unsigned * indexes;   // indexes of noise sym

 public:

  AAF(double v0 = 0);
  AAF(double v0, const double * t1, const unsigned * t2, unsigned T);
  AAF(const AAF & P);
  AAF(Interval iv);
  ~AAF();

  AAF & operator = (const AAF & P);
  AAF operator + (const AAF & P);
  AAF operator - (const AAF & P);
  AAF operator * (const AAF & P);
  AAF operator / (const AAF & P);

  friend AAF sqrt(const AAF & P);
  friend AAF inv(const AAF & P);
  friend AAF sin(const AAF & P);

  AAF operator - ();
  AAF operator * (double);

  friend ostream & operator << (ostream & s, const AAF &P);

  void aafprint() const;
  static void set_default(const unsigned val=0);
  static unsigned inclast();
  int getlength() const;
  double getcenter() const;
  Interval convert() const;
  double rad() const;

};

AAF operator * (double, const AAF);
AAF operator + (double, const AAF);
AAF operator - (double, const AAF);
AAF pow(const AAF & P, int exp);
AAF cos(const AAF & P);
AAF tan(const AAF & P);
AAF cotan(const AAF);


// AAF inline functions

// Create a constant AAF of v0

inline AAF:: AAF(double v0):
     cvalue(v0), length(0),
     coefficients(NULL), indexes(NULL)
{
}


// Says highest symbol in use is val

inline void AAF:: set_default(const unsigned val)
{
  last=val;
}


// Increment the highest symbol in use
// i.e. create a new noise symbol

inline unsigned AAF:: inclast()
{
  return ++last;
}


// Get the length of an AAF
// i.e the number of non-null noise symbols

inline int AAF::getlength() const
{
  return length;
}


#endif  // AA_AAF_H
