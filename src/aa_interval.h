/*
 * aa_interval.h -- A simple interval class
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



#ifndef AA_INTERVAL_H
#define AA_INTERVAL_H

#include <iostream>


// A class for Interval representation
// the class is used by our AAF class

class Interval
{

 private:

  double lo, hi;

 public:

  Interval();
  Interval(double l, double h);
  Interval & operator = (const Interval & I);

  friend istream & operator >> (istream & s, Interval &I);

  double getlo() const;
  double gethi() const;
  void modlo(const double low);
  void modhi(const double high);
  void modlohi(const double low, const double high);
  double mid() const;
  double radius() const;
  double width() const;
  void intvprint() const;

};

ostream & operator << (ostream & s, const Interval &I);
Interval mintrigo( const Interval &I);


// Interval inline functions

// Create an Interval object

inline Interval:: Interval(double l, double h):
     lo(l), hi(h)
{
}


// Get the lower bound of an Interval

inline double Interval::getlo() const
{
  return lo;
}


// Get the higher bound of an Interval

inline double Interval::gethi() const
{
  return hi;
}

// Calculate the width of an Interval

inline double Interval::width() const
{
  return (hi-lo);
}


#endif  // AA_INTERVAL_H
