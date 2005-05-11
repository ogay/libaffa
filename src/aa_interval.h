/*
 * aa_interval.h -- A simple interval class
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



#ifndef AA_INTERVAL_H
#define AA_INTERVAL_H

#include <iostream>


// A class for interval representation
// the class is used by our AAF class

class interval
{

 private:

  double lo, hi;

 public:

  interval();
  interval(double l, double h);
  interval & operator = (const interval & I);

  friend std::istream & operator >> (std::istream & s, interval &I);

  double left() const;
  double right() const;
  void modlo(const double low);
  void modhi(const double high);
  void modlohi(const double low, const double high);
  double mid() const;
  double radius() const;
  double width() const;
  void intvprint() const;

  bool straddles_zero() const;
};

std::ostream & operator << (std::ostream & s, const interval &I);
interval mintrigo( const interval &I);


// interval inline functions

// Create an interval object

inline interval:: interval(double l, double h):
     lo(l), hi(h)
{
}


// Get the lower bound of an interval

inline double interval::left() const
{
  return lo;
}


// Get the higher bound of an interval

inline double interval::right() const
{
  return hi;
}

// Calculate the width of an interval

inline double interval::width() const
{
  return (hi-lo);
}

// Returns true iff the interval contains 0

inline bool interval::straddles_zero() const {
  return (lo <= 0) && (hi >= 0);
}

#endif  // AA_INTERVAL_H
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
