/*
 * aa_aaftrigo.cpp -- Trigonometric operations (all non-affine)
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



#include "aa.h"
#include <algorithm>
#include <cmath>

#define PI (4*atan(1.0))
#define NPTS 8 // number of points for the linear
               // regression approximation


// Sine function
// sine isn't montonic and the second derivative change its sign
// inside the defined interval
// so we use the approximation by least squares


AAF sin(const AAF & P)
{

  double a, b;
  double alpha, dzeta, delta;
  Interval i = P.convert();

  double w = i.width();

  //i = mintrigo(i); // no more needed

  a = i.getlo();
  b = i.gethi();


  // y' = alpha*x+dzeta , the regression line
  // approximate y = sin(x)

  double x[NPTS];
  double y[NPTS];
  double r[NPTS]; // residues, r[i] = y[i]-y'[i]

  double xm = 0;
  double ym = 0;

  if (w >= 2*PI ) // the trivial case, the interval is larger
                  // than 2*PI
    {

      // y' = 0 , delta = 1 cause -1 <= sin(x) <= +1
      alpha = 0.0;
      dzeta = 0.0;
      delta = 1.0;
    }
  else // case of the least squares
    {

      x[0]=a;
      y[0]=sin(a);
      x[NPTS-1]=b;
      y[NPTS-1]=sin(b);

      double pas=w/(NPTS-1);

      for (unsigned i=1; i< NPTS-1; i++)
	{
	  x[i]=x[i-1]+pas;
	  y[i]=sin(x[i]);
	}


      // Calculation of xm and ym , averages of x and y

      for (unsigned i=0; i<NPTS; i++)
	{
	  xm=xm+x[i];
	  ym=ym+y[i];
	}

      xm=xm/NPTS;
      ym=ym/NPTS;



      // Calculation of alpha and dzeta

      double temp1;
      double temp2=0;
      alpha = 0;

      for (unsigned i=0; i<NPTS; i++)
	{
	  temp1=x[i]-xm;
	  alpha+=y[i]*temp1;
	  temp2+=temp1*temp1;
	}


      alpha=alpha/temp2;  // final alpha
      dzeta=ym-alpha*xm; // final dzeta


      // Calculation of the residues
      // We use the absolute value of the residues!

      for (unsigned i=0; i<NPTS; i++)
	{
	  r[i]=fabs(y[i]-(dzeta+alpha*x[i]));
	}


     // The error delta is the maximum
     // of the residues (in absolute values)

      double *ptr;

      ptr = max_element(r, r+NPTS);

      delta = *ptr;

    }



  // z0 = alpha*x0 + dzeta

  AAF Temp(alpha*(P.cvalue)+dzeta);

  Temp.length=(P.length)+1;
  Temp.coefficients = new double [Temp.length];
  Temp.indexes = new unsigned [Temp.length];

  // zi = alpha*xi

  for (unsigned i=0; i < P.length;i++)
    {
      Temp.indexes[i]=P.indexes[i];
      Temp.coefficients[i]=alpha*(P.coefficients[i]);
    }

  // zk = delta

  Temp.indexes[P.length]=Temp.inclast();   // the error indx
  Temp.coefficients[P.length]=delta;

  return Temp;


}


// Cosine function
// we use the identity cos(x)=sin(x+PI/2)

AAF cos(const AAF & P)
{

  AAF Temp = P;
  return sin(Temp+PI/2);

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

AAF cotan(const AAF & P)
{

  return cos(P)/sin(P);

}
