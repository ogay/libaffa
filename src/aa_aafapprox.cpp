/*
 * aa_aafapprox.cpp -- Standart non-affine operations
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


// Operator  *

AAF AAF::operator * (const AAF & P)
{


  unsigned l1 = length;
  unsigned l2 = P.length;

  unsigned * id1=indexes;
  unsigned * id2=P.indexes;

  double * va1=coefficients;
  double * va2=P.coefficients;

  unsigned * pu1=id1;
  unsigned * pu2=id2;

  AAF Temp(cvalue*P.cvalue);  // Create our resulting AAF

  Temp.indexes = new unsigned [l1+l2+1];
  unsigned * idtemp=Temp.indexes;


  // Fill the indexes array

  unsigned * fin = set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp=fin-idtemp;

  Temp.coefficients = new double [ltemp+1];
  double * vatempg=Temp.coefficients;

  Temp.length = ltemp+1;


  // Fill the coefficients array

 for (unsigned i=0;i<ltemp;i++)
    {
      unsigned a=pu1-id1;
      unsigned b=pu2-id2;

      if (a==l1|| id1[a]!=idtemp[i])
	{
	  vatempg[i]=cvalue*va2[b];  // cvalue*va2[b]+(P.cvalue)*0
	  pu2++;
	  continue;
	}

      if (b==l2 || id2[b]!=idtemp[i])
	{
	  vatempg[i]=(P.cvalue)*va1[a];  // cvalue*0+(P.cvalue)*va1[a]
	  pu1++;
	  continue;
	}

      vatempg[i]=cvalue*va2[b]+(P.cvalue)*va1[a];
      pu1++;
      pu2++;
    }


 // Compute the error
 // in a new noise symbol

 Temp.indexes[ltemp]=inclast();
 Temp.coefficients[ltemp]=rad()*(P.rad());

 return Temp;

}


// Operator  /
// It's a non affine-operation
// We use the identity x/y = x * (1/y

AAF AAF::operator / (const AAF & P)
{
  return (*this)*inv(P);
}


// Square root operator
// It's a non affine-operation
// We use the Chebyshev approximation

AAF sqrt(const AAF & P)
{

  double a, b;

  // sqrt(x) is approximated by f(x)=alpha*x+dzeta
  // delta is the maximum absolute error
  double alpha, dzeta, delta;

  double t;  // temporary var

  a=P.convert().getlo(); // [a,b] is our interval
  b=P.convert().gethi();

  t= (sqrt(a)+sqrt(b));


  alpha=1/t; // alpha is the slope of the line r(x) that
             // interpolate (a, sqrt(a)) and (b, f(b))

  // dzeta calculation:
  dzeta= (t/8)+0.5*(sqrt(a*b))/t;

  // Calculation of the error
  delta=(sqrt(b)-sqrt(a))*(sqrt(b)-sqrt(a))/(8*t);


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


 // Compute the error
 // in a new noise symbol

  // zk = delta

  Temp.indexes[P.length]=Temp.inclast();   // the error indx
  Temp.coefficients[P.length]=delta;


  return Temp;

}


// Inverse (1/x) operator
// It's a non-affine operation
// We use mini-range approximation
// cause undershoot can be high with Chebyshev here

AAF inv(const AAF & P)
{


  double a, b;
  double alpha, dzeta, delta;

  double t1, t2;  // temporary var


  a=P.convert().getlo();
  b=P.convert().gethi();


  // a := min(abs(a), abs(b))
  // b := max(abs(a), abs(b))

  t1=fabs(a);
  t2=fabs(b);

  a= t1 <? t2;  // min(t1,t2)
  b= t1 >? t2;  // max(t1,t2)

  // Derivative of 1/x is -1/x*x

  alpha=-1/(b*b);

  Interval i((1/a)-alpha*a,1/b-alpha*b);
  dzeta = i.mid();

  if ((P.convert().getlo()) < 0) dzeta = -dzeta;

  delta=i.radius();



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


 // Compute the error
 // in a new noise symbol

  // zk = delta

  Temp.indexes[P.length]=Temp.inclast();   // the error indx
  Temp.coefficients[P.length]=delta;


  return Temp;

}


// Power function
// only for integer exposents

AAF pow(const AAF & P, int exp)
{

  AAF Temp = P;

  if (!exp)
    {
      Temp = 1;
    }

  else if (exp > 0)
  {
      for (unsigned i=1; i<exp; i++)
	Temp = Temp*P;
  }
  else
  {
      for (unsigned i=1; i<fabs(exp); i++)
	Temp = Temp*P;
      Temp = inv(Temp);
  }

 return Temp;

}
