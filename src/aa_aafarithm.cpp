/*
 * aa_aafarithm.cpp -- Affine arithmetical operations
 * Copyright (c) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
 * Copyright (c) 2004 LIRIS (University Claude Bernard Lyon 1)
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


// Operator +
// Addition is an affine operation

AAF AAF::operator + (const AAF & P)
{
  unsigned l1 = length;
  unsigned l2 = P.length;

  unsigned * id1=indexes;
  unsigned * id2=P.indexes;

  double * va1=coefficients;
  double * va2=P.coefficients;

  unsigned * pu1=id1;
  unsigned * pu2=id2;

  AAF Temp(cvalue+P.cvalue);  // Create our resulting AAF

  Temp.indexes = new unsigned [l1+l2]; // the indexes of the result
  unsigned * idtemp=Temp.indexes;


  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp=fin-idtemp;

  Temp.coefficients = new double [ltemp];
  double * vatempg=Temp.coefficients;

  Temp.length = ltemp;


  // Fill the coefficients array
  // of the resulting AAF

 for (unsigned i=0;i<ltemp;i++)
    {
      unsigned a=pu1-id1;
      unsigned b=pu2-id2;

      if (a==l1|| id1[a]!=idtemp[i])
	{
	  vatempg[i]=va2[b];  // va2[b]+0
	  pu2++;
	  continue;
	}

      if (b==l2 || id2[b]!=idtemp[i])
	{
	  vatempg[i]=va1[a];  // va1[a]+0
	  pu1++;
	  continue;
	}

      vatempg[i]=va1[a]+va2[b];
      pu1++;
      pu2++;
    }

    return Temp;
}


// Operator -

AAF AAF::operator - (const AAF & P)
{

  unsigned l1 = length;
  unsigned l2 = P.length;

  unsigned * id1=indexes;
  unsigned * id2=P.indexes;

  double * va1=coefficients;
  double * va2=P.coefficients;

  unsigned * pu1=id1;
  unsigned * pu2=id2;

  AAF Temp(cvalue-P.cvalue);  // Create our resulting AAF

  Temp.indexes = new unsigned [l1+l2];
  unsigned * idtemp=Temp.indexes;


  // Fill the resulting indexes array
  // by merging the 2 input indexes array

  unsigned * fin = std::set_union(id1,id1+l1,id2,id2+l2,idtemp);
  unsigned ltemp=fin-idtemp;

  Temp.coefficients = new double [ltemp];
  double * vatempg=Temp.coefficients;

  Temp.length = ltemp;


  // Fill the coefficients array
  // of the resulting AAF

 for (unsigned i=0;i<ltemp;i++)
    {
      unsigned a=pu1-id1;
      unsigned b=pu2-id2;

      if (a==l1|| id1[a]!=idtemp[i])
	{
	  vatempg[i]=-va2[b];  // 0-va2[b]
	  pu2++;
	  continue;
	}

      if (b==l2 || id2[b]!=idtemp[i])
	{
	  vatempg[i]=va1[a];  // va1[a]-0
	  pu1++;
	  continue;
	}

      vatempg[i]=va1[a]-va2[b];
      pu1++;
      pu2++;
    }

    return Temp;

}


// Unary operator

AAF AAF::operator - ()
{

  AAF Temp(*this);

  Temp.cvalue=-(Temp.cvalue);
  for (unsigned i=0; i<length; i++)
    Temp.coefficients[i]=-(Temp.coefficients[i]);

  return Temp;

}


// Mul by a constant (on right)
// Affine operation

AAF AAF::operator * (double cst)
{

  AAF Temp(*this);
  Temp.cvalue=cst*cvalue;

  for (unsigned i=0; i<length; i++)
    Temp.coefficients[i]=cst*(Temp.coefficients[i]);

  return Temp;

}


// -- Non member AAF functions --

// Mul by a constant (the left case)

AAF operator * (double cst, AAF P)
{
  AAF Temp(P);
  return Temp*cst;
}


// Add a constant (the left case)

AAF operator + (double cst, AAF P)
{

  AAF Temp = cst;
  return (Temp+P);

}


// Sub a constant (the left case)

AAF operator - (double cst, AAF P)
{

  AAF Temp = cst;
  return (Temp-P);
}
