/*
 * example4.cpp -- Precision performances between the IA and the AA model
 * Copyright (c) 2003 EPFL (Ecole Polytechnique Federale de Lausanne)
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

/* This is a program to compare precison between
 * the AA model and the IA model
 *
 * We use libEasyval by Johan Vervloet and Stefan Becuwe for the IA model 
 *
 * This example is written for a function y=f(x1,x2) but it can be easyly changed
 * to use with a higher number of variables
 *
 * Usage: ./example4 LOWER_BOUND1 UPPER_BOUND1 LOWER_BOUND2 UPPER_BOUND2 BOXN
 *
 * It means [LOWER_BOUNDi,UPPERBOUNDi] is the interval of variables xi divided
 * in BOXN sud-divisions (boxes)
 * 
 * If you want to change the function to display, you can change it in the
 * template fct eval_fct()
 *
 * Compile : g++ -laa lEasyval example4.cpp -o example4
 * (requires libaffa and libEasyval)
 *
 * (C) 2003 Olivier Gay <olivier.gay@a3.epfl.ch>
 */

#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <aa.h>
#include <Easyval.hh>


using namespace std;


#define PI 4*(atan(1.0))
#define PTS 50       // number of points for the interpolated function
                     // 50 is normally enough


// Put here the function you want to evaluate

template <typename TP> TP eval_fct(TP x1, TP x2)
{
  TP y;
  y=1+(x1*x1-2)*x2+x1*x2*x2;;
  return y;
}


int main(int argc, char **argv)
{

  // interval of x1

  double lbound1, ubound1;
  double width1; // subdivision width

  //interval of x2

  double lbound2, ubound2;
  double width2; // subdivision width

  unsigned boxn; // number of boxes

  double ydelta; // the y error

  double ydelta_tot1;// sum of ydelta for aa
  double ydelta_tot2;// sum of ydelta for ia 

  double y; 
  bool fin1, fin2;  
  bool finy, finydelta;

  double min1, max1;
  double min2, max2;

  interval itv1, itv2;
  AAF u1, u2, v;  // v = v(u1,u2)

  Easyval m1,m2, n;  // n = n(m1,m2)


  if (argc != 6) 
    {
      cout << "Usage: " << argv[0] << " LOWER_BOUND1 UPPER_BOUND1 LOWER_BOUND2 UPPER_BOUND2 BOXN" << endl;
      return 0;
    }

  lbound1 = atof(argv[1]);
  ubound1 = atof(argv[2]);

  lbound2 = atof(argv[3]);
  ubound2 = atof(argv[4]);



  boxn = atoi(argv[5]); //

  if ((lbound1 > ubound1)||(lbound2 > ubound2))
    {
      cerr << "The lower bound must be smaller than the upper bound!" << endl;
      return -1;
    }

  width1 = (ubound1-lbound1)/boxn;
  width2 = (ubound2-lbound2)/boxn;

  cout << "Generate statistics:" << endl << endl;

  min1=max1=eval_fct((lbound1+ubound1)/2,(lbound2+ubound2)/2);

  // Affine Arithmetic

  double x1_a, x1_b;
  double x2_a, x2_b;

  fin1 = 1;
  ydelta_tot1 = 0;


  for (unsigned i=0; i <= boxn-1; i++)
    {

      // x1
      x1_a=lbound1+i*width1;
      x1_b=x1_a+width1;
      itv1.modlohi(x1_a,x1_b);
      u1=itv1;

      // x2
      x2_a=lbound2+i*width2;
      x2_b=x2_a+width2;
      itv2.modlohi(x2_a,x2_b);
      u2=itv2;


      v=eval_fct(u1,u2);

      ydelta = v.rad();
      ydelta_tot1 += ydelta;

      y = v.get_center();

      // We check if we calculated a non-finite value
      // and we save the sup box and the inf box

      finy = finite(y);
      finydelta = finite(ydelta);

      if (finy && finydelta)
	{
	  if (min1 > y-ydelta) min1=y-ydelta;
	  if (max1 < y+ydelta) max1=y+ydelta;
	}
      else
	fin1=0;

    }

  // interval Arithmetic

  fin2 = 1;

  ydelta_tot2 = 0;

  min2 = min1;
  max2 = max1;

  for (unsigned i=0; i <= boxn-1; i++)
    {


      // x1
      x1_a=lbound1+i*width1;
      x1_b=x1_a+width1;
      m1.set(x1_a,x1_b);

      // x2
      x2_a=lbound2+i*width2;
      x2_b=x2_a+width2;
      m2.set(x2_a,x2_b);

      n=eval_fct(m1,m2);   

      y = n.midpoint();
      ydelta = n.getSup()-y; 
      ydelta_tot2 += ydelta;

      // We check if we calculated a non-finite value
      // and we save the sup box and the inf box

      finy = finite(y);
      finydelta = finite(ydelta);

      if (finy && finydelta)
	{
	  if (min2 > y-ydelta) min2=y-ydelta;
	  if (max2 < y+ydelta) max2=y+ydelta;
	}
      else
	fin2=0;

    }
 
  double area_AA = (max1-min1);
  double area_IA = (max2-min2);


  cout << "IA: [" << min2 << ":" << max2 << "], " << max2-min2 << ", 100% , 100%" << endl;
  cout << "AA: [" << min1 << ":" << max1 << "], " << max1-min1 << ", " << 100*(area_AA/area_IA) << "% , " << 100*(ydelta_tot1/ydelta_tot2) << "%" << endl;

  if (!(fin1 && fin2))
    {
      cout << "There are some infinite or NaN values so statistics are not reliable" << endl;
    }

  cout << endl;

  exit(0);

}

















