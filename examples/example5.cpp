/*
 * example5.cpp -- Time performances between AA and IA model
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

/* This is a program to compare precison between
 * the AA model and the IA model
 *
 * Results are in seconds, with an accuracy of 10 ms
 *
 * We use libEasyval by Johan Vervloet and Stefan Becuwe for the IA model
 * 
 * This example is written for a function y=f(x1,x2) but it can be easyly changed
 * to use with a higher number of variables
 *
 * Usage: ./example5 LOWER_BOUND1 UPPER_BOUND1 LOWER_BOUND2 UPPER_BOUND2 BOXN
 *
 * It means [LOWER_BOUNDi,UPPERBOUNDi] is the interval of variables xi divided
 * in BOXN sud-divisions (boxes)
 *
 * If you want to change the function to display, you can change it in the
 * template fct eval_fct()
 *
 * Compile : g++ -laa lEasyval example5.cpp -o example5
 * (requires libaa and libEasyval)
 *
 * (C) 2003 Olivier Gay <olivier.gay@epfl.ch>
*/

#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sys/times.h>
#include <unistd.h>
#include <aa.h>
#include <Easyval.hh>

using namespace std;


#define PI 4*(atan(1.0))



// Put here the function you want to evaluate

template <typename TP> TP eval_fct(TP x1, TP x2)
{
  TP y;
  y=1+(x1*x1-2)*x2+x1*x2*x2;;
  return y;
}


int main(int argc, char **argv)
{

  clock_t tstart1, tstart2;
  clock_t tstop1, tstop2;
  struct tms tmsave;

  double total1, total2;

  // interval of x1

  double lbound1, ubound1;
  double width1; // subdivision width

  //interval of x2

  double lbound2, ubound2;
  double width2; // subdivision width

  unsigned boxn; // number of boxes

  double ydelta; // the y error

  double y; 
 

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

  cout << "Time performances between AA and IA models" << endl; 

  // AA

  double x1_a, x1_b;
  double x2_a, x2_b;


  tstart1 = times(&tmsave);

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
      y = v.getcenter();

    }

  tstop1 = times(&tmsave);

  total1 = (double) (tstop1-tstart1)/100;

  // IA

  tstart2 = times(&tmsave);

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
    }

  tstop2 = times(&tmsave);

  total2 = (double) (tstop2-tstart2)/100;
  

  cout << "AA: "<< total1 << "s " << endl;
  cout << "IA: "<< total2 << "s " << endl;

  exit(0);


}





















































