/*
 * example2.cpp -- This example generates datas and commands of the AA model
 * to gnuplot to display
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

/* This is a program to display AA
 *
 * You'll need gnuplot to display the datas
 *
 * It writes 3 files
 * a) gnuplot_commands : several commands to plot the function and the boxes
 * b) data1   : the datas of the evaluated function
 * c) data2   : the datas for the boxes of our AAF
 *
 * Usage: ./example2 LOWER_BOUND UPPER_BOUND BOXN [function]
 *
 * It means we use an [LOWER_BOUND:UPPERBOUND] interval divided
 * in BOXN sud-divisions (boxes) and [function] is an optionnaly
 * function to display in the gnuplot titles
 * 
 * If you want to change the function to display, you can change it in the
 * template fct eval_fct()
 *
 * In the gnuplot display a star after the title indicate
 * that a missing box in the plot is due to a non-finite value
 * (either infinite or nan) during the caculation  
 *
 * Otherwise a missing box indicate a very small height
 * of the box.
 *
 * Compile : g++ -laa example2.cpp -o example2
 *
 * (C) 2003 Olivier Gay <olivier.gay@epfl.ch>
*/

#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <aa.h>

using namespace std;


#define PI 4*(atan(1.0))
#define PTS 50       // number of points for the interpolated function
                     // 50 is normally enough


// Put here the function you want to evaluate

template <typename TP> TP eval_fct(TP x)
{
  TP y;
  y=sqrt(x*x-x+0.5)/sqrt(x*x+0.5); // y(x)
  y=sqrt(y*y-y+0.5)/sqrt(y*y+0.5); // y(y(x))

  return y;
}


int main(int argc, char **argv)
{

  // 0 1 24 ; -1 1

  double lbound, ubound; // lower/upper bound of the interval
  unsigned boxn; // number of boxes

  double width; // width of a box
  double ydelta; // the y error

  //double ydelta_tot;// sum of ydelta for aa

  double x, y; 
  bool fin;  
  bool finy, finydelta;

  double min, max;

  Interval itv;
  AAF u, v;  // v = v(u)


  if (argc < 4) 
    {
      cout << "Usage: " << argv[0] << " LOWER_BOUND UPPER_BOUND BOXN [function]" << endl;
      return 0;
    }

  lbound = atof(argv[1]);
  ubound = atof(argv[2]);
  boxn = atoi(argv[3]);

  if (lbound > ubound)
    {
      cerr << "The lower bound must be smaller than the upper bound!" << endl;
      return -1;
    }

  width = (ubound-lbound)/boxn;

  cout << "Generate datas files and commands file:" << endl << endl;

  ofstream data1("data1");
  ofstream data2("data2");
  ofstream gp("gnuplot_commands");

  // We write the evaluated function
  // Actually we can have similary result by directly ploting
  // the function in gnuplot with "plot" but with this way
  // we can re-use the template of eval_fct()

  data1 << "# Datas of the evaluated function" << endl;
  data1 << "#\tx\ty" << endl;

  for (unsigned i=0; i <= PTS; i++)
    {
      x=lbound+i*(ubound-lbound)/PTS;
      y=eval_fct(x);

      data1 << "\t" << x << "\t";

      // We must test for an overflow or such
      // so if y is either infinite or NaN
      // we put "?" which is understood by gnuplot as
      // missing data

      if (finite(y))
	data1 << y << endl;
      else
	data1 << "?" << endl;

    }

  min=max=y;

  data1.close();

  // We write the error boxes for the AAF

  double x1,x2, xc; 

  fin = 1;

  data2 << "# Datas of the error boxes" << endl;
  data2 << "# for Affine Arithmetic" << endl;
  data2 << "#\tx\ty\tydelta" << endl;

  for (unsigned i=0; i <= boxn-1; i++)
    {
      x1=lbound+i*width;
      x2=x1+width;
      xc=(x1+x2)/2;

      itv.modlohi(x1,x2);
      u=itv;
      v=eval_fct(u);

      ydelta = v.rad();
      y = v.getcenter();

      data2 << "\t" << xc << "\t";

      // We check if we calculated a non-finite value
      // and we save the sup box and the inf box

      finy = finite(y);
      finydelta = finite(ydelta);

      if (finy)
	data2 << y << "\t";
      else
	data2 << "?" << "\t";

      if (finydelta)
	data2 << ydelta << endl;
      else
	  data2 << "?" << endl;

      if (finy && finydelta)
	{
	  if (min > y-ydelta) min=y-ydelta;
	  if (max < y+ydelta) max=y+ydelta;
	}
      else
	fin=0;

    }


  data2.close();

  cout << "AA: [" << min << ":" << max << "], " << max-min << endl << endl;


  // The commands file of gnuplot

  gp << "set grid" << endl;
  gp << "set nokey" << endl;
  gp << "set title \"Affine Arithmetic representation";
  if (!fin) gp << "*";
  if (argc==5) gp << " : " << argv[4];
  gp <<  "\" " << endl;
  gp << "set xrange [" << lbound << ":" << ubound << "]" << endl;
  gp << "plot \"data2\" using 1:2:(" << width/2 << "):3 with boxxyerrorbars 3, ";
  gp << "\"data1\" smooth csplines 1" << endl;  
  gp.close();

  cout << "Execute plot.sh to display the result" << endl;

  exit(0);

}

















