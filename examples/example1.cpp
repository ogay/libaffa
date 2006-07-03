/*
 * example1.cpp -- A simple example of using the lib
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with libaffa; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <aa.h>
#include <cstdio>
#include <iostream>

using namespace std;

int main()
{
  interval u(-2,2);
  interval v(-1,1);
  AAF x = u;
  AAF r = v;
  AAF s = v;

  AAF temp1 = (10+x+r);
  AAF temp2 = (10-x+s);

  AAF z = (10+x+r)*(10-x+s);

  cout << " x" << endl;
  cout << x;
  cout << " r" << endl;
  cout << r;
  cout << " s" << endl;
  cout << s;

  cout << " (10+x+r)" << endl;
  cout << temp1;
  cout << " (10-x+s)" << endl;
  cout << temp2;

  cout << " (10+x+r)*(10-x+s)" << endl;
  cout << z;
  cout << z.convert();

  exit(0);
}

