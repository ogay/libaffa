/* 
 * aa_rounding.h -- Platform dependent rounding
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

#ifndef AA_ROUNDING_H
#define AA_ROUNDING_H

//#ifdef HAVE_IEEEFP_H
//#include <ieeefp.h>
//#define EV_UPWARD FP_RP
//#define EV_DOWNWARD FP_RM
//#define fegetround fpgetround
//#define fesetround fpsetround
//#ifdef HAVE_FP_RND_T
//typedef fp_rnd_t aa_rnd_t;
//#else
//typedef fp_rnd aa_rnd_t;
//#endif
//#else
#include <fenv.h>
#define AA_UPWARD FE_UPWARD
#define AA_DOWNWARD FE_DOWNWARD
typedef unsigned int aa_rnd_t;
//#endif
//#endif

unsigned int aa_fesetround(aa_rnd_t);
aa_rnd_t aa_fegetround(void);

#endif
