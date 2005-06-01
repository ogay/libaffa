/*
 * aa_aafarithm.cpp -- Affine arithmetical operations
 * Copyright (c) 2005 Nathan Hurst
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

#ifndef AA_UTIL
#define AA_UTIL

#include "aa_aaf.h"
#define handle_infinity(x) {\
	if((x).get_special() == AAF_TYPE_NAN) return AAF(AAF_TYPE_NAN); \
	if((x).get_special() == AAF_TYPE_INFINITE) return AAF(interval(-INFINITY, INFINITY)); \
	}

inline static AAF_TYPE binary_special(AAF_TYPE a, AAF_TYPE b) {
    if((a == AAF_TYPE_AFFINE) &&
       (b == AAF_TYPE_AFFINE)) {
        return AAF_TYPE_AFFINE;
    } else if((a == AAF_TYPE_NAN) ||
            (b == AAF_TYPE_NAN)) {
        return AAF_TYPE_NAN;
    } else if((a == (AAF_TYPE)(AAF_TYPE_NAN | AAF_TYPE_AFFINE)) ||
            (b == (AAF_TYPE)(AAF_TYPE_NAN | AAF_TYPE_AFFINE))) {
        return (AAF_TYPE)(AAF_TYPE_NAN | AAF_TYPE_AFFINE);
    }
    return (AAF_TYPE)(a | b);
}

#endif
