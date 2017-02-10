/*
 * Copyright (C) 2014 Tokyo Institute of Technology
 *
 *
 * This file is part of MEGADOCK.
 * MEGADOCK is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MEGADOCK is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MEGADOCK.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

//==================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : (Constant)
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//==================================================================//

#ifndef Constant_h
#define Constant_h 1

#include <math.h>

using namespace std;

const float     EPS         = 1.0e-6;
const float     BIG         = 1.0e+9;
const float     PI      = 4.0*atan(1.0);
const float     RAD         = PI/180.0;
const float     clockpersec_float = (float)CLOCKS_PER_SEC;

#endif
