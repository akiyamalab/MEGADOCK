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

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Parallel
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Parallel_h
#define Parallel_h 1

#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

class Parallel
{
private:
    Parallel(Parallel &c) {}
    size_t       _Nproc2;
    int       _Num_gpu; 
public:
    Parallel(const size_t &nproc2) : _Nproc2(nproc2) {
#ifdef DEBUG
        cout << "Constructing Parallel.\n";
#endif
    }
    virtual   ~Parallel()     {
#ifdef DEBUG
        cout << "Destructing Parallel.\n";
#endif
    }
    void      nproc2(int i)       {
        _Nproc2 = i;
    }
    size_t       nproc2()        {
        return _Nproc2;
    }
    void      num_gpu(int i)       {
        _Num_gpu = i;
    }
    int       num_gpu()        {
        return _Num_gpu;
    }
};

#endif
