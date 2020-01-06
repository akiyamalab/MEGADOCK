/*
 * Copyright (C) 2020 Tokyo Institute of Technology
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
