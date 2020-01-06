/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Docking
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Docking_h
#define Docking_h 1

#include "cpu_time.h"
#include "parallel.h"
#include "parameter.h"
#include "receptor.h"
#include "ligand.h"
#include "fft_process.h"

using namespace std;

class Docking
{
private:
    Docking(Docking &c) {}
    const Docking & operator=(const Docking &c);
    CPUTime   *_cputime;
    Parallel  *_parallel;
    Parameter *_parameter;
    Receptor  *_receptor;
    Ligand    *_ligand;
    FFTProcess    *_fft_process;
    int       _Num_grid;
    float     *_Grid_coord;
    float     **_Mol_coord;
    float     *_Fwork;
    int       *_Iwork;
    size_t       _Memfw;
    size_t       _Memiw;
protected:
    virtual void  maxsize_voxel();
    virtual void  alloc_array(const int &maxatom, const int &nag, const size_t &ng3);
    virtual void  create_voxel(Protein *rprotein, size_t myid2);
    virtual void  ligand_rotationz(float *theta, size_t myid2);
public:
    Docking(CPUTime *pcputime,Parallel *pparallel,Parameter *pparameter,
            Receptor *rreceptor,Ligand *rligand)
        : _cputime(pcputime),_parallel(pparallel),_parameter(pparameter),
          _receptor(rreceptor),_ligand(rligand) {
#ifdef DEBUG
        cout << "Constructing Docking.\n";
#endif
    }
    virtual ~Docking() {
#ifdef DEBUG
        cout << "Destructing Docking.\n";
#endif
        delete [] _Grid_coord;
        delete [] _Mol_coord;
        delete [] _Fwork;
        delete [] _Iwork;
        delete _fft_process;
    }
    virtual void  initialize();
    virtual void  rec_init();
    virtual void  dockz();
    virtual void  dock_memory_free();
    virtual void  output();
    virtual void  output_detail(); // for analysis
    virtual void  output_calc_time_log(); // for analysis
};

#endif
