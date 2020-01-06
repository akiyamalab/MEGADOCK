/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Control
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Control_h
#define Control_h 1

#include "cpu_time.h"
#include "parallel.h"
#include "parameter.h"
#include "receptor.h"
#include "ligand.h"
#include "docking.h"

using namespace std;

class Control
{
private:
    Control(Control &c) {}
    const Control & operator=(const Control &c);
    CPUTime   *_cputime;
    Parallel  *_parallel;
    Parameter *_parameter;
    Receptor  *_receptor;
    Ligand    *_ligand;
    Docking   *_docking;
protected:
    virtual void  gridtable_11base_normal(int &ngrid,vector<int> &ngrid_table);
    virtual void  gridtable_13base_normal(int &ngrid,vector<int> &ngrid_table);
    virtual void  gridtable_07base_normal(int &ngrid,vector<int> &ngrid_table);
    virtual void  gridtable_fftw_custom(int &ngrid,vector<int> &ngrid_table);
    virtual void  gridtable_cufft_custom(int &ngrid,vector<int> &ngrid_table);
    virtual void  autogridr(const int &ngrid,vector<int> &ngrid_table);
    virtual void  autogridl(const int &ngrid,vector<int> &ngrid_table);
    virtual void  checkgridr();
    virtual void  checkgridl();
public:
    Control(CPUTime *pcputime,Parallel *pparallel)
        : _cputime(pcputime),_parallel(pparallel) {
#ifdef DEBUG
        cout << "Constructing Control.\n";
#endif
    }
    virtual ~Control() {
#ifdef DEBUG
        cout << "Destructing Control.\n";
#endif
    }
    virtual void  initialize(int argc,char *argv[]);
    virtual void  execute();
};

#endif
