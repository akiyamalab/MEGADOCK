/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Parameter
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Parameter_h
#define Parameter_h 1

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cassert>
#include <unistd.h>

#include "constant.h"
#include "parallel.h"

using namespace std;

class Parameter
{
private:
    friend class          Control;
    friend class          Docking;
    friend class          FFTProcess;
    friend class          Protein;
    friend class          HeaderDB;
    Parameter(Parameter &c) {}
    const Parameter & operator=(const Parameter &c);
    Parallel          *_parallel;
    string            _RecPDB_file;
    string            _LigPDB_file;
    string            _RLOut_file;
    string            _RLOut_file_detail; 
    string            _RLOut_file_csv;
    string			  calc_id;
    int               _IO_flag[3];
    int               detail_output_flag;
    int               calc_time_log_output_flag;

    int               _Num_grid;
    int               _Num_fft;
    int               _Num_fft_flag;
    int				  _Num_atom_max;
    int               _Num_output;
    int               _Num_output_flag;
    int               _Num_thread_limit;
    int               _Num_GPU_limit;

    int               _Score_func;
    int               _Num_sort;
    float             _Elec_ratio;
    float             _ACE_ratio;
    float             grid_width;
    float             ligand_max_edge;
    int               _Rotation_angle_set;
    int               fft_base_set;
    int               lig_elec_serial_flag;
    int				  fft_library_type;

    int				  tem_flag1;
    int				  tem_flag2;
    int				  tem_flag3;
    int				  tem_flag4;
    int               f1_flag;
    int               f2_flag;
    
    int               _Old_voxel_flag;
    float             _Grid_space_rec;
    float             _Grid_space_lig;
    //rPSC tuning
    float             _rPSC_param_rec_core;
    float             _rPSC_param_lig_core;

    float             *_Zangle;
    int               _Num_rot_angles;
    map<string,float>     _Charmmr;
    map<string,float>     _Charmmc;
    map<string,float>     _ACE;
protected:
    virtual void          usage();
    virtual void          default_param();
    virtual void          pdb_step();
    virtual void          parameter_set();
    virtual void          charmm_radius();
    virtual void          charmm_charge();
    virtual void          ace_value();
    virtual void          dangle_rot1();
    virtual void          dangle_rot3();
    virtual void          dangle_rot24();
    virtual void          dangle_rot360();
    virtual void          dangle_rot3600();
    virtual void          dangle_rot54000();

public:
    Parameter(Parallel *pparallel) : _parallel(pparallel) {
#ifdef DEBUG
        cout << "Constructing Parameter.\n";
#endif
    }
    virtual           ~Parameter() {
#ifdef DEBUG
        cout << "Destructing Parameter.\n";
#endif
        delete [] _Zangle;
    }
    virtual void          initialize(int argc,char *argv[]);
    virtual float         atom_radius(const string &atype);
    virtual float         atom_charge(const string &atype);
    virtual float         atom_ace(const string &atype);
};

#endif
