/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Receptor
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Receptor_h
#define Receptor_h 1

#include "protein.h"
#include "parameter.h"

using namespace std;

class Receptor : public Protein     // Receptor class
{
private:
    //  const Receptor & operator=(const Receptor &c);
public:
    Receptor(string &rinput_file) : Protein(rinput_file) {
#ifdef DEBUG
        cout << "Constructing Receptor.\n";
#endif
    }
    virtual ~Receptor() {
#ifdef DEBUG
        cout << "Destructing Receptor.\n";
#endif
    }
    virtual void  electro(const float &beta,const float &eratio,
                          const int &num_grid,float *grid_coord,float *atom_coord_rotated,
                          int *iwork,float *fwork,const int &old_voxel_flag);
    virtual void  rpscace(const float &aceratio, const int &num_grid, float *grid_coord,
                          float *atom_coord_rotated, int *iwork, float *fwork,
                          const float &param_rec_core, const float &param_lig_core,const int &old_voxel_flag);
    virtual void  precalc(const int &num_grid,int *nearesta,float *rjudge2,
                          float *radius_core2,float *radius_surf2);
};

#endif
