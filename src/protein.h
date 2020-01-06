/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Protein
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Protein_h
#define Protein_h 1

#include "pdb_entry.h"
#include "parameter.h"

using namespace std;

class Protein : public PDBEntry     // Protein class
{
private:
    friend class  Docking;
    friend class  Receptor;
    friend class  Ligand;
    friend class  FFTProcess;
    //  const Protein & operator=(const Protein &c);
    int       _Num_grid;
    float     grid_width;
    float     _Edge[3][2];
    float     _Center[3];
    float     _Angle[3];
    float     *_Radius;
    float     *_Charge;
    float     *_ACE;
    int       *_Atomtype;

protected:
    virtual void  axis_edge();
    virtual void  voxel_init(const int &num_grid,float *grid_coord,
                             float *mol_coord,float *grid_r,float *grid_i,
                             float *distance3,float *xd,float *yd,float *zd);
public:
    Protein(string &rinput_file) : PDBEntry(rinput_file) {
#ifdef DEBUG
        cout << "Constructing Protein.\n";
#endif
    }
    virtual ~Protein() {
#ifdef DEBUG
        cout << "Destructing Protein.\n";
#endif
        delete [] _Radius;
        delete [] _Charge;
        delete [] _ACE;
        delete [] _Atomtype;
    }
    virtual void  initialize(Parameter *rparameter);
    virtual void  shift_center(float *mol_coord);
    virtual float edge(const int &i,int j) {
        return _Edge[i][j];
    }
    virtual float center(const int &i) {
        return _Center[i];
    }
    virtual void  center(float *coord) {
        for( int i = 0 ; i < 3 ; i++ ) {
            _Center[i] = coord[i];
        }
        return;
    }
    virtual float angle(const int &i) {
        return _Angle[i];
    }
    virtual void  angle(float *angle) {
        for( int i = 0 ; i < 3 ; i++ ) {
            _Angle[i] = angle[i];
        }
        return;
    }
    virtual int   num_grid() {
        return _Num_grid;
    }
    virtual void  num_grid(const int &num_grid) {
        _Num_grid = num_grid;
        return;
    }

    virtual void  rpscace(const float &aceratio, const int &num_grid,float *grid_coord,float *atom_coord_rotated,
                          int *iwork,float *fwork,
                          const float &param_rec_core, const float &param_lig_core,const int &old_voxel_flag) {
        ;
    }
    virtual void  electro(const float &beta,const float &eratio,
                          const int &num_grid,float *grid_coord,float *atom_coord_rotated,
                          int *iwork,float *fwork,const int &old_voxel_flag) {
        ;
    }
    virtual void  precalc(const int &num_grid,int *nearesta,float *rjudge2,
                          float *radius_core2,float *radius_surf2) {
        ;
    }
};

#endif
