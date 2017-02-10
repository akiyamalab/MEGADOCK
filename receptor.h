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
