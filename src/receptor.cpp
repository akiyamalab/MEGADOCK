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

#include "receptor.h"
#include <algorithm>
#include <string.h>

//============================================================================//
void Receptor::rpscace(const float &aceratio, const int &num_grid, float *grid_coord,
                       float *atom_coord_rotated, int *iwork, float *fwork,
                       const float &param_rec_core, const float &param_lig_core,const int &old_voxel_flag)
//============================================================================//
{
    const float rho = -3.0;                 // MEGADOCK parameter
    const float epsilon = param_rec_core;   // Tuning Parameter
    const float open_space = -7777.0;       // temporary score for open space

    int       *nearesta;
    int       *surf_atom;
    int       *counter;
    float     *grid_r, *grid_i;
    float     *rjudge2, *radius_core2, *radius_surf2;
    float     *distance3;
    float     *xd, *yd, *zd;

    const int ng1 = num_grid;
    const int ng2 = ng1 * ng1;
    const int ng3 = ng1 * ng1 * ng1;
    const int na  = num_atoms();
    const int nag = num_atoms() * ng1;

    float *voxel_rpscace;
    voxel_rpscace = new float[ng3];

	memset(voxel_rpscace, 0.0, sizeof(float)*ng3);
	memset(fwork, 0.0, sizeof(float)*ng3);

    nearesta  = iwork;
    surf_atom = &iwork[ng3];
    counter   = &iwork[ng3+na];
    grid_r    = fwork;
    grid_i    = &fwork[ng3];
    rjudge2   = &fwork[ng3*2];
    radius_core2  = &fwork[ng3*2+na];
    radius_surf2  = &fwork[ng3*2+na*2];
    distance3     = &fwork[ng3*2+na*3];
    xd        = &fwork[ng3*3+na*3];
    yd        = &fwork[ng3*3+na*3+nag];
    zd        = &fwork[ng3*3+na*3+nag*2];

    precalc(ng1,nearesta,rjudge2,radius_core2,radius_surf2);
    voxel_init(ng1,grid_coord,atom_coord_rotated,grid_r,grid_i,distance3,xd,yd,zd);


// ##############################
// #fill up grid by atom-sphere #
// ##############################

    for( int l = 0 ; l < na ; l++ ) {// atom-grid assignment
        //const int search_range = 2;
        const int search_range = (2.4 + grid_width -0.01) / grid_width;
        const int lc = ng1 * l;
        const int i2 = atom_coord_rotated[3*l  ] / grid_width + ng1 / 2;
        const int j2 = atom_coord_rotated[3*l+1] / grid_width + ng1 / 2;
        const int k2 = atom_coord_rotated[3*l+2] / grid_width + ng1 / 2;
        const int ia = max(i2 - search_range, 0);
        const int ja = max(j2 - search_range, 0);
        const int ka = max(k2 - search_range, 0);
        const int ib = min(i2 + search_range+1, ng1);
        const int jb = min(j2 + search_range+1, ng1);
        const int kb = min(k2 + search_range+1, ng1);
        
        for( int i = ia ; i < ib ; i++ ) {// grid around atom[l]
            if(xd[lc+i] > radius_core2[l]) continue;
            
            for( int j = ja ; j < jb ; j++ ) {
                const float d2 = xd[lc+i]+yd[lc+j];
                if(d2 > radius_core2[l]) continue;
                const int ij=ng2*i+ng1*j;
                
                for( int k = ka ; k < kb ; k++ ) {
                    const float d3 = xd[lc+i] + yd[lc+j] + zd[lc+k];
                    if( d3 < radius_core2[l] ) {// distance(grid-atom) < Van der Waals radius (* core)
                        grid_r[ij+k] = epsilon;    // grid[i] is filled up by atom[l]
                    }
                }
            }
        }
    }

// #####################
// #fill up cavity     #
// #####################
    
    for (int ij=0; ij<ng2; ij++) {// fill up cavity (phase1, 1/2)
        const int i = ij % ng1;
        const int j = ij / ng1;
        
        for (int k=0; k<ng1; k++) {// X, forward
            if(grid_r[i+ng1*j+ng2*k]!=epsilon) {
                grid_r[i+ng1*j+ng2*k]=open_space;
            } else {
                break;
            }
        }
        
        for (int k=ng1-1; k>=0; k--) {// X, backward
            if(grid_r[i+ng1*j+ng2*k]!=epsilon) {
                grid_r[i+ng1*j+ng2*k]=open_space;
            } else {
                break;
            }
        }
        
        for (int k=0; k<ng1; k++) {// Y, forward
            if(grid_r[k+ng1*i+ng2*j]!=epsilon) {
                grid_r[k+ng1*i+ng2*j]=open_space;
            } else {
                break;
            }
        }
        
        for (int k=ng1-1; k>=0; k--) {// Y, backward
            if(grid_r[k+ng1*i+ng2*j]!=epsilon) {
                grid_r[k+ng1*i+ng2*j]=open_space;
            } else {
                break;
            }
        }
        
        for (int k=0; k<ng1; k++) {// Z, forward
            if(grid_r[j+ng1*k+ng2*i]!=epsilon) {
                grid_r[j+ng1*k+ng2*i]=open_space;
            } else {
                break;
            }
        }
        
        for (int k=ng1-1; k>=0; k--) {// Z, backward
            if(grid_r[j+ng1*k+ng2*i]!=epsilon) {
                grid_r[j+ng1*k+ng2*i]=open_space;
            } else {
                break;
            }
        }
    }
    
    for (int ink=0; ink<2; ink++) {// fill up cavity (phase1, 2/2)

        for (int ijk=ng2; ijk<ng3-ng2; ijk++) {
            //const int i = ijk / ng2;
            const int j = (ijk / ng1) % ng1;
            if(j<2||j>ng1-2) continue; // skip border
            const int k = ijk % ng1;
            if(k<2||k>ng1-2) continue; // skip border
            
            if(grid_r[ijk]==open_space) {
                if(grid_r[ijk+1]==0) {
                    grid_r[ijk+1]=open_space;
                }
                if(grid_r[ijk-1]==0) {
                    grid_r[ijk-1]=open_space;
                }
                if(grid_r[ijk+ng1]==0) {
                    grid_r[ijk+ng1]=open_space;
                }
                if(grid_r[ijk-ng1]==0) {
                    grid_r[ijk-ng1]=open_space;
                }
                if(grid_r[ijk+ng2]==0) {
                    grid_r[ijk+ng2]=open_space;
                }
                if(grid_r[ijk-ng2]==0) {
                    grid_r[ijk-ng2]=open_space;
                }
            }
        }
    }
    

    for(int i=0; i<ng3; i++) {// fill up cavity (phase2)
        if(grid_r[i]==0) {
            grid_r[i]=epsilon;
        }
    }
    
    for(int i=0; i<ng3; i++) {// restore open_space to 0 (phase3)
        if(grid_r[i]==open_space) {
            grid_r[i]=0;
        }
    }
    
    
// #################################
// # set surface score of receptor #
// #################################

    for( int l = 0 ; l < na ; l++ ) {// atom-grid assignment
        //const int search_range = 5;
        const int search_range = (12 + grid_width -0.01) / grid_width;
        const int lc = ng1 * l;
        const int i2 = atom_coord_rotated[3*l  ] / grid_width + ng1 / 2;
        const int j2 = atom_coord_rotated[3*l+1] / grid_width + ng1 / 2;
        const int k2 = atom_coord_rotated[3*l+2] / grid_width + ng1 / 2;
        const int ia = max(i2 - search_range, 0);
        const int ja = max(j2 - search_range, 0);
        const int ka = max(k2 - search_range, 0);
        const int ib = min(i2 + search_range+1, ng1);
        const int jb = min(j2 + search_range+1, ng1);
        const int kb = min(k2 + search_range+1, ng1);
        
        for( int i = ia ; i < ib ; i++ ) {// grid around atom[l]
            if(xd[lc+i] > rjudge2[l]) continue;
            
            for( int j = ja ; j < jb ; j++ ) {
                const float d2 = xd[lc+i]+yd[lc+j];
                if(d2 > rjudge2[l]) continue;
                const int ij = ng2*i+ng1*j;
                
                for( int k = ka ; k < kb ; k++ ) {
                    const int ijk = ij+k;
                    const float d3 = xd[lc+i] + yd[lc+j] + zd[lc+k];
                    
                    if( d3 < rjudge2[l]) {// distance(grid-atom) < Van der Waals radius + cutoff D(=3.6ang)
                        //counter[ijk] ++;    // grid[i] is filled up by atom[l]
                        voxel_rpscace[ijk] += 1.0 + _ACE[l] * (-0.8) * aceratio;
                        // rPSC value(1.0) + ACE value (0.8 is weight parameter) * w_h
                    }
                }
            }
        }
    }
    
    for( int ijk = 0 ; ijk < ng3 ; ijk++ ) {
        if( grid_r[ijk]==epsilon ) continue;
        
        grid_r[ijk] = voxel_rpscace[ijk];
    }

    delete[] voxel_rpscace;

    return;
}


//============================================================================//
void Receptor::electro(const float &beta,const float &eratio,
                       const int &num_grid,float *grid_coord,float *atom_coord_rotated,
                       int *iwork,float *fwork,const int &old_voxel_flag)
//============================================================================//
{
    const float   ftr1    = 6.0;      // FTDock parameter
    const float   ftr2    = 8.0;      // FTDock parameter
    const float   er1 = 4.0;          // FTDock parameter
    const float   er2 = 38.0;         // FTDock parameter
    const float   er3 = -224.0;       // FTDock parameter
    const float   er4 = 80.0;         // FTDock parameter

    float     *grid_r;
    float     *grid_i;
    const float   *xd, *yd, *zd;

    const int ng1 = num_grid;
    const int ng2 = ng1 * ng1;
    const int ng3 = ng1 * ng1 * ng1;
    const int na  = num_atoms();
    const int nag = num_atoms() * ng1;

    int search_range;
    float     er, d3;

    grid_r    = fwork;
    grid_i    = &fwork[ng3];
    xd        = &fwork[ng3*3+na*3];
    yd        = &fwork[ng3*3+na*3+nag];
    zd        = &fwork[ng3*3+na*3+nag*2];


    const float grid_width = grid_coord[1] - grid_coord[0];
    const float charge_per_grid = -1 * grid_width * er4 / (beta * eratio);//formula derived from e(r) (r>8A)
    
    for( int l = 0 ; l < na ; l++ ) {
        const float abs_charge = fabs(_Charge[l]);
        if( abs_charge <= EPS ) continue;
        
        if (abs_charge < 0.07) { //hardcode (Gridsize = 1.2A)
            search_range = 4;
        } else if (abs_charge < 0.21) {
            search_range = 5;
        } else {
            search_range = (int)(abs_charge / charge_per_grid);
        }
        
        const int lc = ng1 * l;
        const int i2 = atom_coord_rotated[3*l  ] / grid_width + ng1 / 2; //limit search range
        const int j2 = atom_coord_rotated[3*l+1] / grid_width + ng1 / 2;
        const int k2 = atom_coord_rotated[3*l+2] / grid_width + ng1 / 2;
        const int ia = max(i2 - search_range, 0);
        const int ja = max(j2 - search_range, 0);
        const int ka = max(k2 - search_range, 0);
        const int ib = min(i2 + search_range+1, ng1);
        const int jb = min(j2 + search_range+1, ng1);
        const int kb = min(k2 + search_range+1, ng1);
        
        for( int i = ia ; i < ib ; i++ ) {// grid around atom[l]
            const int ijk_a = ng2 * i;
            
            for( int j = ja ; j < jb ; j++ ) {
                const int ijk_b = ijk_a + ng1 * j;
                const float d2 = xd[lc+i]+yd[lc+j];
                
                for( int k = ka ; k < kb ; k++ ) {
                    const int ijk = ijk_b + k;
                    
                    if( grid_r[ijk] < 0.0 ) continue; //this grid is core of receptor
                    
                    d3 = d2 + zd[lc+k];
                    d3 = sqrt(d3);
                    
                    if( d3 >= ftr2 ) {
                        er = er4;
                    } else if( d3 > ftr1 ) {
                        er = er2*d3 + er3;
                    } else {
                        er = er1;
                    }
                    const float elec = eratio * beta * _Charge[l] / er / d3;
                    grid_i[ijk] += elec;
                }
            }
        }
    }
    
    return;
}

//============================================================================//
void Receptor::precalc(const int &num_grid,int *nearesta,float *rjudge2,
                       float *radius_core2,float *radius_surf2)
//============================================================================//
{
    const float   D = 3.6;            
    const float   rcore2 = 1.5;       
    const float   rsurf2 = 0.8;       

    const int na  = num_atoms();
    const int ng1 = num_grid;
    const int ng3 = ng1 * ng1 * ng1;

    for( int i = 0 ; i < na ; i++ ) {
        rjudge2[i]      = (_Radius[i] + D)*(_Radius[i] + D);
        radius_core2[i] = _Radius[i] * _Radius[i] * rcore2;
        radius_surf2[i] = _Radius[i] * _Radius[i] * rsurf2;
    }

    for( int i = 0 ; i < ng3 ; i++ ) {
        nearesta[i] = -99;
    }

    return;
}
