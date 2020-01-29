/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Ligand
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include "ligand.h"
#include <algorithm> 
#include <string.h>

//============================================================================//
void Ligand::rpscace(const float &aceratio, const int &num_grid,float *grid_coord,
                     float *atom_coord_rotated, int *iwork, float *fwork,
                     const float &param_rec_core, const float &param_lig_core,const int &old_voxel_flag)
//============================================================================//
{
    const float delta = param_lig_core;       // Tuning Parameter
    const float surface = 1.0;        // grid-assignment score (protein surface)
    const float swollen_surface = -8888.0; // temporary score for swollen ligand surface 
    
    int       *nearesta;
    int       *surf_atom;
    float     *grid_r, *grid_i;
    float     *rjudge2, *radius_core2, *radius_surf2;
    float     *distance3;
    float     *xd, *yd, *zd;

    const int ng1 = num_grid;
    const int ng2 = ng1 * ng1;
    const int ng3 = ng1 * ng1 * ng1;
    const int na  = num_atoms();
    const int nag = num_atoms() * ng1;

    nearesta  = iwork;
    surf_atom = &iwork[ng3];
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
//
// Reason for search_range = 2
//   search_range means the range of searching grids(around atom)
//   If search_range = 2 and nearest grid is (10,5,-2), search from (8,3,-4) to (12,7,0) as rectangular
//   In the 3 next grid, distance between grid and atom is 3A (2*1.2A + 1.2A/2)
//   2 is enough because max VdW radius (charmm) is 2.13A and rcore2 = 1.5, so 2.13A * sqrt(1.5) = 2.61A < 3A

    for( int l = 0 ; l < na ; l++ ) {// atom-grid assignment
        const int search_range = 2;
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
            const int ib=ng2*i;
            
            for( int j = ja ; j < jb ; j++ ) {
                const float d2 = xd[lc+i]+yd[lc+j];
                const int ij=ib+ng1*j;
                if(d2 > radius_core2[l]) continue;
                
                for( int k = ka ; k < kb ; k++ ) {
                    const float d3 = d2 + zd[lc+k];
                    if( d3 < radius_core2[l] ) {// distance(grid-atom) < Van der Waals radius (* core)
                        
                        grid_r[ij+k] = delta;    // grid[i] is filled up by atom[l]
                        
                    }
                }
            }
        }
    }

// ############################
// # scrape swollen surface   #
// ############################
//
// 1. scrape surface grid (1 or 2 layer) because protein surface is swollen due to radius_core2
// 2. re-assign score by radius_surf2 in [make protein surface]

    
    for(int i=1; i<ng1-1; i++) {
        const int ib = ng2 * i;
        for(int j=1; j<ng1-1; j++) {
            const int jb = ib + ng1 * j;
            for(int k=1, ijk=jb+1; k<ng1-1; k++,ijk++) {
                if(grid_r[ijk]==delta) {
                    if(grid_r[ijk-1]==0 ||
                       grid_r[ijk+1]==0 ||
                       grid_r[ijk-ng1]==0 ||
                       grid_r[ijk+ng1]==0 ||
                       grid_r[ijk-ng2]==0 ||
                       grid_r[ijk+ng2]==0) {
                        grid_r[ijk]=swollen_surface; 
                    }
                }
            }
        }
    }
    
    
    for(int ijk=ng2; ijk<ng3-ng2; ijk++) {
        if(grid_r[ijk]==swollen_surface) { 
            grid_r[ijk]=0.0;
        }
    }
    
// ##############################
// #make protein surface        #
// ##############################
    
    for( int l = 0 ; l < na ; l++ ) {// for each atom
        const int search_range = 2;
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
            if(xd[lc+i] > radius_surf2[l]) continue;
            const int ib=ng2*i;
	
            for( int j = ja ; j < jb ; j++ ) {
                const float d2 = xd[lc+i]+yd[lc+j];
                const int ij=ib+ng1*j;
                if(d2 > radius_surf2[l]) continue;
                
                for( int k = ka ; k < kb ; k++ ) {
                    const int ijk = ij+k;
                    if(grid_r[ijk]==delta) continue;
                    
                    const float d3 = d2 + zd[lc+k];
                    if( d3 < radius_surf2[l] ) {// distance(grid-atom) < Van der Waals radius (* surf)
                        grid_r[ijk] = surface;    // grid[i] is filled up by atom[l]
                    }
                }
            }
        }
    }
    
    
    if(surface!=delta) {// convert surface grid (next to open space) to [1]
        for(int ijk=ng2; ijk<ng3-ng2; ijk++) {
            const int i = ijk / ng2;
            const int j = (ijk / ng1) % ng1;
            const int k = ijk % ng1;
            
            if(i==0||i==ng1-1||j==0||j==ng1-1||k==0||k==ng1-1)continue; // skip border
            
            if(grid_r[ijk]==delta) {
                if(grid_r[ijk-1]==0 ||
                   grid_r[ijk+1]==0 ||
                   grid_r[ijk-ng1]==0 ||
                   grid_r[ijk+ng1]==0 ||
                   grid_r[ijk-ng2]==0 ||
                   grid_r[ijk+ng2]==0) {
                    grid_r[ijk]=surface;
                }
            }
        }
    }
    
    return;
}
//============================================================================//
void Ligand::electro(const float &beta,const float &eratio,
                     const int &num_grid,float *grid_coord,float *atom_coord_rotated,
                     int *iwork,float *fwork,const int &old_voxel_flag)
//============================================================================//
{
    const int ng1 = num_grid;
    const int ng2 = ng1 * ng1;
    const int ng3 = ng1 * ng1 * ng1;
    const int na = num_atoms();
    const int nag = na * ng1;

    float     *grid_i;
    grid_i    = &fwork[ng3];

    const float pad = (ng1 * grid_width / 2);

    for( int l = 0 ; l < na ; l++ ) { //for each atom, find nearest grid(i,j,k)
        //if(fabs(_Charge[l]) < 0.0001) continue;
        const int l3 = l * 3;
        //const int   i = (atom_coord_rotated[l3  ] + pad) / grid_width;
        //const int   j = (atom_coord_rotated[l3+1] + pad) / grid_width;
        //const int   k = (atom_coord_rotated[l3+2] + pad) / grid_width;
        const int   i = max(0, min(ng1 - 1, (int) ((atom_coord_rotated[l3  ] + pad) / grid_width)));
        const int   j = max(0, min(ng1 - 1, (int) ((atom_coord_rotated[l3+1] + pad) / grid_width)));
        const int   k = max(0, min(ng1 - 1, (int) ((atom_coord_rotated[l3+2] + pad) / grid_width)));
        
        //printf(" [%5d] [x:%12.8f,y:%12.8f,z:%12.8f] [pad:%6.3f], [%3d,%3d,%3d] \n",l,atom_coord_rotated[l3  ],atom_coord_rotated[l3+1],atom_coord_rotated[l3+2],pad,i,j,k);
        //printf(" [%5d] [x:%8.0f,y:%8.0f,z:%8.0f] [pad:%6.3f], [%3d,%3d,%3d] \n",l,atom_coord_rotated[l3  ],atom_coord_rotated[l3+1],atom_coord_rotated[l3+2],pad,i,j,k);
        /*
          if(grid_i[i*ng2+j*ng1+k]!=0){
          printf(" conflict: [%d,%d,%d] (%f > %f)\n",i,j,k,grid_i[i*ng2+j*ng1+k],_Charge[l]);
          }
        */
        grid_i[i*ng2+j*ng1+k] += _Charge[l];
    }
    //for(int i=0;i<ng3;i++) if(grid_i[i]!=0.0) printf(" [%03d,%03d,%03d] :  %6.3f\n",i/ng2,i/ng1%ng1,i%ng1,grid_i[i]);
    //for(int i=0;i<ng3;i++) if(grid_i[i]!=0.0) printf(" [%03d,%03d,%03d] :  %6.3f\n",i/ng2,(i/ng1)%ng1,i%ng1,grid_i[i]);
    
    return;
}

//============================================================================//
void Ligand::precalc(const int &num_grid,int *nearesta,float *rjudge2,
                     float *radius_core2,float *radius_surf2)
//============================================================================//
{
    const float   rcore2 = 1.5;        
    const float   rsurf2 = 1.0;        

    const int na  = num_atoms();
    const int ng1 = num_grid;
    const int ng3 = ng1 * ng1 * ng1;
    
    for( int i = 0 ; i < na ; i++ ) {
        radius_core2[i] = _Radius[i] * _Radius[i] * rcore2;
        radius_surf2[i] = _Radius[i] * _Radius[i] * rsurf2;
    }

    memset(nearesta, -99, sizeof(int)*ng3);

    return;
}

