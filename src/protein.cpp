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

#include "protein.h"
#include "string.h"

//============================================================================//
void Protein::initialize(Parameter *rparameter)
//============================================================================//
{
    pdb_read();
    axis_edge();
    const int na = num_atoms();
    grid_width = rparameter->grid_width;
    

    // for radius
    _Radius = new float[na];

    for( int i = 0 ; i < na ; i++ ) {
        _Radius[i] = rparameter->atom_radius(atom_type(i));
    }

    // for charge
    if( rparameter->_Score_func == 2 || rparameter -> _Score_func == 3 ) {
        _Charge = new float[na];
        for( int i = 0 ; i < na ; i++ ) {
            _Charge[i] = rparameter->atom_charge(atom_type(i));
        }
    }

    // for ACE
    if ( rparameter -> _Score_func == 3 ) {
        _ACE = new float[na];
        for ( int i = 0; i < na; i++) {
            _ACE[i] = rparameter->atom_ace(atom_type(i));
        }
    }
    
    _Atomtype = new int[na];
    for( int i = 0 ; i < na ; i++ ) { // 1:CA, 2:C,N, 0:others
    	if(!strcmp(atom_type(i).substr(0,3).c_str(),"CA ")){
        	_Atomtype[i] = 1;
        }else if( (!strcmp(atom_type(i).substr(0,2).c_str(),"C ")) || (!strcmp(atom_type(i).substr(0,2).c_str(),"N ")) ){
        	_Atomtype[i] = 2;
    	}else{
    		_Atomtype[i] = 0;
    	}
        //printf(" %5d, %10s, %2d\n",i,atom_type(i).c_str(),_Atomtype[i]);
    }

    return;
}

//============================================================================//
void Protein::axis_edge()
//============================================================================//
{
    float coord[3];

    for( int i = 0 ; i < 3 ; i++ ) {
        _Edge[i][0] = BIG;
        _Edge[i][1] = -BIG;
    }

    for( int i = 0 ; i < num_atoms() ; i++ ) {
        for( int j = 0 ; j < 3 ; j++ ) {
            coord[j] = coordinate(i,j);
        }
#ifdef DEBUG2
        cout << "i = " << i << " X = " << coord[0] << " Y = " << coord[1]
             << " Z = " << coord[2] << endl;
#endif
        for( int j = 0 ; j < 3 ; j++ ) {
            if( coord[j] < _Edge[j][0] ) {
                _Edge[j][0] = coord[j];
            }
            if( coord[j] > _Edge[j][1] ) {
                _Edge[j][1] = coord[j];
            }
        }
    }

    for( int i = 0 ; i < 3 ; i++ ) {
        _Center[i] = (_Edge[i][0] + _Edge[i][1]) / 2.0;
    }

    for( int i = 0 ; i < 3 ; i++ ) {
        _Angle[i] = 0.0;
    }

#ifdef DEBUG2
    for( int i = 0 ; i < 3 ; i++ ) {
        cout << "_Edge[" << i << "][0] = " << _Edge[i][0] << endl;
        cout << "_Edge[" << i << "][1] = " << _Edge[i][1] << endl;
    }
#endif

    return;
}

//============================================================================//
void Protein::shift_center(float *mol_coord)
//============================================================================//
{
    float     coord[3];
    const int na = num_atoms();

    for( int i = 0 ; i < 3 ; i++ ) {
        coord[i] = ( _Edge[i][0] + _Edge[i][1] ) / 2.0;
    }

    for( int i = 0 ; i < na ; i++ ) {
        for( int j = 0 ; j < 3 ; j++ ) {
            mol_coord[i*3+j] = coordinate(i,j) - coord[j];
        }
    }

    return;
}

//============================================================================//
void Protein::voxel_init(const int &num_grid,float *grid_coord,float *mol_coord,
                         float *grid_r,float *grid_i,float *distance3,float *xd,
                         float *yd,float *zd)
//============================================================================//
{
    const int na  = num_atoms();
    const int nag = num_atoms() * num_grid;
    const int ng3 = num_grid * num_grid * num_grid;

	memset(grid_r, 0.0, sizeof(float)*ng3);
	memset(grid_i, 0.0, sizeof(float)*ng3);
	memset(distance3, BIG, sizeof(float)*ng3);


    for( int l = 0, k = 0 ; l < na ; l++ ) {
        const int l3 = l * 3;

        for( int i = 0 ; i < num_grid ; i++ ) {
            xd[k  ] = mol_coord[l3  ] - grid_coord[i];
            yd[k  ] = mol_coord[l3+1] - grid_coord[i];
            zd[k++] = mol_coord[l3+2] - grid_coord[i];
        }
    }

    for( int i = 0 ; i < nag ; i++ ) {
        xd[i] *= xd[i];
        yd[i] *= yd[i];
        zd[i] *= zd[i];
    }

    return;
}
