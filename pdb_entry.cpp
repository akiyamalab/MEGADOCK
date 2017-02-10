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
//  Class Name : PDBEntry
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include "pdb_entry.h"

//============================================================================//
void PDBEntry::pdb_read()
//============================================================================//
{
    string    pdb, atomtype;
    vector<float> coordinate0[3];

#ifdef DEBUG
    cout << "Input file = " << _Input_file << endl;
#endif

    ifstream Input(_Input_file.c_str(),ios::in);
    if( !Input ) {
        cerr << "[ERROR] PDB file [" << _Input_file << "] not open!!" << endl;
        exit(1);
    }

    while(1) {
        if( !getline(Input,pdb) ) break;
#ifdef DEBUG2
        cout << pdb << endl;
#endif

        if( pdb.substr(0,6) == "ATOM  " || pdb.substr(0,6) == "HETATM" ) {
            if( pdb.substr(76,2) == " H" ) {
                printf("# This PDB file contains hydrogen atom. MEGADOCK skips hydrogen atom line.\n");
            } else {
                coordinate0[0].push_back(atof(pdb.substr(30,8).c_str()));
                coordinate0[1].push_back(atof(pdb.substr(38,8).c_str()));
                coordinate0[2].push_back(atof(pdb.substr(46,8).c_str()));

                atomtype = pdb.substr(11,10);
                _Atom_type.push_back(erase_space(atomtype));
            }
        }
    }

    Input.close();

    _Num_atoms = coordinate0[0].size();

    _Coordinate = new float[3*_Num_atoms];
    if( !_Coordinate ) {
        cerr << "[ERROR] Out of memory. Number of atomic coordinates*3(x,y,z) = ("
             << 3*_Num_atoms << ") for (_Coordinate) in pdb_entry.cpp!!\n";
        exit(1);
    }

    for( int i = 0 ; i < _Num_atoms ; i++ ) {
        float coord[3];

        for( int j = 0 ; j < 3 ; j++ ) {
            coord[j] = coordinate0[j][i];
        }
        coordinate(i,coord);
    }

    return;
}

//============================================================================//
string PDBEntry::erase_space(const string &s0)
//============================================================================//
{
    int       n;
    string    s1;

    s1 = s0;

    while(1) {
        n = s1.find("  ");

        if( n == (int) std::string::npos ) {
            break;
        } else {
            s1.erase(n,1);
        }
    }

    if( s1.substr(0,1) == " " ) {
        s1.erase(0,1);
    }

    if( s1.substr(s1.size()-1,1) == " " ) {
        s1.erase(s1.size()-1,1);
    }

    return s1;
}
