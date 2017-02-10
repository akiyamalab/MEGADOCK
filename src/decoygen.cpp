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


//==================================================================//
//
//  Software Name : MEGADOCK (decoygen tool)
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//  Last update: December 2, 2014
//
//==================================================================//

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct parse_res {
    float r1, r2, r3, l1, l2, l3, 
          a1, a2, a3, t1, t2, t3, 
          rand1, rand2, rand3, spacing, score;
    int N;
    char receptor[100], ligand[100];
};

void rotateAtom(float, float, float, float *, float *, float *,
                float, float, float );
void createPDB(char *, char *, float, float, float, float, float,
              float, float, float, float, float, float, float, 
              int, int, int, int, float);
struct parse_res parse_out(char *, int);

//--------------------------------------------------
// Function rotateAtom 
//--------------------------------------------------
void rotateAtom (float oldX, float oldY, float oldZ,
                 float *newX, float *newY, float *newZ,
                 float psi, float theta, float phi )
{
    float r11, r21, r31, r12, r22, r32, r13, r23, r33;
    r11 = cos(psi)*cos(phi)  -  sin(psi)*cos(theta)*sin(phi);
    r21 = sin(psi)*cos(phi)  +  cos(psi)*cos(theta)*sin(phi);
    r31 = sin(theta)*sin(phi);
    
    r12 = -cos(psi)*sin(phi)  -  sin(psi)*cos(theta)*cos(phi);
    r22 = -sin(psi)*sin(phi)  +  cos(psi)*cos(theta)*cos(phi);
    r32 = sin(theta)*cos(phi);
    
    r13 = sin(psi)*sin(theta);
    r23 = -cos(psi)*sin(theta);
    r33 = cos(theta);
    
    *newX = r11 * oldX + r12 * oldY + r13 * oldZ;
    *newY = r21 * oldX + r22 * oldY + r23 * oldZ;
    *newZ = r31 * oldX + r32 * oldY + r33 * oldZ;
    
} 

//--------------------------------------------------
// Function createPDB
//--------------------------------------------------
void createPDB(char *out, char *in, float rand1, float rand2, float rand3,
               float r1, float r2, float r3, float l1, float l2, float l3, 
               float a1, float a2, float a3, int t1, int t2, int t3, int N,
               float spacing)
{
    FILE *newfp, *oldfp;
    char liner[85], coord1[10], coord2[10], coord3[10], lastline[100], tmp[85];
    float tx1, ty1, tz1, tx2, ty2, tz2, oldx, oldy, oldz;
    
    newfp = fopen(out, "w");
    oldfp = fopen(in, "r");
    while ((fgets(liner, 31, oldfp)) != NULL){
        strcpy(tmp, liner);
        if(!(strncmp(tmp, "ATOM  ", 6)) || !(strncmp(tmp, "HETATM", 6))){
            fgets(coord1, 9, oldfp); 
            oldx = atof(coord1) - l1;
            fgets(coord2, 9, oldfp); 
            oldy = atof(coord2) - l2;
            fgets(coord3, 9, oldfp); 
            oldz = atof(coord3) - l3;
            fgets(lastline, 60, oldfp);
            rotateAtom(oldx, oldy, oldz, &tx1, &ty1, &tz1, rand1, rand2, rand3);
            rotateAtom(tx1, ty1, tz1, &tx2, &ty2, &tz2, a1, a2, a3);
            // adjust so coordinates are in the box
            if (t1 >= N/2) t1 -= N;
            if (t2 >= N/2) t2 -= N;  
            if (t3 >= N/2) t3 -= N;  
            
            /* write to the file using grid-adjusted coords */
            fprintf (newfp, "%s%8.3f%8.3f%8.3f%s", liner, tx2-t1*spacing+r1, 
                     ty2-t2*spacing+r2, tz2-t3*spacing+r3, lastline);
        }else if(!(strncmp(tmp, "TER   ", 6)) || !(strncmp(tmp, "ENDMDL", 6)) || !(strncmp(tmp, "END   ", 6))) {
            fgets(lastline, 60, oldfp);
            fprintf(newfp, "%s%s", liner, lastline);
        } 
    } 
    fclose(oldfp);
    fclose(newfp);
}

//--------------------------------------------------
// Function parse_out
//--------------------------------------------------
struct parse_res parse_out(char *out_name, int line_num){
    struct parse_res ret;
    FILE *fp;
    fp = fopen(out_name, "r");
    fscanf(fp, "%d %f", &ret.N, &ret.spacing); /* parse 4 lines of .out  */
    fscanf(fp, "%f %f %f", &ret.rand1, &ret.rand2, &ret.rand3);
    fscanf(fp, "%s %f %f %f", ret.receptor, &ret.r1, &ret.r2, &ret.r3);
    fscanf(fp, "%s %f %f %f", ret.ligand, &ret.l1, &ret.l2, &ret.l3);
    int i;
    for (i = 0; i < line_num; i++){
        fscanf(fp, "%f %f %f %f %f %f %f", &ret.a1, &ret.a2, &ret.a3, &ret.t1, &ret.t2, &ret.t3, &ret.score);
    }
    fclose(fp);
    return ret;
}


//--------------------------------------------------
int main(int argc, char **argv )  
{ 

    char new_lig[256], old_lig[256], outfile[256]; 
    int line_num;
    if (argc != 5){
        printf("++ decoygen tool ++\n");
        printf("\n%s [decoy_filename] [used_ligand.pdb] [.outfile] [decoy no.]", argv[0]);
        printf("\n  **Note: files must be in PDB format.\n");
        exit(-1);
    }
    
    strcpy(new_lig, argv[1]);
    strcpy(old_lig, argv[2]);
    strcpy(outfile, argv[3]);
    line_num = atof(argv[4]);
    struct parse_res ret = parse_out(outfile, line_num);
    createPDB(new_lig, old_lig, ret.rand1, ret.rand2, ret.rand3, 
              ret.r1, ret.r2, ret.r3, ret.l1, ret.l2, ret.l3, 
              ret.a1, ret.a2, ret.a3, ret.t1, ret.t2, ret.t3, ret.N, ret.spacing);
    
    return 0;
    
} 
