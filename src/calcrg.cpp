/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */


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

void createRg(char *, float, float, float, float, float, float,
              float, float, float, int, int, int, int, float, int, int);

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
// Function createRg
//--------------------------------------------------
void createRg(float rand1, float rand2, float rand3,
              float r1, float r2, float r3, float l1, float l2, float l3,
              float a1, float a2, float a3, int t1, int t2, int t3,
              int N, float spacing, int isPDB, int nl)
{
    float tx1, ty1, tz1, tx2, ty2, tz2, oldx, oldy, oldz;
    oldx = 0.0; oldy = 0.0; oldz = 0.0;
    rotateAtom(oldx, oldy, oldz, &tx1, &ty1, &tz1, rand1, rand2, rand3);
    rotateAtom(tx1, ty1, tz1, &tx2, &ty2, &tz2, a1, a2, a3);
    if (t1 >= N/2) t1 -= N;
    if (t2 >= N/2) t2 -= N;
    if (t3 >= N/2) t3 -= N;

    if ( isPDB == 0 )
        printf ("%7.3f, %7.3f, %7.3f\n", tx2-t1*spacing+r1, ty2-t2*spacing+r2, tz2-t3*spacing+r3);
    else if ( isPDB == 1 )
        printf ("HETATM %4d  S   RGS P%4d    %8.3f%8.3f%8.3f  1.00  1.00           S\n", nl, nl, 
                tx2-t1*spacing+r1, ty2-t2*spacing+r2, tz2-t3*spacing+r3);

}

//--------------------------------------------------
int main(int argc, char **argv )  
{ 

    char outfile[256]; 
    int isPDBflag;
    if (argc != 3){
        printf("++ center points of ligand ++\n");
        printf("\n%s [.outfile] [isPDB? 0: csv,  1: pseudo PDB]\n", argv[0]);
        exit(-1);
    }
    
    strcpy(outfile, argv[1]);
    isPDBflag = atoi(argv[2]);
    
    struct parse_res ret;
    FILE *fp;
    fp = fopen(outfile, "r");
    fscanf(fp, "%d %f", &ret.N, &ret.spacing); /* parse 4 lines of .out  */
    fscanf(fp, "%f %f %f", &ret.rand1, &ret.rand2, &ret.rand3);
    fscanf(fp, "%s %f %f %f", ret.receptor, &ret.r1, &ret.r2, &ret.r3);
    fscanf(fp, "%s %f %f %f", ret.ligand, &ret.l1, &ret.l2, &ret.l3);

    int i = 1;
    while((fscanf(fp, "%f %f %f %f %f %f %f", &ret.a1, &ret.a2, &ret.a3, &ret.t1, &ret.t2, &ret.t3, &ret.score)) != EOF){
        createRg(ret.rand1, ret.rand2, ret.rand3, ret.r1, ret.r2, ret.r3, ret.l1, ret.l2, ret.l3,
                 ret.a1, ret.a2, ret.a3, ret.t1, ret.t2, ret.t3, ret.N, ret.spacing, isPDBflag, i);
        i++;
        if(i == 10000)
            i = 0;
    }
    fclose(fp);
    return 0;
} 
