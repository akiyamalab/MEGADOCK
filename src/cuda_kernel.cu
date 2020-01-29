/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  cuda_kernel.cu
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include <stdio.h>
#define   FMAX(a,b)  ( ((a)>(b) ) ? (a) : (b) )
#define   FMIN(a,b)  ( ((a)>(b) ) ? (b) : (a) )

__global__ void lig_vox_fill(int ng1
                             ,int na
                             ,float delta
                             ,float *radius2
                             ,float *xd
                             ,float *yd
                             ,float *zd
                             ,float *grid_coord
                             ,float *atom_coord_rotated
                             ,float *grid_r
    						 ,float grid_width)
{
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;
    const int ng2 = ng1 * ng1;
    //const int ng3 = ng2 * ng1;

    if(id < na) {
        //const int search_range = 2;
        const int search_range = (2.4 + grid_width -0.01) / grid_width;
        const int lc = ng1 * id;
        const int id3 = id * 3;
        const int i2 = atom_coord_rotated[id3  ] / grid_width + ng1 / 2;
        const int j2 = atom_coord_rotated[id3+1] / grid_width + ng1 / 2;
        const int k2 = atom_coord_rotated[id3+2] / grid_width + ng1 / 2;
        const int ia = FMAX(i2 - search_range, 0);
        const int ja = FMAX(j2 - search_range, 0);
        const int ka = FMAX(k2 - search_range, 0);
        const int ib = FMIN(i2 + search_range+1, ng1);
        const int jb = FMIN(j2 + search_range+1, ng1);
        const int kb = FMIN(k2 + search_range+1, ng1);

        for( int i = ia ; i < ib ; i++ ) {// grid around atom[l]
            if(xd[lc+i] > radius2[id]) continue;
            for( int j = ja ; j < jb ; j++ ) {
                const float d2 = xd[lc+i]+yd[lc+j];
                if(d2 > radius2[id]) continue;
                const int ij = ng2*i+ng1*j;
                for( int k = ka ; k < kb ; k++ ) {
                    const float d3 = d2 + zd[lc+k];
                    if( d3 < radius2[id] ) {// distance(grid-atom) < Van der Waals radius (* core)
                        grid_r[ij+k] = delta;    // grid[i] is filled up by atom[l]
                    }
                }
            }
        }
    }
    //*/
}


__global__ void lig_rotation(int na, float *theta, float *atom_coord_orig, float *mole_center_coord, float *atom_coord_rotated)
{
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;

    const float r11 = cos(theta[0])*cos(theta[2])  -  sin(theta[0])*cos(theta[1])*sin(theta[2]);
    const float r21 = sin(theta[0])*cos(theta[2])  +  cos(theta[0])*cos(theta[1])*sin(theta[2]);
    const float r31 = sin(theta[1])*sin(theta[2]);
    const float r12 = -cos(theta[0])*sin(theta[2])  -  sin(theta[0])*cos(theta[1])*cos(theta[2]);
    const float r22 = -sin(theta[0])*sin(theta[2])  +  cos(theta[0])*cos(theta[1])*cos(theta[2]);
    const float r32 = sin(theta[1])*cos(theta[2]);
    const float r13 = sin(theta[0])*sin(theta[1]);
    const float r23 = -cos(theta[0])*sin(theta[1]);
    const float r33 = cos(theta[1]);

    if(id < na) {
        const int id3 = id * 3;
        float     x, y, z;

        x = atom_coord_orig[id3  ] - mole_center_coord[0];
        y = atom_coord_orig[id3+1] - mole_center_coord[1];
        z = atom_coord_orig[id3+2] - mole_center_coord[2];
        atom_coord_rotated[id3  ] = r11 * x + r12 * y + r13 * z;
        atom_coord_rotated[id3+1] = r21 * x + r22 * y + r23 * z;
        atom_coord_rotated[id3+2] = r31 * x + r32 * y + r33 * z;
    }
}


__global__ void lig_calc_dis_atomgrid(int na, int ng1, float *xd, float *yd, float *zd, float *grid_coord, float *atom_coord_rotated)
{
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;
    const int nag = na * ng1;
    if(id < nag) {
        const int cur_atom = (id / ng1);
        const int cur_atom3 = cur_atom * 3;
        const int cur_grid = id % ng1;
        xd[id] = atom_coord_rotated[cur_atom3  ] - grid_coord[cur_grid];
        yd[id] = atom_coord_rotated[cur_atom3+1] - grid_coord[cur_grid];
        zd[id] = atom_coord_rotated[cur_atom3+2] - grid_coord[cur_grid];
        xd[id] *= xd[id];
        yd[id] *= yd[id];
        zd[id] *= zd[id];
    }
}

__global__ void lig_vox_init_grid(int ng3,float *grid_r,float *grid_i)
{
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;
    if(id < ng3) { //initialize
        grid_r[id]=0.0;
        grid_i[id]=0.0;
    }
}

__global__ void lig_vox_init_fft(int nf3,cufftComplex *lig_in)
{
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;
    if(id < nf3) { //initialize
        lig_in[id] =  make_cuComplex( 0.0, 0.0);
        //lig_in[id].x=0.0;
        //lig_in[id].y=0.0;
    }
}

__global__ void ligand_voxel_set(int ng1
                                 ,cufftComplex *lig_in
                                 ,float *grid_r
                                 ,float *grid_i)
{
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;
    const int ng2 = ng1 * ng1;
    const int ng3 = ng2 * ng1;
    const int nf1 = ng1 * 2;
    const int nf2 = nf1 * nf1;
    const int ng1_half = ng1 / 2;
    const float surface = 1.0;        // grid-assignment score (protein surface)

    //if(id == 0) for(int i=0;i<ng3;i++) if(grid_i[i]!=0.0) printf(" [%03d,%03d,%03d] :  %6.3f\n",i/ng2,i/ng1%ng1,i%ng1,grid_i[i]);
    //if(id == 0) for(int i=0;i<ng3;i++) printf(" [%03d,%03d,%03d] :  %6.3f\n",i/ng2,(i/ng1)%ng1,i%ng1,grid_i[i]);

    if(id < ng3) {
        const int i = id / ng2;
        const int j = (id / ng1) % ng1;
        const int k = id % ng1;
        const int idoff = (i + ng1_half) * nf2 + (j + ng1_half) * nf1 + (k + ng1_half);

        //*
        if(grid_r[id]==surface) {// this condition judges whether surface(1.0) or temporary score(-8888.0)
            lig_in[idoff] =  make_cuComplex( grid_r[id], grid_i[id]);
        } else {
            lig_in[idoff] =  make_cuComplex( 0.0, grid_i[id]);
        }
        //*
    }
}


__global__ void lig_vox_surface_cut_CtoT(int ng1, float delta, float *grid_r)
{
    // Core score to Temporary score
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;
    const float swollen_surface = -8888.0;
    const int ng2 = ng1 * ng1;
    const int ng3 = ng2 * ng1;
    if(id < ng3) {
        const int i = id / ng2;
        const int j = (id / ng1) % ng1;
        const int k = id % ng1;
        if(i==0||i==ng1-1||j==0||j==ng1-1||k==0||k==ng1-1) { // skip border
        } else {
            if(grid_r[id]==delta) {
                if(grid_r[id-1]==0 ||
                        grid_r[id+1]==0 ||
                        grid_r[id-ng1]==0 ||
                        grid_r[id+ng1]==0 ||
                        grid_r[id-ng2]==0 ||
                        grid_r[id+ng2]==0) {
                    grid_r[id]=swollen_surface; 
                }
            }
        }
    }
}

__global__ void lig_vox_elec(int ng1,int na,float grid_width,float *_Charge,float *atom_coord_rotated,float *grid_i)
{
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;
    const int ng2 = ng1 * ng1;
    const float pad = (ng1 * grid_width / 2);

    //*
    if(id < na) {
        //if(!fabs(_Charge[id]) < 0.0001) continue;
        const int id3 = id * 3;
        //const int   i = (atom_coord_rotated[id3  ] + pad) / grid_width;
        //const int   j = (atom_coord_rotated[id3+1] + pad) / grid_width;
        //const int   k = (atom_coord_rotated[id3+2] + pad) / grid_width;
        const int   i = FMAX(0, FMIN((atom_coord_rotated[id3  ] + pad) / grid_width, ng1 - 1));
        const int   j = FMAX(0, FMIN((atom_coord_rotated[id3+1] + pad) / grid_width, ng1 - 1));
        const int   k = FMAX(0, FMIN((atom_coord_rotated[id3+2] + pad) / grid_width, ng1 - 1));

        //grid_i[i*ng2+j*ng1+k] += _Charge[id];
        //printf(" %08d-1 :  %.2f, %.2f\n",i*ng2+j*ng1+k,grid_i[i*ng2+j*ng1+k],_Charge[id]);
        atomicAdd(&grid_i[i*ng2+j*ng1+k],_Charge[id]);

        //printf(" %08d-2 :  %.2f, %.2f\n",i*ng2+j*ng1+k,grid_i[i*ng2+j*ng1+k],_Charge[id]);
    }
    //*/
}

__global__ void lig_vox_elec_serial(int ng1,int na,float grid_width,float *_Charge,float *atom_coord_rotated,float *grid_i)
{
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;
    const int ng2 = ng1 * ng1;
    const float pad = (ng1 * grid_width / 2);

    if(id==0) {
        for( int l = 0 ; l < na ; l++ ) {
            const int l3 = l*3;
            const int i =(atom_coord_rotated[l3  ] + pad) / grid_width;
            const int j =(atom_coord_rotated[l3+1] + pad) / grid_width;
            const int k =(atom_coord_rotated[l3+2] + pad) / grid_width;
            //printf(" [%5d] [x:%12.8f,y:%12.8f,z:%12.8f] [pad:%6.3f], [%3d,%3d,%3d] \n",l,atom_coord_rotated[l3  ],atom_coord_rotated[l3+1],atom_coord_rotated[l3+2],pad,i,j,k);
            //printf(" [%5d] [x:%8.0f,y:%8.0f,z:%8.0f] [pad:%6.3f], [%3d,%3d,%3d] \n",l,atom_coord_rotated[l3  ],atom_coord_rotated[l3+1],atom_coord_rotated[l3+2],pad,i,j,k);

            //if(grid_i[i*ng2+j*ng1+k]!=0)printf(" Pos : %d, current : %f, new : %f\n",i*ng2+j*ng1+k, grid_i[i*ng2+j*ng1+k], _Charge[l]);

            grid_i[i*ng2+j*ng1+k] += _Charge[l];
        }
    }
}


__device__ void lig_vox_surface_cut_TtoO(int ng3, float delta, float *grid_r)
{
    // Temporary score to Open space score
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;
    const float swollen_surface = -8888.0; // temporary score for swollen ligand surface
    if(id < ng3) {
        if(grid_r[id]==swollen_surface) { 
            grid_r[id]=0.0;
        }
    }
}

__global__ void convolution_gpu(int nf3, float *rec_r, float *rec_i, cufftComplex *lig_out, cufftComplex *lig_in)
{
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;

    if(id<nf3) {
        const float lig_r = lig_out[id].x;
        const float lig_i = lig_out[id].y;

        lig_in[id] =  make_cuComplex( rec_r[id]*lig_r + rec_i[id]*lig_i, rec_r[id]*lig_i - rec_i[id]*lig_r);
        //lig_in[id].x = rec_r[id]*lig_r + rec_i[id]*lig_i;
        //lig_in[id].y = rec_r[id]*lig_i - rec_i[id]*lig_r;
    }
}

__global__ void max_pos_single(int nf3, cufftComplex *out, float *score, int *pos)
{
    //blockDim.x = nThreads
    //score[nBlocks], pos[nBlocks] (nBlocks = nf3 / nThreads)
    //sdata[nThreads]
    extern __shared__ float sdata[];
    float mscore;

    const int thr_id  = threadIdx.x;
    const int nThreads = blockDim.x;
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;

    if(id < nf3) {
        mscore = sdata[thr_id] = cuCrealf(out[id])/nf3;
        __syncthreads();    //all threads set sdata[thr_id]

        //reduction
        for(int offset = nThreads / 2; offset > 0; offset /= 2) {
            if (thr_id < offset) {
                sdata[thr_id] = FMAX(sdata[thr_id],  sdata[thr_id +  offset]);
            }
            __syncthreads();
        }

        if (mscore == sdata[0]) {//mscore specify position of max score
            score[blockIdx.x] = sdata[0];
            pos[blockIdx.x] = id;
            //printf("   BLOCK ID:%d, sdata[0]=%f, pos=%d\n",blockIdx.x,sdata[0],id);
        }
    }
}

__global__ void max_pos_multi_set(int nf3, cufftComplex *out, float *temp_score, int *temp_index)
{
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;
    if(id < nf3) {
        temp_score[id] = out[id].x;
        temp_index[id] = id;
    }
}

//, std::vector<cufftComplex> *temp_result , thrust::vector<cufftComplex> *temp_result
//thrust::device_ptr<cufftComplex> *temp_result cufftComplex *temp_result,thrust::device_ptr<cufftComplex> temp_result
__global__ void max_pos_multi(int nf3, cufftComplex *out, float *score, int *pos,const int num_sort,const int offset)
{
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < offset) {
        if (out[id].x < out[id+offset].x) {
            out[id].x = out[id+offset].x;
            out[id].y = out[id+offset].y;
        }
        /*
        if(id==0) {
            for(int i=0; i<num_sort*2; i++) printf(" id=%d, %f %f\n",i,out[i].x/nf3,out[i].y);
            printf("\n");
        }
        //*/
    }
    //*/
}



/*
__global__ void max_pos_multi(int nf3, cufftComplex *out, float *score, int *pos,const int num_sort, float *temp_score, int *temp_index)
{
    //blockDim.x = nThreads,
    //score[nBlocks], pos[nBlocks] (nBlocks = nf3 / nThreads)
    //sdata[nThreads]
    extern __shared__ float sdata[];
    float mscore;
    int offset;

    const int thr_id  = threadIdx.x;
    const int nThreads = blockDim.x;
    const int id  = blockIdx.x * blockDim.x + threadIdx.x;

    //*
    if(id < nf3) {
        temp_score[id]=cuCrealf(out[id])/nf3;
        temp_index[id]=id;


        /*
        __syncthreads();    //all threads set sdata[thr_id]

        //reduction
        for(offset = nThreads / 2; offset > num_sort; ) {
            offset /= 2;
            if (thr_id < offset) {
                sdata[thr_id] = FMAX(sdata[thr_id],  sdata[thr_id +  offset]);
            }
            //if(id<1)printf(" id=%d, t=%d, off=%d\n",id,num_sort,offset);
            __syncthreads();
        }
        //if(id<1)printf(" [last] id=%d, t=%d, off=%d\n",id,num_sort,offset);

        //thrust::sort(sdata,sdata+10);

        if(id < num_sort) {
            if (mscore == sdata[id]) {//mscore specify position of max score (float equality comparison... amari yokunai)
                score[blockIdx.x] = sdata[0];
                pos[blockIdx.x] = id;
                //printf("   BLOCK ID:%d, sdata[0]=%f, pos=%d\n",blockIdx.x,sdata[0],i);
            }
        }
        //*
        if(temp_score[id] >3000) printf(" id=%d, %f %d\n",id,temp_score[id],temp_index[id]);
    }
    //*
}
//*/










