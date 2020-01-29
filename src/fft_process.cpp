/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : FFTProcess
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include "fft_process.h"

#include <string.h>
#include <unistd.h>

#define NUM_THREADS 512 //should be power of 2

#ifdef CUFFT
#include "cuda_kernel.cu"

#endif


//============================================================================//
void FFTProcess::alloc_array(const int &num_fft)
//============================================================================//
{
    //cout << "FFT::alloc_array |" <<num_fft<< endl; cout.flush();
    _Num_fft = num_fft;

    const size_t nf3       = _Num_fft * _Num_fft * _Num_fft;
    const int num_sort  = _parameter->_Num_sort;
    const int num_angle = _parameter->_Num_rot_angles;
    const int no        = _parameter->_Num_output;
    const size_t nproc2    = _parallel->nproc2();
    int   num_toprank;

    num_toprank = num_angle * num_sort;
    if( no > num_toprank ) num_toprank = no;

    alloc_fft();

    _Select.resize(nproc2);
    for( int i = 0 ; i < nproc2 ; i++ ) {
        _Select[i].resize(num_sort);
    }

    _Top.resize(num_toprank);

    //---------- memory allocation for _Current_rot_angle_num
    _Current_rot_angle_num = new int[nproc2];

    _cputime->record_malloc( sizeof(float)*nf3*2*(1 + nproc2));

    //---------- memory allocation for _FFT_rec_r
    _FFT_rec_r = new float[nf3];
    if( !_FFT_rec_r ) {
        cerr << "[ERROR] Out of memory. Number of listed receptors = ("
             << nf3 << ") for (_FFT_rec_r) in fft_process.cpp!!\n";
        exit(1);
    }

    //---------- memory allocation for _FFT_rec_i
    _FFT_rec_i = new float[nf3];
    if( !_FFT_rec_i ) {
        cerr << "[ERROR] Out of memory. Number of listed receptors = ("
             << nf3 << ") for (_FFT_rec_i) in fft_process.cpp!!\n";
        exit(1);
    }

    return;
}

//============================================================================//
void FFTProcess::alloc_fft()
//============================================================================//
{
    const int nf1 = _Num_fft;
    const size_t nf3 = _Num_fft * _Num_fft * _Num_fft;
    const size_t nproc2  = _parallel->nproc2();
    const int num_gpu = _parallel->num_gpu();
    const int na = _ligand->num_atoms();
    size_t myid2;

#ifdef CUFFT
    const int num_sort = _parameter->_Num_sort;
    const int ng1 = _Num_fft / 2;
    const int ng3 = ng1 * ng1 * ng1;
    const int nag = na * ng1;
    //for ligand voxelization on GPU
    const int nThreads = NUM_THREADS;
    const int nBlocks_nf3 = (nf3 + (nThreads-1)) / nThreads;

    CUFFTin_host  = new cufftComplex[nf3];
    CUFFTout_host = new cufftComplex[nf3];

    _cputime->record_malloc( sizeof(cufftComplex)*nf3*2 ); //_in/outBuf

    //printf(" start: %p\n",&CUFFTin_host[0].x);

    /*
    for( int i = 0 ; i < nf3; i++ ) { // This initialization should be executed only once
        if(i<20)printf(" %p %p\n",&CUFFTin_host[i].x,&CUFFTin_host[i].y);
        if(i>nf3-20)printf(" %p %p\n",&CUFFTin_host[i].x,&CUFFTin_host[i].y);
        CUFFTin_host[i] = make_cuComplex(0.0, 0.0);
        CUFFTin_host[i].x = 0.0;
        CUFFTin_host[i].y = 0.0;
    }
    //*/

    int lenCUFFTin_host = (int)(((long int)&CUFFTin_host[nf3-1].x) - ((long int)&CUFFTin_host[0].x) + sizeof(CUFFTin_host[nf3-1]))/sizeof(CUFFTin_host[nf3-1]);
    if(lenCUFFTin_host !=nf3) printf("# discontinuous memory allocation occurs\n");

    //printf("   end: %ld\n",(long long int)&CUFFTin_host[nf3-1].y - &CUFFTin_host[0].x);

    cufft_plan = new cufftHandle[num_gpu];
    cufft_result = new cufftResult[num_gpu];

    CUFFTin_gpu = new cufftComplex*[num_gpu];
    CUFFTout_gpu = new cufftComplex*[num_gpu];
    _FFT_rec_r_gpu = new float*[num_gpu];
    _FFT_rec_i_gpu = new float*[num_gpu];

    grid_r_gpu = new float*[num_gpu];
    grid_i_gpu = new float*[num_gpu];
    grid_coord_gpu = new float*[num_gpu];
    radius_core2_gpu = new float*[num_gpu];
    radius_surf2_gpu = new float*[num_gpu];
    _Charge_gpu = new float*[num_gpu];
    xd_gpu = new float*[num_gpu];
    yd_gpu = new float*[num_gpu];
    zd_gpu = new float*[num_gpu];
    atom_coord_rotated_gpu = new float*[num_gpu];
    atom_coord_orig_gpu = new float*[num_gpu];
    mole_center_coord_gpu = new float*[num_gpu];
    ligand_rotation_angle_gpu = new float*[num_gpu];
    top_score_gpu = new float*[num_gpu];
    top_index_gpu = new int*[num_gpu];
    top_score_host = new float*[num_gpu];
    top_index_host = new int*[num_gpu];


    for(int gpu_id = 0; gpu_id < num_gpu; gpu_id++) {
        cudaSetDevice(gpu_id);
        cufft_result[gpu_id] = cufftPlan3d(&cufft_plan[gpu_id], nf1, nf1, nf1, CUFFT_C2C);

        checkCudaErrors( cudaMalloc((void **)&CUFFTin_gpu[gpu_id],  sizeof(cufftComplex)*nf3) );
        checkCudaErrors( cudaMalloc((void **)&CUFFTout_gpu[gpu_id], sizeof(cufftComplex)*nf3) );
        checkCudaErrors( cudaMalloc((void **)&_FFT_rec_r_gpu[gpu_id], sizeof(float)*nf3) );
        checkCudaErrors( cudaMalloc((void **)&_FFT_rec_i_gpu[gpu_id], sizeof(float)*nf3) );

        checkCudaErrors( cudaMalloc((void **)&grid_r_gpu[gpu_id],  sizeof(float)*ng3));
        checkCudaErrors( cudaMalloc((void **)&grid_i_gpu[gpu_id],  sizeof(float)*ng3));
        checkCudaErrors( cudaMalloc((void **)&grid_coord_gpu[gpu_id],  sizeof(float)*ng1));
        checkCudaErrors( cudaMalloc((void **)&radius_core2_gpu[gpu_id],  sizeof(float)*na));
        checkCudaErrors( cudaMalloc((void **)&radius_surf2_gpu[gpu_id],  sizeof(float)*na));
        checkCudaErrors( cudaMalloc((void **)&_Charge_gpu[gpu_id],  sizeof(float)*na));
        checkCudaErrors( cudaMalloc((void **)&xd_gpu[gpu_id],  sizeof(float)*nag));
        checkCudaErrors( cudaMalloc((void **)&yd_gpu[gpu_id],  sizeof(float)*nag));
        checkCudaErrors( cudaMalloc((void **)&zd_gpu[gpu_id],  sizeof(float)*nag));
        checkCudaErrors( cudaMalloc((void **)&atom_coord_rotated_gpu[gpu_id],  sizeof(float)*na*3));
        checkCudaErrors( cudaMalloc((void **)&atom_coord_orig_gpu[gpu_id],  sizeof(float)*na*3));
        checkCudaErrors( cudaMalloc((void **)&mole_center_coord_gpu[gpu_id],  sizeof(float)*3));
        checkCudaErrors( cudaMalloc((void **)&ligand_rotation_angle_gpu[gpu_id],  sizeof(float)*3));
        checkCudaErrors( cudaMalloc((void **)&top_score_gpu[gpu_id], sizeof(float)*nBlocks_nf3*num_sort) );
        checkCudaErrors( cudaMalloc((void **)&top_index_gpu[gpu_id], sizeof(int)*nBlocks_nf3*num_sort) );

        top_score_host[gpu_id] = new float[nBlocks_nf3];
        top_index_host[gpu_id] = new int[nBlocks_nf3];

    }

    _cputime->record_malloc( sizeof(float)*nBlocks_nf3*num_gpu + sizeof(int)*nBlocks_nf3*num_gpu );

    //*
    size_t devmem_use, devmem_free, devmem_total;
    cudaMemGetInfo(&devmem_free, &devmem_total);
    devmem_use = devmem_total - devmem_free;
    printf("# GPU Memory : Use %3.1f MB (%4.1f%%), Free %3.1f MB (%4.1f%%), Total %3.1f MB\n",(float)devmem_use/1024.0/1024.0,(float)(100*devmem_use/devmem_total), (float)devmem_free/1024.0/1024.0, (float)(100*devmem_free/devmem_total), (float)devmem_total/1024.0/1024.0);
    //*/

#endif /* CUFFT */

    _FFTWin  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nf3*nproc2);
    _FFTWout = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nf3*nproc2);

    plan_fftw_forward = new fftwf_plan[nproc2];
    plan_fftw_inverse = new fftwf_plan[nproc2];

    #pragma omp parallel private(myid2) num_threads(nproc2) //limit num of procs
    {
        myid2 = omp_get_thread_num();
        #pragma omp for
        for(int id = 0; id < nproc2; id++) {
            #pragma omp critical
            {
                plan_fftw_forward[myid2]=fftwf_plan_dft_3d(nf1,nf1,nf1,&_FFTWin[nf3*myid2],&_FFTWout[nf3*myid2],FFTW_FORWARD,FFTW_ESTIMATE);
                plan_fftw_inverse[myid2]=fftwf_plan_dft_3d(nf1,nf1,nf1,&_FFTWin[nf3*myid2],&_FFTWout[nf3*myid2],FFTW_BACKWARD,FFTW_ESTIMATE);
            }
        }
    }

    _cputime->record_malloc( sizeof(fftwf_complex)*nf3*2*nproc2 );
    return;
}

//============================================================================//
void FFTProcess::receptor_fft(float *grid_r,float *grid_i)
//============================================================================//
{
    const int num_grid= _Num_fft / 2;
    const size_t nf3 = _Num_fft * _Num_fft * _Num_fft;
    const int ndata   = ( _Num_fft - num_grid ) / 2;
    const float   theta   = -2.0 * PI / _Num_fft;

    const int num_gpu = _parallel->num_gpu();

    if(num_gpu > 0) {
#ifdef CUFFT
        int myid2;
        struct timeval et1, et2;
        //memset(CUFFTin_host[0], make_cuComplex(0.0, 0.0), sizeof(cufftComplex)*nf3);
        for( int i = 0 ; i < nf3 ; i++ ) {
            CUFFTin_host[i] = make_cuComplex(0.0, 0.0);
        }

        for( int i = 0, m = 0 ; i < num_grid ; i++ ) {
            const int ic = _Num_fft*_Num_fft*(i+ndata);
            for( int j = 0 ; j < num_grid ; j++ ) {
                const int jc = ic + _Num_fft*(j+ndata);
                for( int k = 0 ; k < num_grid ; k++ ) {
                    CUFFTin_host[jc+k+ndata] = make_cuComplex(grid_r[m  ], grid_i[m]);
                    m++;
                }
            }
        }

        cudaSetDevice(0); //CUFFTin_dev[0] : [0] means 0th GPU

        gettimeofday(&et1,NULL);
        checkCudaErrors( cudaMemcpy(CUFFTin_gpu[0], CUFFTin_host, sizeof(cufftComplex)*nf3, cudaMemcpyHostToDevice) );
        gettimeofday(&et2,NULL);
        _cputime->t6_data_transfer_rec += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

        fft3d(theta,0); // [0] means performed on 0th GPU

        gettimeofday(&et1,NULL);
        checkCudaErrors( cudaMemcpy(CUFFTout_host,CUFFTout_gpu[0],sizeof(cufftComplex)*nf3,cudaMemcpyDeviceToHost) );
        gettimeofday(&et2,NULL);
        _cputime->t6_data_transfer_rec += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

        for( int i = 0 ; i < nf3 ; i++ ) {
            _FFT_rec_r[i] = cuCrealf(CUFFTout_host[i]);
            _FFT_rec_i[i] = cuCimagf(CUFFTout_host[i]);
        }

        gettimeofday(&et1,NULL);

        #pragma omp parallel private(myid2) num_threads(num_gpu)
        {
            myid2 = omp_get_thread_num();
            #pragma omp for
            for(int gpu_id = 0; gpu_id < num_gpu; gpu_id++) {
                cudaSetDevice(myid2);
                checkCudaErrors( cudaMemcpy(_FFT_rec_r_gpu[myid2], _FFT_rec_r, sizeof(float)*nf3, cudaMemcpyHostToDevice) );
                checkCudaErrors( cudaMemcpy(_FFT_rec_i_gpu[myid2], _FFT_rec_i, sizeof(float)*nf3, cudaMemcpyHostToDevice) );
            }
        }

        gettimeofday(&et2,NULL);
        _cputime->t6_data_transfer_rec += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
#endif
    } else {
        memset(_FFTWin, 0.0, sizeof(fftwf_complex)*nf3);

        for( int i = 0, m = 0 ; i < num_grid ; i++ ) {
            const int ic = _Num_fft*_Num_fft*(i+ndata);

            for( int j = 0 ; j < num_grid ; j++ ) {
                const int jc = ic + _Num_fft*(j+ndata);

                for( int k = 0 ; k < num_grid ; k++ ) {
                    _FFTWin[jc+k+ndata][0] = grid_r[m  ];
                    _FFTWin[jc+k+ndata][1] = grid_i[m++];
                }
            }
        }

        fft3d(theta,0);

        for( int i = 0 ; i < nf3 ; i++ ) {
            _FFT_rec_r[i] = _FFTWout[i][0];
            _FFT_rec_i[i] = _FFTWout[i][1];
        }
    }


    return;
}

//============================================================================//
void FFTProcess::ligand_preparation(float *grid_r,float *grid_i, size_t myid2)
//============================================================================//
{
    const int ng1 = _Num_fft / 2;
    const int nf2 = _Num_fft * _Num_fft;
    const size_t nf3 = _Num_fft * _Num_fft * _Num_fft;
    const int ndata   = ( _Num_fft - ng1 ) / 2;
   
    memset(_FFTWin[nf3*myid2], 0.0, sizeof(fftwf_complex)*nf3);
        
    for( int i = 0, m = 0 ; i < ng1 ; i++ ) {
        const int ic = nf2*(i+ndata);

        for( int j = 0 ; j < ng1 ; j++ ) {
            int jc = ic + _Num_fft*(j+ndata);
            
            for( size_t k = 0, myijk=nf3*myid2+jc+ndata ; k < ng1 ; k++, myijk++ ) {
                _FFTWin[myijk][0] = grid_r[m  ];
                _FFTWin[myijk][1] = grid_i[m++];
            }
        }
    }
    
    return;
}

//============================================================================//
void FFTProcess::convolution(size_t myid2)
//============================================================================//
{
    const int nf1 = _Num_fft;
    const int nf2 = nf1*nf1;
    const size_t nf3 = nf1*nf2;

    for( size_t i = 0, j=nf3*myid2 ; i < nf3 ; i++,j++ ) {
      _FFTWin[j][0] = _FFT_rec_r[i]*_FFTWout[j][0] + _FFT_rec_i[i]*_FFTWout[j][1];
      _FFTWin[j][1] = _FFT_rec_r[i]*_FFTWout[j][1] - _FFT_rec_i[i]*_FFTWout[j][0];
    }

    return;
}

//============================================================================//
void FFTProcess::fft3d(const float &theta, size_t myid2)
//============================================================================//
{   
    const size_t nproc2  = _parallel->nproc2();
    const int num_gpu = _parallel->num_gpu();
    struct timeval et3, et4;

    if(myid2 < num_gpu) {
#ifdef CUFFT
        const int nf1 = _Num_fft;
        cufftHandle plan;
        cufftResult res;

        res = cufftPlan3d(&plan, nf1, nf1, nf1, CUFFT_C2C);
        if(!res == CUFFT_SUCCESS) {
            cout << "!fail to plan 3d FFT (DFT):" << res << endl;
            exit(-1);
        }

        if( theta < 0.0 ) {
            res = cufftExecC2C(plan, &CUFFTin_gpu[myid2][0], &CUFFTout_gpu[myid2][0], CUFFT_FORWARD);
        } else {
            res = cufftExecC2C(plan, &CUFFTin_gpu[myid2][0], &CUFFTout_gpu[myid2][0], CUFFT_INVERSE);
        }

        if(!res == CUFFT_SUCCESS) {
            cout << "!fail to exec 3d FFT(in fft3d()):" << res << endl;
            exit(-1);
        }

        res =  cufftDestroy(plan);
#endif
    } else {
        gettimeofday(&et3,NULL);
        if( _parameter->fft_library_type == 2 ) {        
	} else {
            if( theta < 0.0 ) {
                fftwf_execute(plan_fftw_forward[myid2]);
            } else {
                fftwf_execute(plan_fftw_inverse[myid2]);
            }
        }
        gettimeofday(&et4,NULL);
        //printf(" [FFT(host),%s] %10.5f\n\n",((theta<0.0)?"Forward":"Inverse"),(et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6)));
    }

    return;
}

//============================================================================//
void FFTProcess::score_sort(size_t myid2)
//============================================================================//
{
    const int num_sort  = _parameter->_Num_sort;
    const int nf2 = _Num_fft * _Num_fft;
    const int nf3 = _Num_fft * _Num_fft * _Num_fft;
    float temp_top_score;
    int temp_top_index;

    for( int i = 0 ; i < num_sort ; i++ ) {
        _Select[myid2][i].score = -99999.0;
    }

    fftwf_complex *fftout;
    fftout = _FFTWout;
    
    if(num_sort!=1) {
        for( size_t i = 0,myi= nf3*myid2 ; i < nf3 ; i++,myi++ ) {
            const float raw = fftout[myi][0] / nf3;
            if( raw < _Select[myid2][num_sort-1].score) continue;
            for( int j = 0 ; j < num_sort ; j++ ) {
                if( raw > _Select[myid2][j].score ) {
                    for( int k = num_sort-1 ; k > j ; k-- ) {
                        _Select[myid2][k] = _Select[myid2][k-1];
                    }
                    _Select[myid2][j].score    = raw;
                    _Select[myid2][j].index[1] = i / nf2;
                    _Select[myid2][j].index[2] = (i / _Num_fft) % _Num_fft;
                    _Select[myid2][j].index[3] = i % _Num_fft;
                    break;
                }
            }
        }
    } else { // num_sort = 1, take only 1 score per angle
        temp_top_score = 0.0;
        temp_top_index = 0;
        for( size_t i = 0, myi=nf3*myid2 ; i < nf3 ; i++,myi++ ) {
            const float raw = fftout[myi][0];
            if (temp_top_score < raw) {
                temp_top_score = raw;
                temp_top_index = i;
            }
        }
        _Select[myid2][0].score    = temp_top_score / nf3;
        _Select[myid2][0].index[1] = temp_top_index / nf2;
        _Select[myid2][0].index[2] = (temp_top_index / _Num_fft) % _Num_fft;
        _Select[myid2][0].index[3] = temp_top_index % _Num_fft;
    }

    for( int i = 0 ; i < num_sort ; i++ ) {
        //printf(" top %d %f\n",i,_Select[myid2][i].score);
        _Select[myid2][i].index[0] = _Current_rot_angle_num[myid2];
    }

    for( int i = 0 ; i < num_sort ; i++ ) {
        _Top[_Current_rot_angle_num[myid2]*num_sort+i] = _Select[myid2][i];
    }

    return;
}


//============================================================================//
bool ascending_struct(const SortScore &s1,const SortScore &s2)  // for sort
//============================================================================//
{
    return s1.score > s2.score;
}

//============================================================================//
void FFTProcess::sort_index(float *fwork,int *iwork)
//============================================================================//
{
    const int     no  = _parameter->_Num_output;
    const int     nt  = (_Top.size() < no) ? _Top.size() : no;

    partial_sort(_Top.begin(),_Top.begin()+nt,_Top.end(),
                 ascending_struct);

    return;
}

//============================================================================//
void FFTProcess::top_score_clean()
//============================================================================//
{
    const int     num_sort  = _parameter->_Num_sort;
    const int     num_angle  = _parameter->_Num_rot_angles;
    const int     no  = _parameter->_Num_output;
    int   num_toprank;

    num_toprank = num_angle * num_sort;
    if( no > num_toprank ) num_toprank = no;

    for( int j = 0 ; j < num_toprank ; j++ ) {
        _Top[j].score = 0.0;
    }

    return;
}

#ifdef CUFFT
//============================================================================//
void FFTProcess::cuda_fft(float *grid_r,float *grid_i,float *grid_coord,float *atom_coord_rotated,float *theta, size_t myid2)
//============================================================================//
{
    const int nf1 = _Num_fft;
    const int nf2 = nf1 * nf1;
    const size_t nf3 = nf2 * nf1;

    const int num_sort = _parameter->_Num_sort;
    const int na = _ligand->num_atoms();

    struct timeval et1, et2;
    struct timeval et3, et4;
    if(myid2==0) gettimeofday(&et1,NULL);

    float temp_top_score = -999999.0;
    int temp_top_index = -999999;

    const int nThreads = NUM_THREADS;
    const int nBlocks_nf3 = (nf3 + (nThreads-1)) / nThreads;
    if(nBlocks_nf3 * nThreads < nf3) {
        printf(" nf3:%d, nBlocks_nf3:%d, nThreads:%d , nf3=nBlocks_nf3*nThreads\n",nf3,nBlocks_nf3,nThreads);
        fprintf(stderr, " [ERROR] too large FFT size. nf3:%d, nBlocks_nf3:%d\n", nf3, nBlocks_nf3);
        exit(1);
    }

    cudaSetDevice(myid2);
    //printf(" #p10 [myid=%d]\n",myid2);

    ligand_voxelization_on_gpu(theta,myid2);
    checkCudaErrors( cudaDeviceSynchronize() );

    if(myid2==0) gettimeofday(&et2,NULL);
    if(myid2==0) _cputime->t3_1_ligand_voxelization += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
    if(myid2==0) gettimeofday(&et1,NULL);

    cufft_result[myid2] = cufftExecC2C(cufft_plan[myid2], &CUFFTin_gpu[myid2][0], &CUFFTout_gpu[myid2][0], CUFFT_FORWARD);
    if(!cufft_result[myid2] == CUFFT_SUCCESS) {
        cout << "!fail to exec 3d FFT (DFT, Lig):" << cufft_result[myid2] << endl;
        exit(-1);
    }

    //*/
    checkCudaErrors( cudaDeviceSynchronize() );

    if(myid2==0) gettimeofday(&et2,NULL);
    if(myid2==0) _cputime->t3_2_fftprocess_ligand_fft += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    if(myid2==0) gettimeofday(&et1,NULL);
    convolution_gpu<<<nBlocks_nf3, nThreads>>>(nf3, _FFT_rec_r_gpu[myid2], _FFT_rec_i_gpu[myid2], CUFFTout_gpu[myid2], CUFFTin_gpu[myid2]);

    checkCudaErrors( cudaDeviceSynchronize() );

    if(myid2==0) gettimeofday(&et2,NULL);
    if(myid2==0) _cputime->t3_3_fftprocess_convolution += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
    if(myid2==0) gettimeofday(&et1,NULL);

    cufft_result[myid2] = cufftExecC2C(cufft_plan[myid2], &CUFFTin_gpu[myid2][0], &CUFFTout_gpu[myid2][0], CUFFT_INVERSE);
    if(!(cufft_result[myid2] == CUFFT_SUCCESS)) {
        cout << "!fail to exec 3d FFT (IDFT):" << cufft_result[myid2] << endl;
        exit(-1);
    }
    //*
    checkCudaErrors( cudaDeviceSynchronize() );
    if(myid2==0) gettimeofday(&et2,NULL);
    if(myid2==0) _cputime->t3_4_fftprocess_fft_inverse += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
    if(myid2==0) gettimeofday(&et1,NULL);

    // Search max score translation position from CUFFTout_gpu[nf3]

    //printf(" t=%d per angle\n",num_sort);

    for( int i = 0 ; i < num_sort ; i++ ) {
        _Select[myid2][i].score = -99999.0;
    }

    max_pos_single<<<nBlocks_nf3, nThreads, sizeof(float)*nThreads>>>(nf3, CUFFTout_gpu[myid2],  top_score_gpu[myid2], top_index_gpu[myid2]);
    checkCudaErrors( cudaDeviceSynchronize() );

    if(myid2==0) gettimeofday(&et3,NULL);
    checkCudaErrors( cudaMemcpy(top_score_host[myid2],top_score_gpu[myid2],sizeof(float)*nBlocks_nf3,cudaMemcpyDeviceToHost) );
    checkCudaErrors( cudaMemcpy(top_index_host[myid2],top_index_gpu[myid2],sizeof(int)*nBlocks_nf3,cudaMemcpyDeviceToHost) );
    if(myid2==0) gettimeofday(&et4,NULL);
    if(myid2==0) _cputime->t6_data_transfer_in_loop += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));
    checkCudaErrors( cudaDeviceSynchronize() );

    if(num_sort!=1) {
        for(int i=0; i<nBlocks_nf3; i++) {
            if(top_index_host[myid2][i]/nf2 > nf1 || top_index_host[myid2][i] < 0){
                top_score_host[myid2][i] = -99999.99;
                //printf(" error, %d | score, %f \n", top_index_host[myid2][i]/nf2, top_score_host[myid2][i]);
            }
            const float raw = top_score_host[myid2][i];
            if( raw < _Select[myid2][num_sort-1].score) continue;
            for( int j = 0 ; j < num_sort ; j++ ) {
                if( raw > _Select[myid2][j].score ) {
                    for( int k = num_sort-1 ; k > j ; k-- ) {
                        _Select[myid2][k] = _Select[myid2][k-1];
                    }
                    const int index = top_index_host[myid2][i];
                    _Select[myid2][j].score    = raw;
                    _Select[myid2][j].index[1] = index / nf2;
                    _Select[myid2][j].index[2] = (index / _Num_fft) % _Num_fft;
                    _Select[myid2][j].index[3] = index % _Num_fft;
                    break;
                }
            }
        }

    } else { // num_sort = 1, select only 1 score per 1 ligand angle
        for(int i=0; i<nBlocks_nf3; i++) {
            if(top_index_host[myid2][i]/nf2 > nf1 || top_index_host[myid2][i] < 0){
                top_score_host[myid2][i] = -99999.99;
                //printf(" error, %d | score, %f \n", top_index_host[myid2][i]/nf2, top_score_host[myid2][i]);
            }
            if(temp_top_score < top_score_host[myid2][i]) {
                temp_top_score = top_score_host[myid2][i];
                temp_top_index = top_index_host[myid2][i];
            }
        }

        //printf("  m:%f\n\n",temp_top_score);
        //printf("%g (%d) [%d %d %d]\n", temp_top_score, _p, temp_top_index/(n*n),(temp_top_index/n)%n, temp_top_index%n );
        //printf("<%d> %g (%d/%d) %d\n", nBlocks,temp_top_score, temp_top_index, nf3, temp_top_index/nf2);

        _Select[myid2][0].score    = temp_top_score;
        _Select[myid2][0].index[1] = temp_top_index / nf2;
        _Select[myid2][0].index[2] = (temp_top_index / nf1) % nf1;
        _Select[myid2][0].index[3] = temp_top_index % nf1;
        /* / DEBUG
        printf("TEST,  %d\n", _Select[myid2][0].index[1]);
        if ( _Select[myid2][0].index[1] > nf1 ){
            printf(" error, %d\n", _Select[myid2][0].index[1]);
            }*/

    }

    //*** score_sort ***********************************************************

    for( int i = 0 ; i < num_sort ; i++ ) {
        _Select[myid2][i].index[0] = _Current_rot_angle_num[myid2];
        _Top[_Current_rot_angle_num[myid2]*num_sort+i] = _Select[myid2][i];
    }

    //size_t devmem_use, devmem_free, devmem_total;
    //cudaMemGetInfo(&devmem_free, &devmem_total);
    //devmem_use = devmem_total - devmem_free;
    //printf(" [GPU (%d) memory] Use : %10u (%4.1f%%), Free : %10u (%4.1f%%), Total : %10u\n",myid2,devmem_use,(float)(100*devmem_use/devmem_total), devmem_free, (float)(100*devmem_free/devmem_total), devmem_total);


    checkCudaErrors( cudaDeviceSynchronize() );
    if(myid2==0) gettimeofday(&et2,NULL);
    if(myid2==0) _cputime->t3_5_fftprocess_score_sort += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    return;
}


//============================================================================//
void FFTProcess::ligand_voxelization_on_gpu(float *theta, size_t myid2)
//============================================================================//
{
    const int ng1 = _Num_fft / 2;
    const int ng3 = ng1 * ng1 * ng1;
    const int nf1 = _Num_fft;
    const int nf2 = nf1 * nf1;
    const size_t nf3 = nf2 * nf1;

    const float delta = 1.0;
    const float surface = 1.0;
    const float grid_width = _parameter->grid_width;

    const int na = _ligand->num_atoms();
    const int nag = na * ng1;

    struct timeval et1, et2;
    struct timeval et3, et4;

    const int nThreads = NUM_THREADS;
    const int nBlocks_na = (na + (nThreads-1)) / nThreads;
    const int nBlocks_nag = (nag + (nThreads-1)) / nThreads;
    const int nBlocks_ng3 = (ng3 + (nThreads-1)) / nThreads;
    const int nBlocks_nf3 = (nf3 + (nThreads-1)) / nThreads;
    if(nBlocks_nf3 * nThreads < nf3) {
        printf(" nf3:%d, nBlocks_nf3:%d, nThreads:%d , nf3=nBlocks_nf3*nThreads\n",nf3,nBlocks_nf3,nThreads);
        fprintf(stderr, " [ERROR] too large FFT size. nf3:%d, nBlocks_nf3:%d\n", nf3, nBlocks_nf3);
        exit(1);
    }

    //*
    //transfer ligand angle & calc xd,yd,zd,atom_coord_rotated
    if(myid2==0) gettimeofday(&et3,NULL);

    if(myid2==0) gettimeofday(&et1,NULL);
    checkCudaErrors( cudaMemcpy(ligand_rotation_angle_gpu[myid2], theta, sizeof(float)*3, cudaMemcpyHostToDevice) );
    if(myid2==0) gettimeofday(&et2,NULL);
    if(myid2==0) _cputime->t3_1_ligand_voxelization += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    lig_rotation<<<nBlocks_na, nThreads>>>(na, ligand_rotation_angle_gpu[myid2],atom_coord_orig_gpu[myid2], mole_center_coord_gpu[myid2], atom_coord_rotated_gpu[myid2]);
    lig_calc_dis_atomgrid<<<nBlocks_nag, nThreads>>>(na, ng1, xd_gpu[myid2], yd_gpu[myid2], zd_gpu[myid2], grid_coord_gpu[myid2], atom_coord_rotated_gpu[myid2]);
    checkCudaErrors( cudaDeviceSynchronize() );
    if(myid2==0) gettimeofday(&et4,NULL);
    if(myid2==0) _cputime->t3_1_1_ligvoxgpu_copy_htod += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //grid[] initialize
    if(myid2==0) gettimeofday(&et3,NULL);
    lig_vox_init_grid<<<nBlocks_ng3, nThreads>>>(ng3,grid_r_gpu[myid2],grid_i_gpu[myid2]);
    lig_vox_init_fft<<<nBlocks_nf3, nThreads>>>(nf3,CUFFTin_gpu[myid2]);
    checkCudaErrors( cudaDeviceSynchronize() );
    if(myid2==0) gettimeofday(&et4,NULL);
    if(myid2==0) _cputime->t3_1_2_ligvoxgpu_kernel_init += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //atom fill(core)
    if(myid2==0) gettimeofday(&et3,NULL);
    lig_vox_fill<<<nBlocks_ng3, nThreads>>>
    (ng1,na,delta,radius_core2_gpu[myid2],xd_gpu[myid2],yd_gpu[myid2],zd_gpu[myid2],grid_coord_gpu[myid2],atom_coord_rotated_gpu[myid2],grid_r_gpu[myid2], grid_width);
    checkCudaErrors( cudaDeviceSynchronize() );
    if(myid2==0) gettimeofday(&et4,NULL);
    if(myid2==0) _cputime->t3_1_3_ligvoxgpu_kernel_fill_core += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //surface cutting
    if(myid2==0) gettimeofday(&et3,NULL);
    lig_vox_surface_cut_CtoT<<<nBlocks_ng3, nThreads>>>(ng1,delta,grid_r_gpu[myid2]);
    checkCudaErrors( cudaDeviceSynchronize() );
    if(myid2==0) gettimeofday(&et4,NULL);
    if(myid2==0) _cputime->t3_1_4_ligvoxgpu_kernel_cut_surf += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //atom fill(surf)
    if(myid2==0) gettimeofday(&et3,NULL);
    lig_vox_fill<<<nBlocks_ng3, nThreads>>>
    (ng1,na,surface,radius_surf2_gpu[myid2],xd_gpu[myid2],yd_gpu[myid2],zd_gpu[myid2],grid_coord_gpu[myid2],atom_coord_rotated_gpu[myid2],grid_r_gpu[myid2], grid_width);
    checkCudaErrors( cudaDeviceSynchronize() );
    if(myid2==0) gettimeofday(&et4,NULL);
    if(myid2==0) _cputime->t3_1_5_ligvoxgpu_kernel_fill_surf += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //electro
    if(myid2==0) gettimeofday(&et3,NULL);

    if(_parameter->lig_elec_serial_flag == 0) {
        lig_vox_elec<<<nBlocks_ng3, nThreads>>>(ng1, na, grid_width, _Charge_gpu[myid2], atom_coord_rotated_gpu[myid2], grid_i_gpu[myid2]);
    } else {
        lig_vox_elec_serial<<<nBlocks_ng3, nThreads>>>(ng1, na, grid_width, _Charge_gpu[myid2], atom_coord_rotated_gpu[myid2], grid_i_gpu[myid2]);
    }

    /*
    float *tem_grid;
    const int ng2=ng1*ng1;
    tem_grid = new float[ng3];
    checkCudaErrors( cudaMemcpy(tem_grid, grid_i_gpu[myid2], sizeof(float)*ng3, cudaMemcpyDeviceToHost) );
    //for(int i=0;i<ng3;i++) if(tem_grid[i]!=0.0) printf(" [%03d,%03d,%03d] :  %6.3f\n",i/ng2,(i/ng1)%ng1,i%ng1,tem_grid[i]);
    //*/

    checkCudaErrors( cudaDeviceSynchronize() );
    if(myid2==0) gettimeofday(&et4,NULL);
    if(myid2==0) _cputime->t3_1_6_ligvoxgpu_kernel_elec += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //set Voxel grid[ng3] into center of FFT grid[nf3]
    if(myid2==0) gettimeofday(&et3,NULL);
    ligand_voxel_set<<<nBlocks_ng3, nThreads>>>(ng1,CUFFTin_gpu[myid2],grid_r_gpu[myid2],grid_i_gpu[myid2]);
    checkCudaErrors( cudaDeviceSynchronize() );
    if(myid2==0) gettimeofday(&et4,NULL);
    if(myid2==0) _cputime->t3_1_7_ligvoxgpu_kernel_set_array += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    //*/
}


//============================================================================//
void FFTProcess::ligand_data_transfer_gpu(float *grid_coord)
//============================================================================//
{
    const int ng1 = _Num_fft / 2;
    const int na = _ligand->num_atoms();
    const int num_gpu = _parallel->num_gpu();
    const float   rcore2 = 1.5;           // ZDOCK parameter
    const float   rsurf2 = 1.0;           // ZDOCK parameter
    struct timeval et1, et2;

    float *radius_core2;
    float *radius_surf2;
    radius_core2 = new float[na];
    radius_surf2 = new float[na];

    for(int i = 0; i < na; i++) {
        radius_core2[i] = _ligand->_Radius[i] * _ligand->_Radius[i] * rcore2;
        radius_surf2[i] = _ligand->_Radius[i] * _ligand->_Radius[i] * rsurf2;
    }

    gettimeofday(&et1,NULL);
    for(int gpu_id = 0; gpu_id < num_gpu; gpu_id++) {
        cudaSetDevice(gpu_id);
        checkCudaErrors( cudaMemcpy(radius_core2_gpu[gpu_id], radius_core2, sizeof(float)*na, cudaMemcpyHostToDevice) );
        checkCudaErrors( cudaMemcpy(radius_surf2_gpu[gpu_id], radius_surf2, sizeof(float)*na, cudaMemcpyHostToDevice) );
        checkCudaErrors( cudaMemcpy(_Charge_gpu[gpu_id], _ligand->_Charge, sizeof(float)*na, cudaMemcpyHostToDevice) );
        checkCudaErrors( cudaMemcpy(grid_coord_gpu[gpu_id], grid_coord, sizeof(float)*ng1, cudaMemcpyHostToDevice) );
        checkCudaErrors( cudaMemcpy(atom_coord_orig_gpu[gpu_id], _ligand->_Coordinate, sizeof(float)*na*3, cudaMemcpyHostToDevice) );
        checkCudaErrors( cudaMemcpy(mole_center_coord_gpu[gpu_id], _ligand->_Center, sizeof(float)*3, cudaMemcpyHostToDevice) );
    }

    gettimeofday(&et2,NULL);
    _cputime->t6_data_transfer_lig += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    delete[] radius_core2;
    delete[] radius_surf2;

    return;
}

#endif /* CUFFT */



//============================================================================//
void FFTProcess::fft_memory_free()
//============================================================================//
{
    const size_t nproc2    = _parallel->nproc2();
    const int num_gpu = _parallel->num_gpu();
    const size_t nf3 = _Num_fft * _Num_fft * _Num_fft;

    for(int id = 0; id < nproc2; id++) {
        fftwf_destroy_plan(plan_fftw_forward[id]);
        fftwf_destroy_plan(plan_fftw_inverse[id]);
    }

    _cputime->record_free(sizeof(float)*nf3*2*(nproc2));

#ifdef CUFFT

    //const int num_sort = _parameter->_Num_sort;
    const int nThreads = NUM_THREADS;
    const int nBlocks_nf3 = (nf3 + (nThreads-1)) / nThreads;

    for(int gpu_id = 0; gpu_id < num_gpu; gpu_id++) {
        cudaSetDevice(gpu_id);

        cufftDestroy(cufft_plan[gpu_id]);

        checkCudaErrors( cudaFree(CUFFTin_gpu[gpu_id]));
        checkCudaErrors( cudaFree(CUFFTout_gpu[gpu_id]));
        checkCudaErrors( cudaFree(_FFT_rec_r_gpu[gpu_id]));
        checkCudaErrors( cudaFree(_FFT_rec_i_gpu[gpu_id]));

        checkCudaErrors( cudaFree(grid_r_gpu[gpu_id]));
        checkCudaErrors( cudaFree(grid_i_gpu[gpu_id]));
        checkCudaErrors( cudaFree(grid_coord_gpu[gpu_id]));

        checkCudaErrors( cudaFree(radius_core2_gpu[gpu_id]));
        checkCudaErrors( cudaFree(radius_surf2_gpu[gpu_id]));
        checkCudaErrors( cudaFree(_Charge_gpu[gpu_id]));

        checkCudaErrors( cudaFree(xd_gpu[gpu_id]));
        checkCudaErrors( cudaFree(yd_gpu[gpu_id]));

        checkCudaErrors( cudaFree(zd_gpu[gpu_id]));

        checkCudaErrors( cudaFree(atom_coord_rotated_gpu[gpu_id]));
        checkCudaErrors( cudaFree(atom_coord_orig_gpu[gpu_id]));
        checkCudaErrors( cudaFree(mole_center_coord_gpu[gpu_id]));
        checkCudaErrors( cudaFree(ligand_rotation_angle_gpu[gpu_id]));

        checkCudaErrors( cudaFree(top_score_gpu[gpu_id]));
        checkCudaErrors( cudaFree(top_index_gpu[gpu_id]));

        delete [] top_score_host[gpu_id];
        delete [] top_index_host[gpu_id];

    }

    _cputime->record_free( sizeof(float)*nBlocks_nf3*num_gpu + sizeof(int)*nBlocks_nf3*num_gpu );

#endif

    return;
}
