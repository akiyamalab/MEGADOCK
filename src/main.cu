/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : (main)
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include <string.h>
#include "cpu_time.h"
#include "control.h"

#ifdef CUFFT
#include <helper_cuda.h>
#define VERSION "4.1.3 for GPU & "
#else
#define VERSION "4.1.3 for CPU & "
#endif

#ifdef MPI_DP
#define VTEXT "multiple nodes"
#else
#define VTEXT "single node"
#endif

#define LASTUPDATED "26 March, 2019"

//============================================================================//
#ifdef MPI_DP
int application(int argc,char *argv[])
#else
int main(int argc, char *argv[])
#endif
//============================================================================//
{
    Parallel  *_parallel;
    CPUTime   *_cputime;
    Control   *_control;

    struct timeval et1, et2;
    struct timeval et3, et4;
    int nproc2 = 0;
    int device_count_gpu = 0;

    gettimeofday(&et1,NULL);
    gettimeofday(&et3,NULL);


    cout << " MEGADOCK ver. "<< VERSION << VTEXT <<  endl;
    cout << "      megadock@bi.c.titech.ac.jp   lastupdated: " << LASTUPDATED << endl;
    cout << endl;

    _cputime = new CPUTime();
    _cputime->initialize();

#ifdef _OPENMP
    #pragma omp parallel
    {
        nproc2 = omp_get_num_threads();
        if(omp_get_thread_num() == 0) {
            cout << "# Using OpenMP parallelization: " << nproc2 << " threads." << endl;
        }
    }
    //printf("#OpenMP version %d\n", _OPENMP);
#else
    nproc2 = 1;
#endif //#ifdef _OPENMP

#ifdef CUFFT
    int nogpu_flag = 0;
    for (int num = 0; num < (argc-1); ++num) {
        if(!strcmp(argv[num], "-G")) {
            if(argv[num+1] != NULL) {
                if(atoi(argv[num+1]) == 0) {
                    nogpu_flag = 1;
                }
            }
        }
    }

    if(nogpu_flag != 1) {
        checkCudaErrors( cudaGetDeviceCount(&device_count_gpu) );
        if (device_count_gpu == 0) {
            fprintf(stderr, "GPU Error: no devices supporting CUDA.\n");
            exit(-1);
        }

        cudaDeviceProp deviceProp;
        checkCudaErrors( cudaGetDeviceProperties(&deviceProp, 0));
        if (deviceProp.major < 1) {
            fprintf(stderr, "GPU Error: device does not support CUDA.\n");
            exit(-1);
        }

        cudaSetDeviceFlags(cudaDeviceMapHost);
        fprintf(stdout, "# Using CUDA device %d: %s\n", 0, deviceProp.name);
        cudaSetDevice(0);
        //fprintf(stdout, "# Init CUDA device OK.\n");

        int cufft_version;
        cufftGetVersion(&cufft_version);
        printf("# CUFFT version : %d\n", cufft_version);
    }
#endif

    _parallel = new Parallel(nproc2);
    _parallel->num_gpu(device_count_gpu); 

#ifdef CUFFT
    printf("# Number of available [threads / GPUs] : [%d / %d]\n",nproc2,device_count_gpu);
#endif

    gettimeofday(&et4,NULL);
    _cputime->t1_initialize += (et4.tv_sec-et3.tv_sec + (float)((et4.tv_usec-et3.tv_usec)*1e-6));

    _control = new Control(_cputime,_parallel);
    _control->initialize(argc,argv);
    _control->execute();

    delete _control;
    delete _parallel;

    _cputime->output();

    delete _cputime;

    gettimeofday(&et2,NULL);

    const float elapsed_time = (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
    printf("\n");
    printf("Elapsed time                  = %8.2f sec.\n",elapsed_time);

    return 0;
}
