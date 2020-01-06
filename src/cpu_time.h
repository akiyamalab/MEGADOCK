/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : CPUTime
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#ifndef Cpu_time_h
#define Cpu_time_h 1

#include <stdio.h>
#include <string>
#include <iostream>

using namespace std;

class CPUTime
{
private:
    friend class main;
    friend class Control;
    friend class Docking;
    friend class FFTProcess;
    CPUTime(CPUTime &c) {}
    const CPUTime & operator=(const CPUTime &c);

public:
    float t1_initialize;
    float t2_receptor_process;
    float t3_docking_total;
    float t4_docking_output_detail;
    float t5_docking_output;

    float t3_1_ligand_voxelization;
    float t3_2_fftprocess_ligand_fft;
    float t3_3_fftprocess_convolution;
    float t3_4_fftprocess_fft_inverse;
    float t3_5_fftprocess_score_sort;
    float t3_6_fftprocess_sort_index;
    
	float t3_1_1_ligvoxgpu_copy_htod;
	float t3_1_2_ligvoxgpu_kernel_init;
	float t3_1_3_ligvoxgpu_kernel_fill_core;
	float t3_1_4_ligvoxgpu_kernel_cut_surf;
	float t3_1_5_ligvoxgpu_kernel_fill_surf;
	float t3_1_6_ligvoxgpu_kernel_elec;
	float t3_1_7_ligvoxgpu_kernel_set_array;
	
	float t6_data_transfer_rec;
	float t6_data_transfer_lig;
	float t6_data_transfer_in_loop;
	
	float t7_offload_transfer;
	float t7_offload_calc;

	long long m1_current_malloc_size;
	long long m2_maximum_malloc_size;

    CPUTime() {
#ifdef DEBUG
        cout << "Constructing CPUTime.\n";
#endif
    }
    virtual   ~CPUTime()      {
#ifdef DEBUG
        cout << "Destructing CPUTime\n";
#endif
    }
    virtual void      initialize();
    virtual void      output();
    virtual void      record_malloc(const int &size) {
        m1_current_malloc_size += (long long)size;
        if(m1_current_malloc_size > m2_maximum_malloc_size){
        	m2_maximum_malloc_size = m1_current_malloc_size;
        }
        //printf(" # %11d bytes allocated, Current alloc size : %11lld [Bytes], Max alloc size : %11lld [Bytes] (%.2f GB)\n"
        //, size, m1_current_malloc_size, m2_maximum_malloc_size, (float)(m2_maximum_malloc_size/1024.0/1024.0/1024.0));
    }
    virtual void      record_free(const int &size) {
        m1_current_malloc_size -= (long long)size;
        if(m1_current_malloc_size < 0){
        	//printf(" # Allocate memory size error\n");
        }
        //printf(" # %11d bytes freed    , Current alloc size : %11lld [Bytes], Max alloc size : %11lld [Bytes] (%.2f GB)\n"
        //, size, m1_current_malloc_size, m2_maximum_malloc_size, (float)(m2_maximum_malloc_size/1024.0/1024.0/1024.0));
    }
};

#endif
