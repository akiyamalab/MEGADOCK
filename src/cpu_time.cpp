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
//  Class Name : CPUTime
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include "cpu_time.h"

//============================================================================//
void CPUTime::initialize()
//============================================================================//
{

    t1_initialize = 0.0;
    t2_receptor_process = 0.0;
    t3_docking_total = 0.0;
    t4_docking_output_detail = 0.0;
    t5_docking_output = 0.0;

    t3_1_ligand_voxelization = 0.0;
    t3_2_fftprocess_ligand_fft = 0.0;
    t3_3_fftprocess_convolution = 0.0;
    t3_4_fftprocess_fft_inverse = 0.0;
    t3_5_fftprocess_score_sort = 0.0;
    t3_6_fftprocess_sort_index = 0.0;

    t3_1_1_ligvoxgpu_copy_htod = 0.0;
    t3_1_2_ligvoxgpu_kernel_init = 0.0;
    t3_1_3_ligvoxgpu_kernel_fill_core = 0.0;
    t3_1_4_ligvoxgpu_kernel_cut_surf = 0.0;
    t3_1_5_ligvoxgpu_kernel_fill_surf = 0.0;
    t3_1_6_ligvoxgpu_kernel_elec = 0.0;
    t3_1_7_ligvoxgpu_kernel_set_array = 0.0;
    
    t6_data_transfer_rec = 0.0;
    t6_data_transfer_lig = 0.0;
    t6_data_transfer_in_loop = 0.0;
	
    t7_offload_transfer = 0.0;
    t7_offload_calc = 0.0;

    m1_current_malloc_size = 0;
    m2_maximum_malloc_size = 0;

    return;
}

//============================================================================//
void CPUTime::output()
//============================================================================//
{

    float total_time =  t1_initialize+
                        t2_receptor_process+
                        t3_docking_total+
                        t4_docking_output_detail+
                        t5_docking_output;

    float total_transfer =  t6_data_transfer_rec+
                            t6_data_transfer_lig+
                            t6_data_transfer_in_loop;

    //correct the calc-time gap due to parallel calc
    const float cufft_ligvox_on_GPU_total_idzero = 0.00000000001
    	                                          + t3_1_1_ligvoxgpu_copy_htod
                                                  + t3_1_2_ligvoxgpu_kernel_init
                                                  + t3_1_3_ligvoxgpu_kernel_fill_core
                                                  + t3_1_4_ligvoxgpu_kernel_cut_surf
                                                  + t3_1_5_ligvoxgpu_kernel_fill_surf
                                                  + t3_1_6_ligvoxgpu_kernel_elec
                                                  + t3_1_7_ligvoxgpu_kernel_set_array;

    //printf(" max:%f, one:%f",t3_1_ligand_voxelization,cufft_ligvox_on_GPU_total_idzero);
    t3_1_1_ligvoxgpu_copy_htod            *= (t3_1_ligand_voxelization / cufft_ligvox_on_GPU_total_idzero);
    t3_1_2_ligvoxgpu_kernel_init          *= (t3_1_ligand_voxelization / cufft_ligvox_on_GPU_total_idzero);
    t3_1_3_ligvoxgpu_kernel_fill_core     *= (t3_1_ligand_voxelization / cufft_ligvox_on_GPU_total_idzero);
    t3_1_4_ligvoxgpu_kernel_cut_surf      *= (t3_1_ligand_voxelization / cufft_ligvox_on_GPU_total_idzero);
    t3_1_5_ligvoxgpu_kernel_fill_surf     *= (t3_1_ligand_voxelization / cufft_ligvox_on_GPU_total_idzero);
    t3_1_6_ligvoxgpu_kernel_elec          *= (t3_1_ligand_voxelization / cufft_ligvox_on_GPU_total_idzero);
    t3_1_7_ligvoxgpu_kernel_set_array     *= (t3_1_ligand_voxelization / cufft_ligvox_on_GPU_total_idzero);
    
    float t3_docking_total_idzero = 0.00000000001
                                       + t3_1_ligand_voxelization
                                       + t3_2_fftprocess_ligand_fft
                                       + t3_3_fftprocess_convolution
                                       + t3_4_fftprocess_fft_inverse
                                       + t3_5_fftprocess_score_sort
                                       + t3_6_fftprocess_sort_index;
    //printf(" max:%f, one:%f",t3_docking_total,t3_docking_total_idzero);

    printf("\n");
    printf("CPU time (initialize)         = %8.2f sec. (%4.1f%%)\n",t1_initialize                     , 100.0*t1_initialize                     /total_time);
    printf("CPU time (receptor process)   = %8.2f sec. (%4.1f%%)\n",t2_receptor_process               , 100.0*t2_receptor_process               /total_time);
    printf("CPU time (docking)            = %8.2f sec. (%4.1f%%)\n",t3_docking_total                  , 100.0*t3_docking_total                  /total_time);
    printf("  | Ligand voxelization       = %8.2f sec. (%4.1f%%)\n",t3_1_ligand_voxelization          , 100.0*t3_1_ligand_voxelization          /total_time);
#ifdef CUFFT
if(cufft_ligvox_on_GPU_total_idzero>0.1){
    printf("    | ligvox_copy_htod        = %8.2f sec. (%4.1f%%)\n",t3_1_1_ligvoxgpu_copy_htod        , 100.0*t3_1_1_ligvoxgpu_copy_htod        /total_time);
    printf("    | ligvox_kernel_init      = %8.2f sec. (%4.1f%%)\n",t3_1_2_ligvoxgpu_kernel_init      , 100.0*t3_1_2_ligvoxgpu_kernel_init      /total_time);
    printf("    | ligvox_kernel_fill_core = %8.2f sec. (%4.1f%%)\n",t3_1_3_ligvoxgpu_kernel_fill_core , 100.0*t3_1_3_ligvoxgpu_kernel_fill_core /total_time);
    printf("    | ligvox_kernel_cut_surf  = %8.2f sec. (%4.1f%%)\n",t3_1_4_ligvoxgpu_kernel_cut_surf  , 100.0*t3_1_4_ligvoxgpu_kernel_cut_surf  /total_time);
    printf("    | ligvox_kernel_fill_surf = %8.2f sec. (%4.1f%%)\n",t3_1_5_ligvoxgpu_kernel_fill_surf , 100.0*t3_1_5_ligvoxgpu_kernel_fill_surf /total_time);
    printf("    | ligvox_kernel_elec      = %8.2f sec. (%4.1f%%)\n",t3_1_6_ligvoxgpu_kernel_elec      , 100.0*t3_1_6_ligvoxgpu_kernel_elec      /total_time);
    printf("    | ligvox_kernel_set_array = %8.2f sec. (%4.1f%%)\n",t3_1_7_ligvoxgpu_kernel_set_array , 100.0*t3_1_7_ligvoxgpu_kernel_set_array /total_time);
}
#endif
    printf("  | Ligand FFT                = %8.2f sec. (%4.1f%%)\n",t3_2_fftprocess_ligand_fft        , 100.0*t3_2_fftprocess_ligand_fft        /total_time);
    printf("  | Convolution               = %8.2f sec. (%4.1f%%)\n",t3_3_fftprocess_convolution       , 100.0*t3_3_fftprocess_convolution       /total_time);
    printf("  | Inverse FFT               = %8.2f sec. (%4.1f%%)\n",t3_4_fftprocess_fft_inverse       , 100.0*t3_4_fftprocess_fft_inverse       /total_time);
    printf("  | Sort score (each core)    = %8.2f sec. (%4.1f%%)\n",t3_5_fftprocess_score_sort        , 100.0*t3_5_fftprocess_score_sort        /total_time);
    printf("  | Sort score (merge)        = %8.2f sec. (%4.1f%%)\n",t3_6_fftprocess_sort_index        , 100.0*t3_6_fftprocess_sort_index        /total_time);
#ifdef PHI
    printf("  | Data transfer             = %8.2f sec. (%4.1f%%)\n",t7_offload_transfer        , 100.0*t7_offload_transfer        /total_time);
#endif
    printf("CPU time (detailed output)    = %8.2f sec. (%4.1f%%)\n",t4_docking_output_detail          , 100.0*t4_docking_output_detail          /total_time);
    printf("CPU time (output)             = %8.2f sec. (%4.1f%%)\n",t5_docking_output                 , 100.0*t5_docking_output                 /total_time);

    return;
}
