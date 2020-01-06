/*
 * Copyright (C) 2020 Tokyo Institute of Technology
 */

//============================================================================//
//
//  Software Name : MEGADOCK
//
//  Class Name : Control
//
//  Contact address : Tokyo Institute of Technology, AKIYAMA Lab.
//
//============================================================================//

#include "control.h"

//============================================================================//
void Control::initialize(int argc,char *argv[])
//============================================================================//
{
    int ngrid;
    vector<int> ngrid_table;

    struct timeval et1, et2;
    gettimeofday(&et1,NULL);

    // Parameter
    _parameter = new Parameter(_parallel);
    _parameter->initialize(argc,argv);
    _cputime->record_malloc( sizeof(float)*_parameter->_Num_rot_angles*3 + sizeof(map<string,float>)*(_parameter->_Charmmr.size() + _parameter->_Charmmc.size() + _parameter->_ACE.size()) ); //Rotation angles[], Atom radius, charge, ACE[]

    // Number of processors limitation
    const int thread_limit = _parameter->_Num_thread_limit;
    const int gpu_limit = _parameter->_Num_GPU_limit;

    if(_parallel->nproc2() > thread_limit) {
        _parallel->nproc2(thread_limit);
    }

    if(_parallel->num_gpu() > gpu_limit || _parallel->num_gpu() > _parallel->nproc2()) {
        _parallel->num_gpu( min(gpu_limit, (int)_parallel->nproc2()) );
    }
    printf("# Using %3d CPU cores, %d GPUs\n", _parallel->nproc2(), _parallel->num_gpu());

    // Receptor
    _receptor = new Receptor(_parameter->_RecPDB_file);
    _receptor->initialize(_parameter);
    _cputime->record_malloc( sizeof(float)*_receptor->num_atoms()*3 ); //Atom coordinate

    // Ligand
    _ligand = new Ligand(_parameter->_LigPDB_file);
    _ligand->initialize(_parameter);
    _cputime->record_malloc( sizeof(float)*_ligand->num_atoms()*3 ); //Atom coordinate

    if( !_parameter->_Num_fft_flag ) {
        switch (_parameter->fft_base_set) {
        case 13:
            gridtable_13base_normal(ngrid,ngrid_table);
            break;
        case 7:
            gridtable_07base_normal(ngrid,ngrid_table);
            break;
        case 11:
            gridtable_11base_normal(ngrid,ngrid_table);
            break;
        case 0:
            gridtable_fftw_custom(ngrid,ngrid_table);
            break;
        case 1:
            gridtable_cufft_custom(ngrid,ngrid_table);
            break;
        }
        autogridr(ngrid,ngrid_table);
        autogridl(ngrid,ngrid_table);
    } else {
        checkgridr();
        checkgridl();
    }

    // Docking
    _docking = new Docking(_cputime,_parallel,_parameter,_receptor,_ligand);
    _docking->initialize();

    gettimeofday(&et2,NULL);
    _cputime->t1_initialize += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    return;
}

//============================================================================//
void Control::gridtable_11base_normal(int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{

    int grid_table[] = {2,3,4,5,6,7,8,9,10,11,12,14,15,16,18,20,21,22,24,25,27,
                        28,30,32,33,35,36,40,42,44,45,48,49,50,54,55,56,60,63,
                        64,66,70,72,75,77,80,81,84,88,90,96,98,99,100,105,108,
                        110,112,120,121,125,126,128,132,135,140,144,147,150,
                        154,160,162,165,168,175,176,180,189,192,196,198,200,
                        210,216,220,224,225,231,240,242,243,245,250,252,256,
                        264,270,275,280,288,294,297,300,308,315,320,324,330,
                        336,343,350,352,360,363,375,378,384,385,392,396,400,
                        405,420,432,440,441,448,450,462,480,484,486,490,495,
                        500,504,512,525,528,539,540,550,560,567,576,588,594,
                        600,605,616,625,630,640,648,660,672,675,686,693,700,
                        704,720,726,729,735,750,756,768,770,784,792,800,810,
                        825,840,847,864,875,880,882,891,896,900,924,945,960,
                        968,972,980,990,1000,1008,1024
                       };

    ngrid = sizeof(grid_table)/sizeof(grid_table[0]);
    ngrid_table.resize(ngrid);
    if( ngrid_table.size() < ngrid ) {
        cerr << "[ERROR] Out of memory. Number of FFT table lists = ("
             << ngrid << ") for (ngrid_table) in control.cpp!!\n";
        exit(1);
    }

    for( int i = 0 ; i < ngrid ; i++ ) {
        ngrid_table[i] = grid_table[i];
    }

    return;
}

//============================================================================//
void Control::gridtable_13base_normal(int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{

    int grid_table[] = { 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 ,
                         12 , 13 , 14 , 15 , 16 , 18 , 20 , 21 , 22 , 24 ,
                         25 , 26 , 27 , 28 , 30 , 32 , 33 , 35 , 36 , 39 ,
                         40 , 42 , 44 , 45 , 48 , 49 , 50 , 52 , 54 , 55 ,
                         56 , 60 , 63 , 64 , 65 , 66 , 70 , 72 , 75 , 77 ,
                         78 , 80 , 81 , 84 , 88 , 90 , 91 , 96 , 98 , 99 ,
                         100 , 104 , 105 , 108 , 110 , 112 , 117 , 120 , 121 , 125 ,
                         126 , 128 , 130 , 132 , 135 , 140 , 143 , 144 , 147 , 150 ,
                         154 , 156 , 160 , 162 , 165 , 168 , 169 , 175 , 176 , 180 ,
                         182 , 189 , 192 , 195 , 196 , 198 , 200 , 208 , 210 , 216 ,
                         220 , 224 , 225 , 231 , 234 , 240 , 242 , 243 , 245 , 250 ,
                         252 , 256 , 260 , 264 , 270 , 273 , 275 , 280 , 286 , 288 ,
                         294 , 297 , 300 , 308 , 312 , 315 , 320 , 324 , 325 , 330 ,
                         336 , 338 , 343 , 350 , 351 , 352 , 360 , 363 , 364 , 375 ,
                         378 , 384 , 385 , 390 , 392 , 396 , 400 , 405 , 416 , 420 ,
                         429 , 432 , 440 , 441 , 448 , 450 , 455 , 462 , 468 , 480 ,
                         484 , 486 , 490 , 495 , 500 , 504 , 507 , 512 , 520 , 525 ,
                         528 , 539 , 540 , 546 , 550 , 560 , 567 , 572 , 576 , 585 ,
                         588 , 594 , 600 , 605 , 616 , 624 , 625 , 630 , 637 , 640 ,
                         648 , 650 , 660 , 672 , 675 , 676 , 686 , 693 , 700 , 702 ,
                         704 , 715 , 720 , 726 , 728 , 729 , 735 , 750 , 756 , 768 ,
                         770 , 780 , 784 , 792 , 800 , 810 , 819 , 825 , 832 , 840 ,
                         845 , 847 , 858 , 864 , 875 , 880 , 882 , 891 , 896 , 900 ,
                         910 , 924 , 936 , 945 , 960 , 968 , 972 , 975 , 980 , 990 ,
                         1000 , 1001 , 1008 , 1014 , 1024
                       };

    ngrid = sizeof(grid_table)/sizeof(grid_table[0]);
    ngrid_table.resize(ngrid);
    if( ngrid_table.size() < ngrid ) {
        cerr << "[ERROR] Out of memory. Number of FFT table lists = ("
             << ngrid << ") for (ngrid_table) in control.cpp!!\n";
        exit(1);
    }

    for( int i = 0 ; i < ngrid ; i++ ) {
        ngrid_table[i] = grid_table[i];
    }

    return;
}

//============================================================================//
void Control::gridtable_07base_normal(int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{

    int grid_table[] = { 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 12 ,
                         14 , 15 , 16 , 18 , 20 , 21 , 24 , 25 , 27 , 28 ,
                         30 , 32 , 35 , 36 , 40 , 42 , 45 , 48 , 49 , 50 ,
                         54 , 56 , 60 , 63 , 64 , 70 , 72 , 75 , 80 , 81 ,
                         84 , 90 , 96 , 98 , 100 , 105 , 108 , 112 , 120 , 125 ,
                         126 , 128 , 135 , 140 , 144 , 147 , 150 , 160 , 162 , 168 ,
                         175 , 180 , 189 , 192 , 196 , 200 , 210 , 216 , 224 , 225 ,
                         240 , 243 , 245 , 250 , 252 , 256 , 270 , 280 , 288 , 294 ,
                         300 , 315 , 320 , 324 , 336 , 343 , 350 , 360 , 375 , 378 ,
                         384 , 392 , 400 , 405 , 420 , 432 , 441 , 448 , 450 , 480 ,
                         486 , 490 , 500 , 504 , 512 , 525 , 540 , 560 , 567 , 576 ,
                         588 , 600 , 625 , 630 , 640 , 648 , 672 , 675 , 686 , 700 ,
                         720 , 729 , 735 , 750 , 756 , 768 , 784 , 800 , 810 , 840 ,
                         864 , 875 , 882 , 896 , 900 , 945 , 960 , 972 , 980 , 1000 ,
                         1008 , 1024
                       };


    ngrid = sizeof(grid_table)/sizeof(grid_table[0]);
    ngrid_table.resize(ngrid);
    if( ngrid_table.size() < ngrid ) {
        cerr << "[ERROR] Out of memory. Number of FFT table lists = ("
             << ngrid << ") for (ngrid_table) in control.cpp!!\n";
        exit(1);
    }

    for( int i = 0 ; i < ngrid ; i++ ) {
        ngrid_table[i] = grid_table[i];
    }

    return;
}

//============================================================================//
void Control::gridtable_fftw_custom(int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{

    int grid_table[] = { 8,10,16,18,21,22,24,25,32,33,35,
                         36,39,40,42,44,45,49,50,52,55,56,60,64,70,75,78,84,91,
                         100,105,110,117,126,130,140,147,150,154,156,162,165,169,175,189,200,210,216,220,225,
                         231 , 234 , 240 , 242 , 243 , 245 , 250 ,
                         252 , 256 , 260 , 264 , 270 , 273 , 275 , 280 , 286 , 288 ,
                         294 , 297 , 300 , 308 , 312 , 315 , 320 , 324 , 325 , 330 ,
                         336 , 338 , 343 , 350 , 351 , 352 , 360 , 363 , 364 , 375 ,
                         378 , 384 , 385 , 390 , 392 , 396 , 400 , 405 , 416 , 420 ,
                         429 , 432 , 440 , 441 , 448 , 450 , 455 , 462 , 468 , 480 ,
                         484 , 486 , 490 , 495 , 500 , 504 , 507 , 512 , 520 , 525 ,
                         528 , 539 , 540 , 546 , 550 , 560 , 567 , 572 , 576 , 585 ,
                         588 , 594 , 600 , 605 , 616 , 624 , 625 , 630 , 637 , 640 ,
                         648 , 650 , 660 , 672 , 675 , 676 , 686 , 693 , 700 , 702 ,
                         704 , 715 , 720 , 726 , 728 , 729 , 735 , 750 , 756 , 768 ,
                         770 , 780 , 784 , 792 , 800 , 810 , 819 , 825 , 832 , 840 ,
                         845 , 847 , 858 , 864 , 875 , 880 , 882 , 891 , 896 , 900 ,
                         910 , 924 , 936 , 945 , 960 , 968 , 972 , 975 , 980 , 990 ,
                         1000 , 1001 , 1008 , 1014 , 1024
                       }; //selected the just values to have faster FFT calculation time (in n < 226)


    ngrid = sizeof(grid_table)/sizeof(grid_table[0]);
    ngrid_table.resize(ngrid);
    if( ngrid_table.size() < ngrid ) {
        cerr << "[ERROR] Out of memory. Number of FFT table lists = ("
             << ngrid << ") for (ngrid_table) in control.cpp!!\n";
        exit(1);
    }

    for( int i = 0 ; i < ngrid ; i++ ) {
        ngrid_table[i] = grid_table[i];
    }

    return;
}

//============================================================================//
void Control::gridtable_cufft_custom(int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{

    int grid_table[] = { 16,32,36,40,42,48,50,64,72,80,81,84,96,98,100,108,112,128,144,160,162,168,200,216,224,240 ,
                         243 , 245 , 250 , 252 , 256 , 270 , 280 , 288 , 294 ,
                         300 , 315 , 320 , 324 , 336 , 343 , 350 , 360 , 375 , 378 ,
                         384 , 392 , 400 , 405 , 420 , 432 , 441 , 448 , 450 , 480 ,
                         486 , 490 , 500 , 504 , 512 , 525 , 540 , 560 , 567 , 576 ,
                         588 , 600 , 625 , 630 , 640 , 648 , 672 , 675 , 686 , 700 ,
                         720 , 729 , 735 , 750 , 756 , 768 , 784 , 800 , 810 , 840 ,
                         864 , 875 , 882 , 896 , 900 , 945 , 960 , 972 , 980 , 1000 ,
                         1008 , 1024
                       }; //selected the just values to have faster FFT calculation time (in n < 226)


    ngrid = sizeof(grid_table)/sizeof(grid_table[0]);
    ngrid_table.resize(ngrid);
    if( ngrid_table.size() < ngrid ) {
        cerr << "[ERROR] Out of memory. Number of FFT table lists = ("
             << ngrid << ") for (ngrid_table) in control.cpp!!\n";
        exit(1);
    }

    for( int i = 0 ; i < ngrid ; i++ ) {
        ngrid_table[i] = grid_table[i];
    }

    return;
}



//============================================================================//
void Control::autogridr(const int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{
    int       num_grid = 1;
    float     size, size_rec = 0.0;

    for( int i = 0 ; i < 3 ; i++ ) {
        size = _receptor->edge(i,1) - _receptor->edge(i,0);
        
        //printf(" %f, %f\n",_receptor->edge(i,1),_receptor->edge(i,0));

        if( size > size_rec ) {
            size_rec = size;
        }
    }

    cout << "\nReceptor = " << _receptor->input_file() << endl;
    cout << "Receptor max size = " << size_rec << endl;

    size_rec += 2.0 * _parameter->_Grid_space_rec;
    cout << "Required voxel size = " << size_rec << endl;

    num_grid = 1 + int(size_rec / _parameter->grid_width);

    for( int i = 0 ; i < ngrid ; i++ ) {
        if( ngrid_table[i] >= num_grid ) {
            num_grid = ngrid_table[i];
            break;
        }
    }

    _receptor->num_grid(num_grid);

    cout << "Number of grid = " << num_grid << endl;
    cout << "FFT N = " << num_grid*2 << endl;

    return;
}

//============================================================================//
void Control::autogridl(const int &ngrid,vector<int> &ngrid_table)
//============================================================================//
{
    int       num_grid = 1;
    float     size_lig = 0.0;
    float     x1, y1, z1, x2, y2, z2, d2;
    const int na  = _ligand->num_atoms();

    for( int i = 0 ; i < na-1 ; i++ ) {
        x1 = _ligand->coordinate(i,0);
        y1 = _ligand->coordinate(i,1);
        z1 = _ligand->coordinate(i,2);

        for( int j = i+1 ; j < na ; j++ ) {
            x2 = _ligand->coordinate(j,0);
            y2 = _ligand->coordinate(j,1);
            z2 = _ligand->coordinate(j,2);

            d2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);

            if( d2 > size_lig ) {
                size_lig = d2;
            }
        }
    }

    size_lig = sqrt(size_lig);

    cout << "\nLigand = " << _ligand->input_file() << endl;
    cout << "Ligand max size = " << size_lig << endl;

    size_lig += 2.0 * _parameter->_Grid_space_lig;
    
    _parameter->ligand_max_edge = size_lig;

    cout << "Required voxel size = " << size_lig << endl;

    num_grid = 1 + int(size_lig / _parameter->grid_width);

    for( int i = 0 ; i < ngrid ; i++ ) {
        if( ngrid_table[i] >= num_grid ) {
            num_grid = ngrid_table[i];
            break;
        }
    }

    _ligand->num_grid(num_grid);

    cout << "Number of grid = " << num_grid << endl;
    cout << "FFT N = " << num_grid*2 << endl;

    return;
}

//============================================================================//
void Control::checkgridr()
//============================================================================//
{
    float     size, size_rec  = 0.0;
    const int num_grid    = _parameter->_Num_grid;
    const float   search_length   = _parameter->grid_width * num_grid;

    for( int i = 0 ; i < 3 ; i++ ) {
        size = _receptor->edge(i,1) - _receptor->edge(i,0);

        if( size > size_rec ) {
            size_rec = size;
        }
    }

    cout << "\nReceptor max size = " << size_rec << endl;

    size_rec += 2.0*_parameter->_Grid_space_rec;

    cout << "Required voxel size = " << size_rec << endl;

    if( size_rec > search_length ) {
        cerr << "[ERROR] Receptor data is too big!!\n";
        exit(1);
    }

    _receptor->num_grid(num_grid);

    cout << "\n(Receptor)\n";
    cout << "Number of grid = " << num_grid << endl;
    cout << "FFT N = " << num_grid*2 << endl;
    cout << "Grid size = " << _parameter->grid_width << endl;

    return;
}

//============================================================================//
void Control::checkgridl()
//============================================================================//
{
    float     size_lig = 0.0;
    float     x1, y1, z1, x2, y2, z2, d2;
    const int na  = _ligand->num_atoms();
    const int num_grid    = _parameter->_Num_grid;
    const float   search_length   = _parameter->grid_width * num_grid;

    for( int i = 0 ; i < na-1 ; i++ ) {
        x1 = _ligand->coordinate(i,0);
        y1 = _ligand->coordinate(i,1);
        z1 = _ligand->coordinate(i,2);

        for( int j = i+1 ; j < na ; j++ ) {
            x2 = _ligand->coordinate(j,0);
            y2 = _ligand->coordinate(j,1);
            z2 = _ligand->coordinate(j,2);

            d2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);

            if( d2 > size_lig ) {
                size_lig = d2;
            }
        }
    }

    size_lig = sqrt(size_lig);
    cout << "\nLigand max size = " << size_lig << endl;

    size_lig += 2.0*_parameter->_Grid_space_lig;
    cout << "Required voxel size = " << size_lig << endl;

    if( size_lig > search_length ) {
        cerr << "[ERROR] Ligand data is too big!!\n";
        exit(1);
    }

    _ligand->num_grid(num_grid);

    cout << "\n(Ligand)\n";
    cout << "Number of grid = " << num_grid << endl;
    cout << "FFT N = " << num_grid*2 << endl;
    cout << "Grid size = " << _parameter->grid_width << endl;

    return;
}

//============================================================================//
void Control::execute()
//============================================================================//
{
    struct timeval et1, et2;

    cout << "\n---------- Start docking calculations" << endl;

    gettimeofday(&et1,NULL); // Receptor process (voxelization, forward FFT of Receptor)
    _docking->rec_init();
    gettimeofday(&et2,NULL);
    _cputime->t2_receptor_process += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    gettimeofday(&et1,NULL); // docking
    _docking->dockz();
    gettimeofday(&et2,NULL);
    _cputime->t3_docking_total += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    if(_parameter->detail_output_flag == 1) { // detailed result output
        gettimeofday(&et1,NULL);
        _docking->output_detail();
        gettimeofday(&et2,NULL);
        _cputime->t4_docking_output_detail += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));
    }

    if(_parameter->calc_time_log_output_flag >= 1) { // calculation info
        _docking->output_calc_time_log();
    }

    gettimeofday(&et1,NULL); // normal result output
    _docking->output();
    gettimeofday(&et2,NULL);
    _cputime->t5_docking_output += (et2.tv_sec-et1.tv_sec + (float)((et2.tv_usec-et1.tv_usec)*1e-6));

    _docking->dock_memory_free();

    const int ng1 = _parameter->_Num_grid;
    const int ng3 = ng1*ng1*ng1;
    const int nf1 = ng1*2;
    const int nf3 = nf1*nf1*nf1;
    const int nproc2 = _parallel->nproc2();
    const int natom = _parameter->_Num_atom_max;
    const int nag = natom * ng1;
    const size_t _Memfw = ng3*3+natom*3+nag*3;
    const size_t _Memiw = ng3*2+natom*4;

    //delete docking include delete fft_process, _FFT_rec_r/i[nf3], _FFTWin/out[nf3*nproc2]
    _cputime->record_free( sizeof(float)*nf3*2 + sizeof(fftwf_complex)*nf3*2*nproc2);
#ifdef CUFFT
    _cputime->record_free( sizeof(cufftComplex)*nf3*2 ); //_in/outBuf
#endif
    _cputime->record_free( sizeof(float)*_Memfw*nproc2 + sizeof(int)*_Memiw*nproc2 ); //_F/Iwork
    delete _docking;
    _cputime->record_free( sizeof(float)*_ligand->num_atoms()*3 );
    delete _ligand;
    _cputime->record_free( sizeof(float)*_receptor->num_atoms()*3 );
    delete _receptor;
    _cputime->record_free( sizeof(float)*_parameter->_Num_rot_angles*3 + sizeof(map<string,float>)*(_parameter->_Charmmr.size() + _parameter->_Charmmc.size() + _parameter->_ACE.size()) );
    delete _parameter;


    return;
}
