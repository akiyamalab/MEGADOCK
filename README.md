# MEGADOCK 4.0

**MEGADOCK 4.0** is a structural bioinformatics software for FFT-grid-based protein-protein docking that takes advantage of the massively parallel CUDA architechture of NVIDIA GPUs and multiple computation nodes. 

## Requirements
- [FFTW3](http://www.fftw.org)
  - `--enable-float` flag must be specified when you compile FFTW3 (see also FAQ).
- [OpenMPI](http://www.open-mpi.org) (if you use MPI)
- [CUDA Toolkit](https://developer.nvidia.com/cuda-zone) ver >= 5.0 (if you use GPU)
- [GPU Computing SDK code samples](https://developer.nvidia.com/cuda-zone) (if you use GPU)
  - CUDA SDK, same version as CUDA Toolkit

## Installation
You can use MEGADOCK 4.0 on various environment.
- (a) GPU cluster
- (b) CPU cluster
- (c) GPU single node 
- (d) CPU single node

Please refer to the following URL and see the appropriate instructions.  
http://www.bi.cs.titech.ac.jp/megadock/install.html

## Reference
Masahito Ohue, Takehiro Shimoda, Shuji Suzuki, Yuri Matsuzaki, Takashi Ishida, Yutaka Akiyama. MEGADOCK 4.0: an ultra-high-performance protein-protein docking software for heterogeneous supercomputers, Bioinformatics, 30(22): 3281-3283, 2014.

## Lisence
MEGADOCK is open source licensed under the GNU General Public License, version 3 or later. 

## Fundings
This work is partially supported by JSPS Grant-in-Aid for Scientific Research (KAKENHI) (A) Grant Number 24240044.

----
Copyright © 2014-2015 Akiyama Laboratory, Tokyo Institute of Technology, All Rights Reserved.

