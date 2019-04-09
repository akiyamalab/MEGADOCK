# MEGADOCK

**MEGADOCK** is a structural bioinformatics software for FFT-grid-based protein-protein using MPI/OpenMP/GPU parallelization.

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![Build Status](https://travis-ci.org/akiyamalab/MEGADOCK.svg?branch=master)](https://travis-ci.org/akiyamalab/MEGADOCK)


## Target Environments
| Type | Environment     | Parallelization   |
|:----:|-----------------|-------------------|
|  (a) | GPU cluster     | GPU, MPI + OpenMP |
|  (b) | CPU cluster     | MPI + OpenMP      |
|  (c) | GPU single node | GPU               |
|  (d) | CPU single node | OpenMP            |

## Requirements
| Libraries                                                       | GPU cluster | CPU cluster | GPU | CPU | Notes |
|:----------------------------------------------------------------|:-----------:|:-----------:|:---:|:---:|:------|
| [FFTW3](http://www.fftw.org)                                    | x           | x           | x   | x   | `--enable-float` flag required (see [FAQ](http://www.bi.cs.titech.ac.jp/megadock/faq.html))|
| [OpenMPI](http://www.open-mpi.org)                              | x           | x           |     |     |       |
| [CUDA Toolkit](https://developer.nvidia.com/cuda-zone)          | x           |             | x   |     | `ver >= 5.0` required |
| [CUDA SDK code samples](https://developer.nvidia.com/cuda-zone) | x           |             | x   |     | **same version as CUDA Toolkit** |

## Installation
For installation details, please read appropriate section on followings:
- [doc/BUILD.md](./doc/BUILD.md)
- [doc/README.md](./doc/README.md)

### MEGADOCK for Docker
- Prebuild image : [akiyamalab/megadock](https://hub.docker.com/r/akiyamalab/megadock/) (Docker Hub)
- Build from Dockerfile : [Dockerfiles/README.md](Dockerfiles/README.md)


## Reference
Masahito Ohue, Takehiro Shimoda, Shuji Suzuki, Yuri Matsuzaki, Takashi Ishida, Yutaka Akiyama. **MEGADOCK 4.0: an ultra-high-performance protein-protein docking software for heterogeneous supercomputers**, *Bioinformatics*, 30(22): 3281-3283, 2014. http://dx.doi.org/10.1093/bioinformatics/btu532

Masahito Ohue, Yuri Matsuzaki, Nobuyuki Uchikoga, Takashi Ishida, Yutaka Akiyama. **MEGADOCK: An all-to-all protein-protein interaction prediction system using tertiary structure data**, *Protein and Peptide Letters*, 21(8): 766-778, 2014. http://eurekaselect.com/112757


## Older Versions
For older versions are available here.    
[http://www.bi.cs.titech.ac.jp/megadock/archives/](http://www.bi.cs.titech.ac.jp/megadock/archives/)

## Lisence
MEGADOCK is open source licensed under the GNU General Public License, version 3 or later. 

## Fundings
This work is partially supported by JSPS Grant-in-Aid for Scientific Research (KAKENHI) (A) Grant Number 24240044.

----
Copyright © 2014-2018 Akiyama Laboratory, Tokyo Institute of Technology, All Rights Reserved.

