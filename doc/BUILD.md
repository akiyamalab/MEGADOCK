# Build Documentation

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
| [FFTW3](http://www.fftw.org)                                    | x           | x           | x   | x   | `--enable-float` flag required in compile (see [FAQ](http://www.bi.cs.titech.ac.jp/megadock/faq.html))|
| [OpenMPI](http://www.open-mpi.org)                              | x           | x           |     |     |       |
| [CUDA Toolkit](https://developer.nvidia.com/cuda-zone)          | x           |             | x   |     | `ver >= 5.0` required |
| [CUDA SDK code samples](https://developer.nvidia.com/cuda-zone) | x           |             | x   |     | **same version as CUDA Toolkit required** |

#### Notes
For more detailed information, please reffer to the [FAQ](http://www.bi.cs.titech.ac.jp/megadock/faq.html) page.

--------------------------------------------------------


## 1. Get source code and change directory
```sh
# from GitHub (clone)
git clone https://github.com/akiyamalab/MEGADOCK.git
cd MEGADOCK

# from GitHub (zip)
wget https://github.com/akiyamalab/MEGADOCK/archive/master.zip
unzip master.zip
cd MEGADOCK-master

# from webpage (tar)
wget http://www.bi.cs.titech.ac.jp/megadock/archives/megadock-4.1.0.tgz
tar xvzf megadock-4.1.0.tgz
cd megadock-4.1.0
```


## 2-a: Compile for GPU, MPI & OpenMP hybrid parallelization

### Edit Makefile (PATH, options)

#### PATH
```Makefile
CUDA_INSTALL_PATH ?= your/cuda/toolkit/install/path # default: /usr/local/cuda
CUDA_SAMPLES_PATH ?= your/cuda/sdk/install/path     # default: /usr/local/cuda/samples
FFTW_INSTALL_PATH ?= your/fftw/library/install/path # default: /usr/local
CPPCOMPILER       ?= g++                            # default: g++ (or icpc, others)
MPICOMPILER       ?= mpicxx                         # default: mpicxx (or others)
OPTIMIZATION      ?= -O3                            # 
OMPFLAG           ?= -fopenmp                       # default: -fopenmp (g++), or -openmp, -qopenmp (intel)
```
#### Options
```
USE_GPU := 1
USE_MPI := 1
```

### make
After completion of `make` command, a binary file `megadock-gpu-dp` will be generated.
```sh
make
```

--------------------------------------------------------

## 2-b: Compile for MPI & OpenMP hybrid parallelization (no use GPU)

### Edit Makefile (PATH, options)

#### PATH
```Makefile
FFTW_INSTALL_PATH ?= your/fftw/library/install/path # default: /usr/local
CPPCOMPILER       ?= g++                            # default: g++ (or icpc, others)
MPICOMPILER       ?= mpicxx                         # default: mpicxx (or others)
OPTIMIZATION      ?= -O3                            # 
OMPFLAG           ?= -fopenmp                       # default: -fopenmp (g++), or -openmp, -qopenmp (intel)
```
#### Options
```
USE_GPU := 0
USE_MPI := 1
```

### make
After completion of `make` command, a binary file `megadock-dp` will be generated.
```sh
make
```

--------------------------------------------------------

## 2-c: Compile for GPU parallelization (on single node)

### Edit Makefile (PATH, options)

#### PATH
```Makefile
CUDA_INSTALL_PATH ?= your/cuda/toolkit/install/path # default: /usr/local/cuda
CUDA_SAMPLES_PATH ?= your/cuda/sdk/install/path     # default: /usr/local/cuda/samples
FFTW_INSTALL_PATH ?= your/fftw/library/install/path # default: /usr/local
CPPCOMPILER       ?= g++                            # default: g++ (or icpc, others)
OPTIMIZATION      ?= -O3                            # 
OMPFLAG           ?= -fopenmp                       # default: -fopenmp (g++), or -openmp, -qopenmp (intel)
```
#### Options
```
USE_GPU := 1
USE_MPI := 0
```

### make
After completion of `make` command, a binary file `megadock-gpu` will be generated.
```sh
make
```

--------------------------------------------------------

## 2-d: Compile for CPU single node (only thread parallelization)

### Edit Makefile (PATH, options)

#### PATH
```Makefile
FFTW_INSTALL_PATH ?= your/fftw/library/install/path # default: /usr/local
CPPCOMPILER       ?= g++                            # default: g++ (or icpc, others)
OPTIMIZATION      ?= -O3                            # 
OMPFLAG           ?= -fopenmp                       # default: -fopenmp (g++), or -openmp, -qopenmp (intel)
```
#### Options
```
USE_GPU := 0
USE_MPI := 0
```

### make
After completion of `make` command, a binary file `megadock` will be generated.
```sh
make
```

--------------------------------------------------------

## Clean up all files

```sh
make allclean
```


