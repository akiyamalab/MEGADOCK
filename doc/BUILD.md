# Installation Guide


## Target Environments

| Type | Target Env.     | Approach           |
|:----:|-----------------|--------------------|
|  (a) | GPU cluster     | GPU + OpenMP + MPI |
|  (b) | CPU cluster     | OpenMP + MPI       |
|  (c) | GPU node        | GPU + OpenMP       |
|  (d) | CPU node        | OpenMP             |


## Requirement

| Library                                                         | GPU cluster | CPU cluster | GPU node | CPU node | Note |
|:----------------------------------------------------------------|:-----------:|:-----------:|:---:|:---:|:------|
| [FFTW3](http://www.fftw.org)                                    | yes           | yes           | yes   | yes   | `--enable-float` required ([FAQ](http://www.bi.cs.titech.ac.jp/megadock/faq.html))|
| MPI (e.g. [OpenMPI](http://www.open-mpi.org))                   | yes           | yes           |     |     |       |
| [CUDA Toolkit](https://developer.nvidia.com/cuda-zone)          | yes           |             | yes   |     | `ver >= 5.0` required |
| [CUDA SDK code samples](https://developer.nvidia.com/cuda-zone) | yes           |             | yes   |     | same version as toolkit |

For more detailed build information, please reffer to the [FAQ](http://www.bi.cs.titech.ac.jp/megadock/faq.html) page.


## Section Link

Please select appropriate section for your environment:

- [(a) Compile for GPU cluster (GPU & MPI)](#a-compile-for-gpu-cluster-gpu--mpi)
- [(b) Compile for MPI cluster (MPI)](#b-compile-for-mpi-cluster-mpi)
- [(c) Compile for GPU node (GPU)](#c-compile-for-gpu-node-gpu)
- [(d) Compile for CPU node (only thread parallelization)](#d-compile-for-cpu-node-only-thread-parallelization)
- [Note for compilation](#note-for-compilation)
----

# (a) Compile for GPU cluster (GPU & MPI)

## 1. Get source code and change directory

```sh
# from GitHub (clone)
git clone https://github.com/akiyamalab/MEGADOCK.git
cd MEGADOCK

# from GitHub (zip)
wget https://github.com/akiyamalab/MEGADOCK/archive/master.zip
unzip master.zip
cd MEGADOCK-master
```

## 2. Edit Makefile (PATH, options)

```Makefile
CUDA_INSTALL_PATH ?= your/cuda/toolkit/install/path # default: /usr/local/cuda
CUDA_SAMPLES_PATH ?= your/cuda/sdk/install/path     # default: /usr/local/cuda/samples
FFTW_INSTALL_PATH ?= your/fftw/library/install/path # default: /usr/local
CPPCOMPILER       ?= g++                            # default: g++ (or icpc, others)
MPICOMPILER       ?= mpicxx                         # default: mpicxx (or others)
OPTIMIZATION      ?= -O3                            # 
OMPFLAG           ?= -fopenmp                       # default: -fopenmp (g++), or -openmp, -qopenmp (intel)
```

```Makefile
USE_GPU := 1
USE_MPI := 1
```

## 3. generate binary

After completion of `make` command, a binary file `megadock-gpu-dp` will be generated.

```sh
# build a binary
make

# see help
megadock-gpu-dp -h
```

----

# (b) Compile for MPI cluster (MPI)


## 1. Get source code and change directory

```sh
# from GitHub (clone)
git clone https://github.com/akiyamalab/MEGADOCK.git
cd MEGADOCK

# from GitHub (zip)
wget https://github.com/akiyamalab/MEGADOCK/archive/master.zip
unzip master.zip
cd MEGADOCK-master
```

## 2. Edit Makefile (PATH, options)

```Makefile
FFTW_INSTALL_PATH ?= your/fftw/library/install/path # default: /usr/local
CPPCOMPILER       ?= g++                            # default: g++ (or icpc, others)
MPICOMPILER       ?= mpicxx                         # default: mpicxx (or others)
OPTIMIZATION      ?= -O3                            # 
OMPFLAG           ?= -fopenmp                       # default: -fopenmp (g++), or -openmp, -qopenmp (intel)
```

```Makefile
USE_GPU := 0
USE_MPI := 1
```

## 3. generate binary

After completion of `make` command, a binary file `megadock-dp` will be generated.

```sh
# build a binary
make

# see help
megadock-dp -h
```

----

# (c) Compile for GPU node (GPU)

## 1. Get source code and change directory
```sh
# from GitHub (clone)
git clone https://github.com/akiyamalab/MEGADOCK.git
cd MEGADOCK

# from GitHub (zip)
wget https://github.com/akiyamalab/MEGADOCK/archive/master.zip
unzip master.zip
cd MEGADOCK-master
```

## 2. Edit Makefile (PATH, options)

```Makefile
CUDA_INSTALL_PATH ?= your/cuda/toolkit/install/path # default: /usr/local/cuda
CUDA_SAMPLES_PATH ?= your/cuda/sdk/install/path     # default: /usr/local/cuda/samples
FFTW_INSTALL_PATH ?= your/fftw/library/install/path # default: /usr/local
CPPCOMPILER       ?= g++                            # default: g++ (or icpc, others)
OPTIMIZATION      ?= -O3                            # 
OMPFLAG           ?= -fopenmp                       # default: -fopenmp (g++), or -openmp, -qopenmp (intel)
```

```Makefile
USE_GPU := 1
USE_MPI := 0
```

## 3. generate binary

After completion of `make` command, a binary file `megadock-gpu` will be generated.

```sh
# build a binary
make

# see help
megadock-gpu -h
```

--------------------------------------------------------

# (d) Compile for CPU node (only thread parallelization)

## 1. Get source code and change directory
```sh
# from GitHub (clone)
git clone https://github.com/akiyamalab/MEGADOCK.git
cd MEGADOCK

# from GitHub (zip)
wget https://github.com/akiyamalab/MEGADOCK/archive/master.zip
unzip master.zip
cd MEGADOCK-master
```

## 2. Edit Makefile (PATH, options)

```Makefile
FFTW_INSTALL_PATH ?= your/fftw/library/install/path # default: /usr/local
CPPCOMPILER       ?= g++                            # default: g++ (or icpc, others)
OPTIMIZATION      ?= -O3                            # 
OMPFLAG           ?= -fopenmp                       # default: -fopenmp (g++), or -openmp, -qopenmp (intel)
```

```
USE_GPU := 0
USE_MPI := 0
```

## 3. generate binary

After completion of `make` command, a binary file `megadock` will be generated.

```sh
# build a binary
make

# see help
megadock -h
```

----

## Clean up

```sh
# remove generated files
make allclean
```

----

# Note for compilation

- `USE_MPI` and `USE_GPU` flags in Makefile should be `1` or `0` for all compilations. If the white spaces exist at end of the line, those flags will be ignored.
- For GPU use, `SM_VERSIONS` can specity the target NVIDIA GPU architectures of the generated binary. Please change it for your target system. (Default: `sm_60`)
