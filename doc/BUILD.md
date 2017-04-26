# Build Documentation

## Required Tools
* FFTW3
  - http://www.fftw.org
* OpenMPI (if you use MPI)
  - http://www.open-mpi.org
* CUDA Toolkit ver >= 5.0 (if you use GPU)
  - https://developer.nvidia.com/cuda-zone
* GPU Computing SDK code samples (CUDA SDK, **same version as CUDA Toolkit**, if you use GPU)
  - https://developer.nvidia.com/cuda-zone


--------------------------------------------------------


## Compile (a): GPU, MPI & OpenMP hybrid parallelization

### 1. Extract tarball contents
```sh
tar xzf megadock-4.0.tgz
cd megadock-4.0
```

### 2. Edit Makefile
```Makefile
CUDA_INSTALL_PATH ?= your/cuda/toolkit/install/path (default: /usr/local/cuda )
CUDA_SAMPLES_PATH ?= your/cuda/sdk/install/path     (default: ${HOME}/samples )
FFTW_INSTALL_PATH ?= your/fftw/library/install/path (default: /usr/local      )
CPPCOMPILER       ?= icpc, g++ or others
MPICOMPILER       ?= mpicxx or others
OPTIMIZATION      ?= -O3
OMPFLAG           ?= -openmp (intel) or -fopenmp (g++)
```

### 3. Compile 
```sh
make
```

A binary file `megadock-gpu-dp` will be generated. 


--------------------------------------------------------


## Compile (b): MPI & OpenMP hybrid parallelization (no use GPU)

### 1. Extract tarball contents
```sh
tar xzf megadock-4.0.tgz
cd megadock-4.0
```

### 2. Edit Makefile
```Makefile
FFTW_INSTALL_PATH ?= your/fftw/library/install/path (default: /usr/local )
CPPCOMPILER       ?= icpc, g++ or others
MPICOMPILER       ?= mpicxx or others
OPTIMIZATION      ?= -O3
OMPFLAG           ?= -openmp (intel) or -fopenmp (g++)

USE_GPU := 0
```

### 3. Compile
```sh
make
```

A binary file `megadock-dp` will be generated.


--------------------------------------------------------


## Compile (c): GPU parallelization (on single node)

### 1. Extract tarball contents
```sh
tar xzf megadock-4.0.tgz
cd megadock-4.0
```

### 2. Edit Makefile
```Makefile
CUDA_INSTALL_PATH ?= your/cuda/toolkit/install/path (default: /usr/local/cuda )
CUDA_SAMPLES_PATH ?= your/cuda/sdk/install/path     (default: ${HOME}/samples )
FFTW_INSTALL_PATH ?= your/fftw/library/install/path (default: /usr/local      )
CPPCOMPILER       ?= icpc, g++ or others
OPTIMIZATION      ?= -O3
OMPFLAG           ?= -openmp (intel) or -fopenmp (g++)

USE_MPI := 0
```

### 3. Compile
```sh
make
```

A binary file `megadock-gpu` will be generated.


--------------------------------------------------------



## Compile (d): CPU single node (only thread parallelization)

### 1. Extract tarball contents
```sh
tar xzf megadock-4.0.tgz
cd megadock-4.0
```

### 2. Edit Makefile
```Makefile
FFTW_INSTALL_PATH ?= your/fftw/library/install/path (default: /usr/local )
CPPCOMPILER       ?= icpc, g++ or others
OPTIMIZATION      ?= -O3
OMPFLAG           ?= -openmp (intel) or -fopenmp (g++)

USE_GPU := 0
USE_MPI := 0
```

### 3. Compile
```sh
make
```

A binary file `megadock` will be generated.


--------------------------------------------------------


## Clean up all files

```sh
make allclean
```


