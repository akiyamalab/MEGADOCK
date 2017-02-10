
Required Tools
===================================================================================

* FFTW3 <http://www.fftw.org>
* OpenMPI <http://www.open-mpi.org> (if you use MPI)
* CUDA Toolkit ver >= 5.0 <https://developer.nvidia.com/cuda-zone> (if you use GPU)
* GPU Computing SDK code samples (CUDA SDK, same version as CUDA Toolkit)
   <https://developer.nvidia.com/cuda-zone> (if you use GPU)


Compile (a): GPU, MPI & OpenMP hybrid parallelization
===================================================================================
* Extract tarball contents
	$ tar xzf megadock-4.0.tgz
	$ cd megadock-4.0

* Edit Makefile
    CUDA_INSTALL_PATH ?= your/cuda/toolkit/install/path (default: /usr/local/cuda )
    CUDA_SAMPLES_PATH ?= your/cuda/sdk/install/path     (default: ${HOME}/samples )
    FFTW_INSTALL_PATH ?= your/fftw/library/install/path (default: /usr/local      )
    CPPCOMPILER       ?= icpc, g++ or others
    MPICOMPILER       ?= mpicxx or others
    OPTIMIZATION      ?= -O3
    OMPFLAG           ?= -openmp (intel) or -fopenmp (g++)

* Compile 
	$ make

  A binary file 'megadock-gpu-dp' will be generated.


Compile (b): MPI & OpenMP hybrid parallelization (no use GPU)
===================================================================================
* Extract tarball contents
    $ tar xzf megadock-4.0.tgz
    $ cd megadock-4.0

* Edit Makefile
    FFTW_INSTALL_PATH ?= your/fftw/library/install/path (default: /usr/local )
    CPPCOMPILER       ?= icpc, g++ or others
    MPICOMPILER       ?= mpicxx or others
    OPTIMIZATION      ?= -O3
    OMPFLAG           ?= -openmp (intel) or -fopenmp (g++)

    USE_GPU := 0

* Compile
    $ make

  A binary file 'megadock-dp' will be generated.


Compile (c): GPU parallelization (on single node)
===================================================================================
* Extract tarball contents
    $ tar xzf megadock-4.0.tgz
    $ cd megadock-4.0

* Edit Makefile
    CUDA_INSTALL_PATH ?= your/cuda/toolkit/install/path (default: /usr/local/cuda )
    CUDA_SAMPLES_PATH ?= your/cuda/sdk/install/path     (default: ${HOME}/samples )
    FFTW_INSTALL_PATH ?= your/fftw/library/install/path (default: /usr/local      )
    CPPCOMPILER       ?= icpc, g++ or others
    OPTIMIZATION      ?= -O3
    OMPFLAG           ?= -openmp (intel) or -fopenmp (g++)

    USE_MPI := 0

* Compile
    $ make

  A binary file 'megadock-gpu' will be generated.


Compile (d): CPU single node (only thread parallelization)
===================================================================================
* Extract tarball contents
    $ tar xzf megadock-4.0.tgz
    $ cd megadock-4.0

* Edit Makefile
    FFTW_INSTALL_PATH ?= your/fftw/library/install/path (default: /usr/local )
    CPPCOMPILER       ?= icpc, g++ or others
    OPTIMIZATION      ?= -O3
    OMPFLAG           ?= -openmp (intel) or -fopenmp (g++)

    USE_GPU := 0
    USE_MPI := 0

* Compile
    $ make

  A binary file 'megadock' will be generated.



Clean up all files
===================================================================================
    $ make allclean


