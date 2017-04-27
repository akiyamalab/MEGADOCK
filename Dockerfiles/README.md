# Dockerfiles

Docker images can be downloaded from DockerHub.  
[https://hub.docker.com/r/akiyamalab/megadock/](https://hub.docker.com/r/akiyamalab/megadock/)

## Build Requirements
| Requirement Tools                                         | GPU cluster | CPU cluster | GPU | CPU | Notes       |
|:----------------------------------------------------------|:-----------:|:-----------:|:---:|:---:|:------------|
| [Docker](https://docs.docker.com/engine/installation/)    | N/A         | x           | x   | x   |             |
| [nvidia-docker](https://github.com/NVIDIA/nvidia-docker)  | N/A         |             | x   |     | for GPU use |
| [CUDA Toolkit](https://developer.nvidia.com/cuda-zone)    | N/A         |             | x   |     | for GPU use |


## (a): GPU, MPI & OpenMP hybrid parallelization

**Note: Currently NOT Available**

## (b): MPI & OpenMP hybrid parallelization (no use GPU)

**Note: Currently DO NOT support multi-node execution**  
If you use this Docker image for running on multi-node, please configure `sshd` and use container orchestration tools. (e.g. Docker Swarm, Kubernetes, Apache Mesos, etc.)

### 1. build Docker image
```sh
# on MEGADOCK_ROOT dir
docker build . -f Dockerfiles/Dockerfile.mpi -t megadock:mpi
```

### 2. run sample
```sh
docker run -it  \
    megadock:mpi  \
    mpirun -n 4 megadock-dp -tb SAMPLE.table

# run with your pdb ( ${DATA_PATH} = your pdb-data directory abs path  )
docker run -it  \
    -v ${DATA_PATH}:/opt/MEGADOCK/data \
    megadock:mpi  \
    mpirun -n 4 megadock-dp -tb ${PDB_TABLE}.table
```


## (c): GPU parallelization (on single node)

**[nvidia-docker](https://github.com/NVIDIA/nvidia-docker) is required.**

### 1. build Docker image
```sh
# on MEGADOCK_ROOT dir
nvidia-docker build . -f Dockerfiles/Dockerfile.gpu -t megadock:gpu
```

### 2. run sample
```sh
nvidia-docker run -it  \
    megadock:gpu  \
    megadock-gpu -R 1gcq_r.pdb -L 1gcq_r.pdb

# run with your pdb ( ${DATA_PATH} = your pdb-data directory abs path  )
nvidia-docker run -it  \
    -v ${DATA_PATH}:/opt/MEGADOCK/data \
    megadock:gpu  \
    megadock-gpu -R ${RECEPTOR}.pdb -L ${LIGAND}.pdb
```


## (d): CPU single node (only thread parallelization)

### 1. build Docker image
```sh
# on MEGADOCK_ROOT dir
docker build . -f Dockerfiles/Dockerfile.cpu -t megadock:cpu
```

### 2. run sample
```sh
docker run -it  \
    megadock:cpu  \
    megadock -R 1gcq_r.pdb -L 1gcq_r.pdb

# run with your pdb ( ${DATA_PATH} = your pdb-data directory abs path  )
docker run -it \
    -v ${DATA_PATH}:/opt/MEGADOCK/data \
    megadock:cpu  \
    megadock -R ${RECEPTOR}.pdb -L ${LIGAND}.pdb
```

