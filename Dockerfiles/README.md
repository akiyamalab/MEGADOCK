# Dockerfiles

## Requirement
* [Docker](https://docs.docker.com/engine/installation/)
* [nvidia-docker](https://github.com/NVIDIA/nvidia-docker) ( if you use GPU)


## Compile (a): GPU, MPI & OpenMP hybrid parallelization

**Note: Currently NOT Available**

## Compile (b): MPI & OpenMP hybrid parallelization (no use GPU)

**Note: Currently DO NOT support multi-node execution**  
If you use this Docker image on clusters, please configure `sshd` and use container orchestration tools. (e.g. Docker Swarm, Kubernetes, Apache Mesos, etc.)

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


## Compile (c): GPU parallelization (on single node)

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


## Compile (d): CPU single node (only thread parallelization)

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

