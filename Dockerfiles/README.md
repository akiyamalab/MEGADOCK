# Dockerfiles

Docker images can be downloaded from DockerHub.  
[https://hub.docker.com/r/akiyamalab/megadock/](https://hub.docker.com/r/akiyamalab/megadock/)

## Build Requirements
| Requirement Tools                                         | GPU cluster | CPU cluster | GPU | CPU | Notes       |
|:----------------------------------------------------------|:-----------:|:-----------:|:---:|:---:|:------------|
| [Docker](https://docs.docker.com/engine/installation/)    | N/A         | N/A           | x   | x   |             |
| [nvidia-docker](https://github.com/NVIDIA/nvidia-docker)  | N/A         | N/A           | x   |     | > 2.0 |
| [CUDA Toolkit](https://developer.nvidia.com/cuda-zone)    | N/A         | N/A            | x   |     |  |

## Simple Example
```sh
# CPU single node (OpenMP parallelization)
docker run -it akiyamalab/megadock:cpu  megadock -R 1gcq_r.pdb -L 1gcq_r.pdb

# GPU single node (GPU parallelization)
docker run -it --runtime=nvidia akiyamalab/megadock:gpu  megadock-gpu -R 1gcq_r.pdb -L 1gcq_r.pdb
```

----

## CPU single node (OpenMP parallelization)

### 1. build Docker image
```sh
# on ${MEGADOCK_ROOT} dir
docker build . -f Dockerfiles/cpu/Dockerfile -t akiyamalab/megadock:cpu
```

### 2. run sample
```sh
docker run -it akiyamalab/megadock:cpu \
    megadock -R 1gcq_r.pdb -L 1gcq_r.pdb

# run with your pdb ( ${DATA_PATH} = your pdb-data directory abs path  )
docker run -it -v ${DATA_PATH}:/opt/MEGADOCK/data akiyamalab/megadock:cpu \
    megadock -R ${RECEPTOR}.pdb -L ${LIGAND}.pdb
```

## GPU single node (GPU parallelization)

**[nvidia-docker](https://github.com/NVIDIA/nvidia-docker) is required.**

### 1. build Docker image
```sh
# on ${MEGADOCK_ROOT} dir
docker build . -f Dockerfiles/gpu/Dockerfile -t akiyamalab/megadock:gpu
```

### 2. run sample
```sh
docker run -it --runtime=nvidia akiyamalab/megadock:gpu \
    megadock-gpu -R 1gcq_r.pdb -L 1gcq_r.pdb

# run with your pdb ( ${DATA_PATH} = your pdb-data directory abs path  )
docker run -it --runtime=nvidia -v ${DATA_PATH}:/opt/MEGADOCK/data akiyamalab/megadock:gpu \
    megadock-gpu -R ${RECEPTOR}.pdb -L ${LIGAND}.pdb
```
