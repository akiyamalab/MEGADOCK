# MEGADOCK with Docker

Docker container images are available on [`akiyamalab/megadock`](https://hub.docker.com/r/akiyamalab/megadock/) repository (DockerHub).
- [`akiyamalab/megadock:cpu`](https://hub.docker.com/r/akiyamalab/megadock/)
- [`akiyamalab/megadock:gpu`](https://hub.docker.com/r/akiyamalab/megadock/)

## Requirements
| Requirement Tools                                         | GPU | CPU | Notes       |
|:----------------------------------------------------------|:---:|:---:|:------------|
| [Docker](https://docs.docker.com/engine/installation/)    | x   | x   |             |
| [nvidia-docker](https://github.com/NVIDIA/nvidia-docker)  | x   |     | > 2.0 |

## Quick Example
```sh
# CPU single node (OpenMP parallelization)
docker run akiyamalab/megadock:cpu megadock -R data/1gcq_r.pdb -L data/1gcq_l.pdb -o data/1gcq_r-1gcq_r.out

# GPU single node (GPU parallelization)
docker run --runtime=nvidia akiyamalab/megadock:gpu megadock-gpu -R data/1gcq_r.pdb -L data/1gcq_l.pdb -o data/1gcq_r-1gcq_r.out
```

----

## Build Docker Container for MEGADOCK (CPU)

### 1. build Docker image
```sh
# on ${REPOSITORY_ROOT} dir
docker build . -f Dockerfiles/cpu/Dockerfile -t akiyamalab/megadock:cpu
```

### 2. run sample
```sh
docker run akiyamalab/megadock:cpu megadock -R data/1gcq_r.pdb -L data/1gcq_l.pdb

# optional) start interactive shell
docker run -it akiyamalab/megadock:cpu

# optional) run with your pdb (e.g. ${DATA_PATH} = your pdb-data directory abs path  )
docker run -v ${DATA_PATH}:/opt/MEGADOCK/data akiyamalab/megadock:cpu megadock -R data/${RECEPTOR}.pdb -L data/${LIGAND}.pdb
```

## Build Docker Container for MEGADOCK (GPU)

**[nvidia-docker](https://github.com/NVIDIA/nvidia-docker) is required.**

### 1. build Docker image
```sh
# on ${REPOSITORY_ROOT} dir
docker build . -f Dockerfiles/gpu/Dockerfile -t akiyamalab/megadock:gpu
```

### 2. run sample
```sh
docker run -it --runtime=nvidia akiyamalab/megadock:gpu megadock-gpu -R data/1gcq_r.pdb -L data/1gcq_l.pdb

# optional) start interactive shell
docker run -it --runtime=nvidia akiyamalab/megadock:gpu

# optional) run with your pdb (e.g. ${DATA_PATH} = your pdb-data directory abs path  )
docker run -it --runtime=nvidia -v ${DATA_PATH}:/opt/MEGADOCK/data akiyamalab/megadock:gpu megadock-gpu -R data/${RECEPTOR}.pdb -L data/${LIGAND}.pdb
```
