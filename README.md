# MEGADOCK

**MEGADOCK** is an ultra-high-performance protein-protein prediction software for heterogeneous supercomputers using FFT-grid-based docking with MPI/OpenMP/GPU parallelization.

[![License: CC BY-NC 4.0](https://licensebuttons.net/l/by-nc/4.0/80x15.png)](https://creativecommons.org/licenses/by-nc/4.0/) 
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

![Build Status](https://github.com/akiyamalab/MEGADOCK/workflows/Build%20containers/badge.svg?branch=master)


## Target Environments

| Type | Target Env.     | Approach           |
|:----:|-----------------|--------------------|
|  (a) | GPU cluster     | GPU + OpenMP + MPI |
|  (b) | CPU cluster     | OpenMP + MPI       |
|  (c) | GPU node        | GPU + OpenMP       |
|  (d) | CPU node        | OpenMP             |


## Installation and Command Details

For installation and command details, please read appropriate section on followings:
- Read command and script details
  - [doc/README.md](./doc/README.md)
- Build a binary from source code
  - [doc/BUILD.md](./doc/BUILD.md)
- Build a docker container image
  - [doc/README_for_docker.md](doc/README_for_docker.md)
    - [akiyamalab/megadock](https://hub.docker.com/r/akiyamalab/megadock/) (Docker Hub)


## Reference

Masahito Ohue, Takehiro Shimoda, Shuji Suzuki, Yuri Matsuzaki, Takashi Ishida, Yutaka Akiyama. **MEGADOCK 4.0: an ultra-high-performance protein-protein docking software for heterogeneous supercomputers**, *Bioinformatics*, 30(22): 3281-3283, 2014. http://doi.org/10.1093/bioinformatics/btu532

Masahito Ohue, Yuri Matsuzaki, Nobuyuki Uchikoga, Takashi Ishida, Yutaka Akiyama. **MEGADOCK: An all-to-all protein-protein interaction prediction system using tertiary structure data**, *Protein and Peptide Letters*, 21(8): 766-778, 2014. https://doi.org/10.2174/09298665113209990050


## Older Versions

[http://www.bi.cs.titech.ac.jp/megadock/archives/](http://www.bi.cs.titech.ac.jp/megadock/archives/)


## License

MEGADOCK is licensed by CC BY-NC 4.0. (See [LICENSE](./LICENSE))
This software and derivatives are NOT allowed for any commercial use without formal prior authorization by Tokyo Institute of Technology. Please contact us about commercial use. If you are considering commercial use, please contact us as there are different charged use options licensed by the Tokyo Institute of Technology.


## Fundings

This work is partially supported by JSPS Grant-in-Aid for Scientific Research (KAKENHI) (A) Grant Number 24240044.

----
Copyright © 2014-2022 Akiyama Laboratory, Tokyo Institute of Technology, All Rights Reserved.
