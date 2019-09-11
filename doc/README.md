# Command Guide

This documentation show the command usage and examples for protein-protein docking calculation using MEGADOCK.
If you do not have any MEGADOCK binary, please read the build documentation ([doc/BUILD.md](./BUILD.md)).

> Note: The sample input files (`.pdb`, `.table`) are stored in `data` directory.


## Section Link

Please select your target environment and check the command example.

| Type | Target Env.     | Binary             | Note 
|:----:|-----------------|--------------------| :-
|  (a) | [GPU cluster](#command-for-multiple-nodes-environment) | `megadock-gpu-dp`  | multiple nodes
|  (b) | [CPU cluster](#command-for-multiple-nodes-environment) | `megadock-dp`      | multiple nodes
|  (c) | [GPU node](#command-for-single-node-environment)      | `megadock-gpu`     | single node
|  (d) | [CPU node](#command-for-single-node-environment)      | `megadock`         | single node

Other Commands:

| Command    | Note
| :-         | :-
| `decoygen` | [Generate decoy pdb files from docking output](#Generate-decoy-from-docking-output)
| `block`    | [Blocking target residues](#Blocking-target-residues)
| `ppiscore` | [Protein-Protein interaction prediction](#Protein-Protein-interaction-prediction)

----

# Command for Single Node Environment

Target Environments:
  - (c) Compile for GPU node (GPU)
  - (d) Compile for CPU node (only thread parallelization)

## Command example

```sh
# (c) GPU single node
./megadock-gpu -R receptor.pdb -L ligand.pdb -o receptor-ligand.out

# (d) CPU single node
./megadock -R receptor.pdb -L ligand.pdb -o receptor-ligand.out
```

### Run small sample

```sh
# (c) GPU single node
./megadock-gpu -R data/1gcq_r.pdb -L data/1gcq_l.pdb -o data/1gcq_r-1gcq_l.out

# (d) CPU single node
./megadock -R data/1gcq_r.pdb -L data/1gcq_l.pdb -o data/1gcq_r-1gcq_l.out
```

## Parameters for docking calculation

| Required            | Note
| :-------------------| :-
| `-R [filename]`     | receptor pdb file (input)
| `-L [filename]`     | ligand pdb file (input)


| Optional                  | default       | Note
| :------------------------ | :------------ | :-- 
| `-o [filename]`           | $R-$L.out     | output filename (output)
| `-O`                      |               | output docking detail files
| `-N [integer]`            | 2000          | set the number of output predictions
| `-t [integer]`            | 1             | set the number of predictions per each rotation
| `-F [integer]`            |               | set the number of FFT point (default: none)
| `-v [float]`              | 1.2 (angstrom)| set the voxel pitch size
| `-D`                      |               | set the 6 deg. (54000 angles) of rotational sampling (default to none, 15 deg. (3600 angles) of rotational sampling)
| `-r [integer]`            | 3600          | set the number of rotational sampling angles (54000: 54000 angles, 1: 1 angles, 24: 24 angles, default to 3600 angles)
| `-e [float]`              | 1.0           | set the electrostatics term ratio
| `-d [float]`              | 1.0           | set the hydrophobic term ratio
| `-a [float]`              | -45.0         | set the rPSC receptor core penalty
| `-b [float]`              | 1.0           | set the rPSC ligand core penalty
| `-f [1/2/3]`              | 3             | set function (1: rPSC only, 2: rPSC+Elec, 3: rPSC+Elec+RDE, default to 3)
| `-h`                      |               | show command help


----

# Command for Multiple Nodes Environment

Target Environments:
- (a) Compile for GPU cluster (GPU & MPI)
- (b) Compile for CPU cluster (MPI)

## Command example

```sh
# (a) GPU & MPI (e.g. 4 MPI processes)
mpirun -n 4 ./megadock-gpu-dp -tb data/SAMPLE.table

# (b) CPU & MPI (e.g. 4 MPI processes)
mpirun -n 4 ./megadock-dp -tb data/SAMPLE.table
```


## Parameters for parallel execution

| Required                     | Note
| :----------------------------| :-
| `-tb [filename] `            | docking table (.table)


| Optional                  | default       | Note
| :------------------------ | :------------ | :-- 
| `-lg [filename]`          | master.log    | log filename
| `-rt [integer]`           | 0             | the number of retries


## Example of docking table (.table)

Parameters and docking target list should be written in a text file.
A table file should have `TITLE` and `PARAM` lines followed by parameters listed by the same order as in the `PARAM` line.
`SAMPLE.table` file in this package shows an example.

- [data/SAMPLE.table](data/SAMPLE.table)

```
TITLE=sample jobs
PARAM=-R $1 -L $2 -O
data/1gcq_r.pdb	data/1gcq_r.pdb
data/1gcq_r.pdb	data/1gcq_l.pdb
data/1gcq_l.pdb	data/1gcq_r.pdb
data/1gcq_l.pdb	data/1gcq_l.pdb
```

You can specify docking parameters loaded by MEGADOCK by setting `PARAM` in the table file.
Lines follows the `PARAM` lines specifies parameters for each docking job which will be distributed to available nodes by MPI.

### More details about parameters for each docking calculation

- [See above section](#Parametersfor-docking-calculation)

----

# Thread parallelization (OpenMP)

MEGADOCK parallelizes rotation calculations by using OpenMP.
You can specify the number of OpenMP threads for parallel calculations by environmental variable such as `$OMP_NUM_THREADS`.

```sh
# e.g.) megadock binary using 8 OpenMP threads
export OMP_NUM_THREADS=8
./megadock -R data/receptor.pdb -L data/ligand.pdb -o data/receptor_ligand.out

# e.g.) each MPI process uses 8 OpenMP threads
mpirun -n 4 -x OMP_NUM_THREADS=8  ./megadock-dp -tb data/SAMPLE.table
```

----


# Generate decoy from docking output

Each docking job generates docking output file in which rotation angles (x, y, z) of ligand from the initial structure, numbers of voxel translation (x, y, z) from receptor center coordinate and docking scores are listed.  

If you want to generate decoy pdbfiles, please use `decoygen` command. It is generated together with MEGADOCK binary.

## Command Usage

```sh
./decoygen [decoy_filename] [used_ligand.pdb] [.outfile] [decoy no.]
```


## Command Example

For example, the megadock output file was generated as `dock.out` using `rec.pdb` as receputer input and `lig.pdb` as ligand output, and you want to **generate 1st ranked decoy**.

```sh
# megadock docking calculation
megadock -R rec.pdb -L lig.pdb -o dock.out

# generate 1st ranked decoy
./decoygen lig.1.pdb lig.pdb dock.out 1
cat rec.pdb lig.1.pdb > decoy.1.pdb

# lig.1.pdb   : rotated and translated lig.pdb
# decoy.1.pdb : complex pdb file
```

If you want to generate all decoys;

```sh
# megadock docking calculation
megadock -R rec.pdb -L lig.pdb -o dock.out

# generate all decoys
for i in `seq 1 2000`; do ./decoygen lig.${i}.pdb lig.pdb dock.out $i; cat rec.pdb lig.${i}.pdb > decoy.${i}.pdb; done
```

----

# Blocking target residues

`block` tool (python script) is to block some residues of a receptor protein from docking.
If you know that some residues are not in the binding sides, please list their residue numbers and input `block` tool.
This program changes the name of residues to "BLK" and prints the new pdb on the screen.

```sh
./block [pdbfile] [chain] [target residue list]
# ex)
# ./block 1gcq_r.pdb B 182-186,189,195-198,204 > blocked.pdb
```

Target residues list is separated by commas and no spaces.
You can also use `-` (hyphen): "182-186" means blocking residues of 182, 183, ..., 186.
Blocked residues are substituted for 'BLK'.
Updated PDB coordinates are written to "standard output". MEGADOCK can only block receptor residues.

----

# Protein-Protein interaction prediction

You can calculates the evaluated value of Protein-Protein Interaction scores (PPI score) by using `ppiscore` command which is written in perl script.

## Command Example

```sh
# 1. Docking calculation using MEGADOCK with `-t=3` option
megadock -R rec.pdb -L lig.pdb -o dock.out -t 3 -N 10800

# -t: the number of predictions per each rotation
# -N: the number of output predictions


# 2. Output PPI score
./ppiscore dock.out 10800
```

### Note

`ppiscore` is compatible with re-ranking tools to get better docking candidates (e.g. [ZRANK1](http://zdock.umassmed.edu/software/)) and those tools can be used for the output files of docking calculation using megadock.
If you want to use such reranking tools, please refer to [our website documentation](http://www.bi.cs.titech.ac.jp/megadock/ppi.html).

For more details about `ppiscore`, please check [our webpage documentation](http://www.bi.cs.titech.ac.jp/megadock/ppi.html) or following references.

----


## References

- Masahito Ohue, Takehiro Shimoda, Shuji Suzuki, Yuri Matsuzaki, Takashi Ishida, Yutaka Akiyama. MEGADOCK 4.0: an ultra-high-performance protein-protein docking software for heterogeneous supercomputers. Bioinformatics, 30(22), 3281-3283, 2014.

- Masahito Ohue, Yuri Matsuzaki, Nobuyuki Uchikoga, Takashi Ishida, Yutaka Akiyama. MEGADOCK: An all-to-all protein-protein interaction prediction system using tertiary structure data. Protein and Peptide Letters, 21(8), 766-778, 2014.

- Takehiro Shimoda, Takashi Ishida, Shuji Suzuki, Masahito Ohue, Yutaka Akiyama. MEGADOCK-GPU: An accelerated protein-protein docking calculation on GPUs. In Proc. ACM-BCB 2013 (ParBio Workshop 2013), 884-890, 2013.
 
- Yuri Matsuzaki, Nobuyuki Uchikoga, Masahito Ohue, Takehiro Shimoda, Toshiyuki Sato, Takashi Ishida, Yutaka Akiyama. MEGADOCK 3.0: A high-performance protein-protein interaction prediction software using hybrid parallel computing for petascale supercomputing environments. Source Code for Biology and Medicine, 8(1): 18, 2013.

- Masahito Ohue, Yuri Matsuzaki, Takashi Ishida, Yutaka Akiyama. Improvement of the Protein-Protein Docking Prediction by Introducing a Simple Hydrophobic Interaction Model: an Application to Interaction Pathway Analysis. Lecture Note in Bioinformatics 7632 (In Proc. of PRIB 2012), 178-187, Springer Heidelberg, 2012.

----

## Contact

- Email : megadock@bi.c.titech.ac.jp
- URL : http://www.bi.cs.titech.ac.jp/megadock
