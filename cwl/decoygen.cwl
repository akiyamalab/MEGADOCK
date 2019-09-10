#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool 

baseCommand: decoygen

hints: 
  DockerRequirement:
    dockerPull: akiyamalab/megadock:4.1.1-cpu

inputs:
  ligand:
    type: File
  dockings:
    type: File
  rank:
    type: int
    default: 1

arguments: # ./decoygen lig.1.pdb lig.pdb dock.out 1
 - $(inputs.ligand.nameroot).$(inputs.rank).pdb
 - $(inputs.ligand)
 - $(inputs.dockings)
 - $(inputs.rank)

outputs:
  decoys:
    type: File
    outputBinding:
      glob: $(inputs.ligand.nameroot).$(inputs.rank).pdb
