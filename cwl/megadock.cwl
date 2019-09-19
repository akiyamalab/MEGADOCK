#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: megadock

hints:
  DockerRequirement:
    dockerPull: akiyamalab/megadock:4.1.1-cpu

inputs:
  receptor:
    type: File
  ligand:
    type: File
  preds_per_rots:
    type: int
    default: 1
  num_preds:
    type: int
    default: 2000

arguments:
 - -R
 - $(inputs.receptor)
 - -L
 - $(inputs.ligand)
 - -t
 - $(inputs.preds_per_rots)
 - -N
 - $(inputs.num_preds)

outputs:
  dockings:
    type: File
    outputBinding:
      glob: $(inputs.receptor.nameroot)-$(inputs.ligand.nameroot).out
