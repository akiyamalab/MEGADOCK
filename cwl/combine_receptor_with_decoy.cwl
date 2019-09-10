#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool


inputs:
  receptor:
    type: File
  decoy:
    type: File
  rank:
    type: int
    default: 1

baseCommand: cat

arguments:
 - $(inputs.receptor)
 - $(inputs.decoy)

stdout: $(inputs.receptor.nameroot).$(inputs.rank).decoy$(inputs.decoy.basename)

outputs:
  receptor_decoy_complex:
    type: stdout 
