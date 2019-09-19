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

arguments:
 - -R
 - $(inputs.receptor)
 - -L
 - $(inputs.ligand)

outputs:
  dockings:
    type: File
    outputBinding:
      glob: $(inputs.receptor.nameroot)-$(inputs.ligand.nameroot).out
