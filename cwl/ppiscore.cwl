#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool 

baseCommand: ppiscore

hints: 
  DockerRequirement:
    dockerPull: akiyamalab/megadock:4.1.1-ppiscore

inputs:
  dockings:
    type: File

arguments: # ppiscore dock.out
  - $(inputs.dockings)

outputs:
  ppiscore:
    type: stdout
stdout: $(inputs.dockings.nameroot).score
