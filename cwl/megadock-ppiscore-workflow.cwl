#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}

inputs:
  ligands: File[]
  receptors: File[]

steps:
  dock_proteins:
    run: megadock.cwl
    in:
      ligand: ligands 
      receptor: receptors
    scatter: [ ligand, receptor ]
    scatterMethod: dotproduct
    out: [ dockings ]
  calc_ppiscore:
    run: ppiscore.cwl
    in:
      dockings: dock_proteins/dockings
    scatter: [ dockings ]
    out: [ ppiscore ]

outputs:
  dockings:
    type: File[]
    outputSource: dock_proteins/dockings
  ppiscore:
    type: File[]
    outputSource: calc_ppiscore/ppiscore
