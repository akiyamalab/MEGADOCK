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
  generate_decoys:
    run: decoygen.cwl
    in:
      ligand: ligands
      dockings: dock_proteins/dockings
    scatter: [ ligand, dockings ]
    scatterMethod: dotproduct
    out: [ decoys ]
  combine_receptor_with_decoy:
    run: combine_receptor_with_decoy.cwl
    in:
      receptor: receptors
      decoy: generate_decoys/decoys
    scatter: [ receptor, decoy ]
    scatterMethod: dotproduct
    out: [ receptor_decoy_complex ]

outputs:
  dockings:
    type: File[]
    outputSource: dock_proteins/dockings
  receptor_decoy_complexes:
    type: File[]
    outputSource: combine_receptor_with_decoy/receptor_decoy_complex
