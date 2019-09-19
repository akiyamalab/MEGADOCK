# MEGADOCK for Common Workflow Language (CWL)

Installation to see the documents from following website:

- https://www.commonwl.org/
- https://github.com/common-workflow-language/cwltool

## Command Example

### Generate docking pose (decoy pdb)

```sh
# megadock -> decoygen (1st rank: default)
cwl-runner cwl/megadock-decoygen-workflow.cwl --ligand data/1gcq_l.pdb --receptor data/1gcq_r.pdb
ls *decoy.*.pdb

# with inputs file in YAML format
cwl-runner cwl/megadock-decoygen-workflow.cwl cwl/inputs.yml
ls *decoy.*.pdb
```

### Calculate docking score (PPI score)

```sh
# megadock -> ppiscore
cwl-runner cwl/megadock-ppiscore-workflow.cwl --ligand data/1gcq_l.pdb --receptor data/1gcq_r.pdb
cat *.score

# with inputs file in YAML format
cwl-runner cwl/megadock-ppiscore-workflow.cwl cwl/inputs.yml
cat *.score
```
