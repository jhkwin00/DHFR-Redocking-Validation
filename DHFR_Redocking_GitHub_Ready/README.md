# DHFR Redocking Validation with AutoDock Vina

## Overview

This project evaluates protein–ligand pose recovery in an AutoDock Vina redocking workflow using six dihydrofolate reductase (DHFR)–ligand crystal structures.

The analysis evaluates not only the top-ranked docking pose, but also the complete set of reported Vina poses for each structure. This makes it possible to distinguish between:

1. **Pose-ranking limitation:** a native-like pose is sampled but is not ranked first.
2. **Pose-sampling limitation under the current setup:** no reported pose satisfies the native-like reference criterion.

During quality control, an inconsistency was identified between legacy converted SDF files and the coordinates in the original Vina PDBQT `MODEL` records. Final conclusions in this repository are therefore based only on a model-by-model validation workflow using the original Vina docking outputs.

> **Scope note:** This is a six-structure redocking case study. Its conclusions apply to the receptor and ligand preparation, docking-box settings, and AutoDock Vina configuration retained in this repository.

## Research Questions

- Can AutoDock Vina recover native-like ligand poses for the selected DHFR crystal structures?
- How frequently does the top-ranked Vina pose satisfy a conventional RMSD-based pose-recovery reference?
- When the top-ranked pose is not native-like, was a native-like pose sampled at a lower rank?
- Are pose identities and coordinates preserved through PDBQT-to-SDF conversion and downstream RMSD analysis?

## Structures Evaluated

| PDB ID | RCSB PDB entry | Reported Vina poses analyzed | Top-1 RMSD (Å) | Best RMSD across reported poses (Å) | Rank of best-RMSD pose | Interpretation |
|---|---|---:|---:|---:|---:|---|
| 1KMV | [RCSB PDB](https://www.rcsb.org/structure/1KMV) | 20 | 0.212 | 0.212 | 1 | Top-1 native-like pose recovered |
| 3FS6 | [RCSB PDB](https://www.rcsb.org/structure/3FS6) | 20 | 1.226 | 1.226 | 1 | Top-1 native-like pose recovered |
| 3NTZ | [RCSB PDB](https://www.rcsb.org/structure/3NTZ) | 20 | 4.294 | 2.118 | 4 | No native-like pose among reported poses under the current docking setup |
| 3NXO | [RCSB PDB](https://www.rcsb.org/structure/3NXO) | 20 | 1.723 | 1.641 | 2 | Top-1 native-like pose recovered |
| 4M6K | [RCSB PDB](https://www.rcsb.org/structure/4M6K) | 20 | 1.036 | 0.989 | 2 | Top-1 native-like pose recovered |
| 5HSR | [RCSB PDB](https://www.rcsb.org/structure/5HSR) | 19 | 2.229 | 0.956 | 7 | Native-like pose sampled but not ranked first |

> A direct RMSD below **2.0 Å** is used here as a conventional descriptive reference for native-like pose recovery. Raw RMSD values, Vina ranks, affinities, and MCS coverage are retained in the final result tables.

## Key Results

- **Top-1 pose recovery:** 4 of 6 structures (**66.7%**) had a top-ranked Vina pose with RMSD < 2.0 Å.
- **Reported-pose ensemble recovery:** 5 of 6 structures (**83.3%**) contained at least one pose with RMSD < 2.0 Å.
- **Ranking limitation observed for 5HSR:** The top-ranked pose had RMSD = **2.229 Å**, whereas the lowest-RMSD pose had RMSD = **0.956 Å** at Vina rank **7**.
- **No pose below the 2.0 Å reference for 3NTZ:** Among 20 analyzed poses, the lowest RMSD was **2.118 Å** at rank **4**.
- **Coordinate-identity validation was necessary:** Legacy top-1 SDF files did not correspond to the original Vina `MODEL 1` coordinates. Final conclusions were recalculated from individually extracted original PDBQT models.

These results show why docking-pose generation and docking-pose ranking should be evaluated separately. In particular, the 5HSR result demonstrates that a native-like pose may be present in the reported docking ensemble while not receiving the best Vina rank.

## Workflow

1. Selected six DHFR–ligand crystal structures from the RCSB Protein Data Bank.
2. Extracted native ligand coordinates from the experimental structures.
3. Prepared receptor and ligand input files for AutoDock Vina.
4. Performed redocking with structure-specific docking-box configurations.
5. Retained original Vina PDBQT outputs, docking configurations, preparation logs, and docking logs.
6. Extracted each Vina PDBQT `MODEL` independently.
7. Converted individual docking models to SDF files for downstream analysis.
8. Calculated receptor-frame, symmetry-aware direct RMSD values between native ligands and docked poses.
9. Recorded Vina rank, affinity, RMSD, MCS coverage, and native-like-pose status for every analyzed pose.
10. Audited the coordinate identity of legacy top-1 files against the original Vina PDBQT model records.
11. Restricted final conclusions to the independently validated Stage 7 outputs.

## Pose-Comparison Metric

Pose similarity was assessed using a **direct, receptor-frame RMSD** metric.

- No spatial fitting or ligand superposition was applied before RMSD calculation.
- Hydrogen atoms were removed before comparison.
- Atom correspondence was determined with an RDKit maximum-common-substructure (MCS) procedure.
- Symmetry-equivalent MCS mappings were enumerated, and the minimum direct RMSD across tested mappings was retained.
- The number of MCS atoms, MCS coverage fraction, and number of mapping pairs tested are reported for each pose.

Therefore, the reported RMSD values represent ligand-coordinate displacement in the receptor coordinate frame rather than RMSD after optimal ligand alignment.

## Validation and Quality Control

Initial analyses used legacy top-1 SDF files generated during an earlier conversion workflow. A subsequent audit showed that these files did not match the coordinates of the original Vina `MODEL 1` records.

To address this issue:

- Original Vina multi-model PDBQT output files were retained as the primary docking record.
- Each PDBQT `MODEL` was extracted independently.
- The Vina affinity associated with each original model was parsed from its corresponding `REMARK VINA RESULT` record.
- Each independently extracted model was converted and evaluated separately.
- Final per-pose RMSD values and per-structure conclusions were generated from this validated model-level workflow.

> **Important:** Files under `04_legacy_results_not_for_final_conclusions/` are retained for provenance and auditability only. They must **not** be used as final redocking conclusions.

## Final Result Files

Use the following files for final interpretation.

| File | Description |
|---|---|
| `01_final_validated_results/stage7_individual_model_summary.tsv` | Per-structure summary of top-1 recovery, best RMSD, rank of the best-RMSD pose, and native-like-pose counts |
| `01_final_validated_results/stage7_individual_model_pose_rmsd.tsv` | Per-pose Vina rank, PDBQT model number, Vina affinity, direct RMSD, MCS coverage, and native-like status |
| `01_final_validated_results/stage5_vs_stage7_rank1_coordinate_audit.tsv` | Coordinate-identity audit comparing legacy Stage 5 files with Stage 7 rank-1 records |
| `01_final_validated_results/stage5_legacy_pose_model_identity.tsv` | Mapping between legacy Stage 5 pose files and original Vina PDBQT model identities |

A field-level definition of all result-table columns is provided in [`docs/data_dictionary.md`](docs/data_dictionary.md).

## Repository Structure

```text
.
├── README.md
├── .gitignore
├── SHA256SUMS.txt
├── FILE_MANIFEST.txt
│
├── 01_final_validated_results/
├── 02_reproducibility_inputs/
├── 03_stage7_scripts_and_pose_files/
├── 04_legacy_results_not_for_final_conclusions/
├── 05_project_documentation/
└── docs/
    └── data_dictionary.md
```

## Reproducibility Materials

This repository retains the primary materials needed to inspect the docking and validation workflow:

- Native ligand coordinate files
- Receptor and ligand PDBQT input files
- AutoDock Vina configuration files
- Original Vina docking outputs
- Docking-preparation and Vina log files
- Model-level pose-extraction and RMSD-analysis scripts
- Converted pose files
- Final result tables
- File manifest and SHA-256 checksums

The archived analysis scripts document the validated analysis logic. Their project-relative paths reflect the original execution layout and may require path updates before re-execution in a new local environment.

## Software Used

- AutoDock Vina
- Meeko
- RDKit
- Open Babel
- Python 3
- NumPy

Exact package and environment versions were not recorded in this handoff package.

## Limitations

- This project is a six-structure redocking case study, not a large-scale docking benchmark.
- Results apply to the current receptor preparation, ligand preparation, docking-box settings, and Vina configuration.
- The 2.0 Å RMSD value is used as a conventional descriptive reference, not as an absolute measure of docking validity.
- RMSD atom correspondence is MCS-based; MCS coverage varies by ligand and is reported in the per-pose results.
- Vina affinity values are docking scores and were not calibrated as experimental binding affinities.
- Molecular-dynamics refinement, free-energy calculations, and machine-learning-based rescoring were **not** performed in this project.

## Data Source and Citation

The experimental protein–ligand structures used in this project were obtained from the [RCSB Protein Data Bank](https://www.rcsb.org/).

When using or discussing these structures, cite the relevant RCSB PDB accession entries (1KMV, 3FS6, 3NTZ, 3NXO, 4M6K, and 5HSR) and their associated primary publications where appropriate.
