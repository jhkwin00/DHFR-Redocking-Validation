# Data Dictionary

This document defines the columns in the final validated result tables under `01_final_validated_results/`.

## General conventions

- Distances and RMSD values are in **angstroms (Å)**.
- `vina_rank` denotes the rank reported by AutoDock Vina; rank 1 is the top-ranked pose.
- `native_like_at_2A = yes` means `direct_native_rmsd_angstrom < 2.0`.
- “Native” refers to the ligand coordinates extracted from the corresponding experimental PDB structure.
- The final redocking conclusions are based on the two Stage 7 result files. The two Stage 5 files are coordinate-identity audit records, not final pose-recovery results.

## `stage7_individual_model_summary.tsv`

One row per PDB structure. This is the main structure-level summary.

| Column | Description |
|---|---|
| `pdbid` | Four-character RCSB PDB identifier for the DHFR–ligand structure. |
| `status` | Processing status recorded by the validated workflow. |
| `total_models_analyzed` | Number of individual Vina PDBQT `MODEL` records analyzed for the structure. |
| `top1_model_number` | Original PDBQT model number for the Vina rank-1 pose. |
| `top1_affinity_kcal_mol` | Vina affinity reported for the rank-1 pose, in kcal/mol. This is a docking score, not an experimental binding affinity. |
| `top1_native_rmsd_angstrom` | Direct receptor-frame RMSD between the rank-1 pose and the native ligand. |
| `best_native_rmsd_among_models_angstrom` | Lowest direct receptor-frame RMSD among all analyzed Vina models for the structure. |
| `rank_of_best_native_rmsd` | Vina rank of the pose with the lowest native RMSD. |
| `model_of_best_native_rmsd` | Original PDBQT model number of the pose with the lowest native RMSD. |
| `affinity_of_best_native_rmsd_kcal_mol` | Vina affinity of the pose with the lowest native RMSD, in kcal/mol. |
| `native_like_model_count_at_2A` | Number of analyzed poses with direct native RMSD below 2.0 Å. |
| `first_native_like_rank_at_2A` | Lowest Vina rank at which a pose with direct native RMSD below 2.0 Å occurred. |
| `interpretation` | Structure-level result label produced by the workflow. |

## `stage7_individual_model_pose_rmsd.tsv`

One row per individually extracted Vina pose. This is the main pose-level result table.

| Column | Description |
|---|---|
| `pdbid` | Four-character RCSB PDB identifier. |
| `vina_rank` | Rank assigned by AutoDock Vina; 1 is the top-ranked pose. |
| `pdbqt_model_number` | Corresponding `MODEL` number in the original multi-model Vina PDBQT output. |
| `vina_affinity_kcal_mol` | Vina affinity for the pose, in kcal/mol. It is a docking score, not an experimentally calibrated binding affinity. |
| `direct_native_rmsd_angstrom` | Direct receptor-frame heavy-atom RMSD between the docked pose and the native ligand, without spatial fitting. |
| `native_like_at_2A` | `yes` when `direct_native_rmsd_angstrom < 2.0`; otherwise `no`. The threshold is used as a conventional descriptive reference. |
| `mcs_atom_count` | Number of heavy atoms in the maximum common substructure used to define atom correspondence between docked and native ligands. |
| `mcs_fraction_of_smaller_ligand` | MCS atom count divided by the heavy-atom count of the smaller of the two compared ligand graphs. |
| `mapping_pairs_tested` | Number of symmetry-equivalent atom-mapping pairs evaluated before retaining the minimum direct RMSD. |
| `legacy_stage5_coordinate_rmsd_angstrom` | Coordinate difference between the Stage 7 pose and its associated legacy Stage 5 record. This is an audit field, not a native-pose RMSD. |
| `status` | Processing status recorded by the validated workflow. |

## `stage5_vs_stage7_rank1_coordinate_audit.tsv`

One row per PDB structure. This file audits whether the legacy Stage 5 top-1 coordinate file corresponds to the validated Stage 7 rank-1 record.

| Column | Description |
|---|---|
| `pdbid` | Four-character RCSB PDB identifier. |
| `status` | Processing status recorded by the audit workflow. |
| `stage5_heavy_atom_count` | Heavy-atom count in the legacy Stage 5 top-1 file. |
| `stage7_rank1_heavy_atom_count` | Heavy-atom count in the validated Stage 7 rank-1 pose. |
| `same_atom_order_coordinate_rmsd_A` | Coordinate RMSD computed by comparing atoms in the stored order of the two files. This is an identity audit, not a native-pose RMSD. |
| `best_element_matched_coordinate_rmsd_A` | Best coordinate RMSD after element-matched mapping, where such enumeration was performed. This is an identity audit, not a native-pose RMSD. |
| `element_matching_method` | Method or reason recorded for element-matched coordinate comparison. |
| `conclusion` | Audit conclusion concerning coordinate agreement between legacy Stage 5 and validated Stage 7 rank-1 files. |

## `stage5_legacy_pose_model_identity.tsv`

One row per PDB structure. This file identifies the original Vina model closest to each legacy Stage 5 pose file.

| Column | Description |
|---|---|
| `pdbid` | Four-character RCSB PDB identifier. |
| `status` | Processing status recorded by the audit workflow. |
| `closest_vina_rank_to_legacy_stage5_file` | Vina rank of the original model with coordinates closest to the legacy Stage 5 file. |
| `closest_pdbqt_model_number` | PDBQT `MODEL` number of the original model closest to the legacy Stage 5 file. |
| `legacy_to_closest_model_rmsd_angstrom` | Coordinate RMSD between the legacy Stage 5 file and the closest original Vina model. This is an identity-audit metric, not a native-pose RMSD. |
| `legacy_identity_conclusion` | Audit conclusion concerning whether the legacy Stage 5 file exactly matches an original Vina model. |
