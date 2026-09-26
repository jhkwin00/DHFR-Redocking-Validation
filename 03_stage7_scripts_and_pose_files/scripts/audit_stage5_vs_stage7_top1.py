#!/usr/bin/env python3
"""
Compare legacy Stage 5 top-1 SDF against rank 1 in the Stage 7 all-pose SDF.

This audit determines whether the earlier Stage 5 conversion used the same
Vina MODEL 1 geometry that Stage 7 uses.
"""

from pathlib import Path
import csv
import itertools
import math

import numpy as np
from rdkit import Chem


PROJECT = Path(__file__).resolve().parents[3]

STAGE5_DIR = PROJECT / "results" / "stage5_rmsd_qc" / "docked_top1"
STAGE7_DIR = (
    PROJECT
    / "results"
    / "stage7_topN_pose_analysis"
    / "converted_sdf"
)
OUT_TSV = (
    PROJECT
    / "results"
    / "stage7_topN_pose_analysis"
    / "reports"
    / "stage5_vs_stage7_rank1_coordinate_audit.tsv"
)

PDBIDS = ["1KMV", "3FS6", "3NTZ", "3NXO", "4M6K", "5HSR"]


def load_first_sdf_molecule(path):
    supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)

    for mol in supplier:
        if mol is not None and mol.GetNumConformers() > 0:
            mol = Chem.RemoveHs(mol, sanitize=False)
            return mol

    raise ValueError(f"No readable coordinate-bearing SDF molecule: {path}")


def heavy_atom_indices(mol):
    return [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() > 1
    ]


def coordinates(mol, indices):
    conf = mol.GetConformer()
    return np.asarray(
        [
            [
                conf.GetAtomPosition(i).x,
                conf.GetAtomPosition(i).y,
                conf.GetAtomPosition(i).z,
            ]
            for i in indices
        ],
        dtype=float,
    )


def rmsd(xyz_a, xyz_b):
    return float(np.sqrt(np.mean(np.sum((xyz_a - xyz_b) ** 2, axis=1))))


def best_element_matched_rmsd(stage5, stage7):
    """
    For each chemical element, enumerate coordinate assignments independently.
    This distinguishes atom-order changes from true coordinate changes.

    The ligands are small; a safety cap remains in place for large symmetry
    permutations.
    """
    idx5 = heavy_atom_indices(stage5)
    idx7 = heavy_atom_indices(stage7)

    if len(idx5) != len(idx7):
        raise ValueError(
            f"Heavy-atom count differs: stage5={len(idx5)}, stage7={len(idx7)}"
        )

    groups5 = {}
    groups7 = {}

    for idx in idx5:
        atomic_num = stage5.GetAtomWithIdx(idx).GetAtomicNum()
        groups5.setdefault(atomic_num, []).append(idx)

    for idx in idx7:
        atomic_num = stage7.GetAtomWithIdx(idx).GetAtomicNum()
        groups7.setdefault(atomic_num, []).append(idx)

    if set(groups5) != set(groups7):
        raise ValueError("Element compositions differ between the two SDFs.")

    # Avoid combinatorial explosion; coordinate order RMSD remains available.
    if any(len(indices) > 8 for indices in groups5.values()):
        return float("nan"), "not_enumerated_element_group_too_large"

    stage5_order = []
    all_permutations = []

    for atomic_num in sorted(groups5):
        group5 = groups5[atomic_num]
        group7 = groups7[atomic_num]

        if len(group5) != len(group7):
            raise ValueError(
                f"Element {atomic_num} count differs between the two SDFs."
            )

        stage5_order.extend(group5)
        all_permutations.append(list(itertools.permutations(group7)))

    xyz5 = coordinates(stage5, stage5_order)

    best = math.inf
    for selected_groups in itertools.product(*all_permutations):
        stage7_order = [
            idx
            for group_permutation in selected_groups
            for idx in group_permutation
        ]
        xyz7 = coordinates(stage7, stage7_order)
        best = min(best, rmsd(xyz5, xyz7))

    return best, "enumerated_by_element"


def main():
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    rows = []

    for pdbid in PDBIDS:
        stage5_path = STAGE5_DIR / f"{pdbid}_top1_heavy_from_pdbqt.sdf"
        stage7_path = STAGE7_DIR / f"{pdbid}_all_vina_poses.sdf"

        try:
            stage5 = load_first_sdf_molecule(stage5_path)
            stage7_rank1 = load_first_sdf_molecule(stage7_path)

            idx5 = heavy_atom_indices(stage5)
            idx7 = heavy_atom_indices(stage7_rank1)

            if len(idx5) != len(idx7):
                raise ValueError(
                    f"Heavy-atom count differs: stage5={len(idx5)}, "
                    f"stage7={len(idx7)}"
                )

            atom_order_rmsd = rmsd(
                coordinates(stage5, idx5),
                coordinates(stage7_rank1, idx7),
            )
            element_rmsd, method = best_element_matched_rmsd(
                stage5,
                stage7_rank1,
            )

            if atom_order_rmsd < 1e-4:
                conclusion = "identical_coordinates_and_atom_order"
            elif element_rmsd < 1e-4:
                conclusion = "same_coordinate_set_but_atom_order_differs"
            else:
                conclusion = "different_pose_coordinates_or_conversion_result"

            rows.append(
                {
                    "pdbid": pdbid,
                    "status": "ok",
                    "stage5_heavy_atom_count": len(idx5),
                    "stage7_rank1_heavy_atom_count": len(idx7),
                    "same_atom_order_coordinate_rmsd_A": f"{atom_order_rmsd:.6f}",
                    "best_element_matched_coordinate_rmsd_A": (
                        f"{element_rmsd:.6f}"
                        if not math.isnan(element_rmsd)
                        else ""
                    ),
                    "element_matching_method": method,
                    "conclusion": conclusion,
                }
            )

            print(
                f"[OK] {pdbid}: order-RMSD={atom_order_rmsd:.6f} A; "
                f"element-matched={element_rmsd:.6f} A; {conclusion}"
            )

        except Exception as exc:
            print(f"[ERROR] {pdbid}: {exc}")
            rows.append(
                {
                    "pdbid": pdbid,
                    "status": f"error: {exc}",
                }
            )

    fields = [
        "pdbid",
        "status",
        "stage5_heavy_atom_count",
        "stage7_rank1_heavy_atom_count",
        "same_atom_order_coordinate_rmsd_A",
        "best_element_matched_coordinate_rmsd_A",
        "element_matching_method",
        "conclusion",
    ]

    with OUT_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote: {OUT_TSV}")


if __name__ == "__main__":
    main()
