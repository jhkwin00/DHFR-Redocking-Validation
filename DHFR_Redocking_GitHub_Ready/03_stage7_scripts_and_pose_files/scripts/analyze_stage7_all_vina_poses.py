#!/usr/bin/env python3
"""
Stage 7: Analyze all Vina output poses rather than top-1 alone.

For every Vina MODEL:
  1. Convert multi-model PDBQT output to a multi-record SDF using Open Babel.
  2. Calculate a symmetry-aware, receptor-frame RMSD to the native ligand.
  3. Preserve Vina rank and affinity.
  4. Report the best RMSD observed among all generated poses.

Interpretation:
- Native-like pose sampled at a non-top1 rank:
    scoring/ranking failure candidate.
- No native-like pose among all output poses:
    native-like pose was not sampled under the current docking setup.

The 2.0 A threshold is retained only as a conventional descriptive reference;
all raw RMSD values and ranks are also written.
"""

from pathlib import Path
import csv
import itertools
import math
import re
import shutil
import subprocess
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS


PROJECT = Path(__file__).resolve().parents[3]

NATIVE_DIR = PROJECT / "results" / "stage2_ligand_extraction" / "ligand_only"
VINA_DIR = PROJECT / "results" / "stage4_vina_redocking" / "docking_out"

CONVERTED_DIR = PROJECT / "results" / "stage7_topN_pose_analysis" / "converted_sdf"
REPORT_DIR = PROJECT / "results" / "stage7_topN_pose_analysis" / "reports"

POSE_TSV = REPORT_DIR / "stage7_pose_rmsd_by_rank.tsv"
SUMMARY_TSV = REPORT_DIR / "stage7_topN_summary.tsv"

PDBIDS = ["1KMV", "3FS6", "3NTZ", "3NXO", "4M6K", "5HSR"]

# Conventional pose-recovery reference threshold.
NATIVE_LIKE_RMSD_THRESHOLD_A = 2.0

# Protect against pathological symmetry enumeration.
MAX_MATCHES_PER_MOL = 200

VINA_AFFINITY_PATTERN = re.compile(
    r"^REMARK VINA RESULT:\s*"
    r"(-?\d+(?:\.\d+)?)"
)


def prepare_mol(mol, label):
    """Remove hydrogens and initialize RDKit property/ring information."""
    mol = Chem.RemoveHs(mol, sanitize=False)

    if mol is None or mol.GetNumConformers() == 0:
        raise ValueError(f"{label}: no usable conformer after H removal.")

    mol.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(mol)
    return mol


def load_native_ligand(path):
    mol = Chem.MolFromPDBFile(
        str(path),
        sanitize=False,
        removeHs=False,
        proximityBonding=True,
    )
    if mol is None:
        raise ValueError(f"Cannot read native ligand PDB: {path}")
    return prepare_mol(mol, "Native ligand")


def load_sdf_poses(path):
    poses = []

    supplier = Chem.SDMolSupplier(
        str(path),
        sanitize=False,
        removeHs=False,
    )

    for index, mol in enumerate(supplier, start=1):
        if mol is None:
            raise ValueError(
                f"Open Babel SDF record {index} could not be parsed: {path}"
            )
        poses.append(prepare_mol(mol, f"Docked pose {index}"))

    if not poses:
        raise ValueError(f"No poses could be read from converted SDF: {path}")

    return poses


def heavy_atom_count(mol):
    return sum(atom.GetAtomicNum() > 1 for atom in mol.GetAtoms())


def coordinates(mol, atom_indices):
    conformer = mol.GetConformer()
    return np.asarray(
        [
            [
                conformer.GetAtomPosition(atom_index).x,
                conformer.GetAtomPosition(atom_index).y,
                conformer.GetAtomPosition(atom_index).z,
            ]
            for atom_index in atom_indices
        ],
        dtype=float,
    )


def raw_rmsd(native_xyz, docked_xyz):
    """RMSD without any alignment; receptor-frame pose metric."""
    return float(
        np.sqrt(np.mean(np.sum((native_xyz - docked_xyz) ** 2, axis=1)))
    )


def find_mcs_query(native, docked):
    """
    Match by element and exact bond order. The Stage 6 audit showed MCS covers
    all heavy atoms in the top-1 poses for this project.
    """
    result = rdFMCS.FindMCS(
        [native, docked],
        timeout=30,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrderExact,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
    )

    if result.numAtoms < 3:
        raise ValueError(f"MCS too small: {result.numAtoms} atoms.")

    query = Chem.MolFromSmarts(result.smartsString)
    if query is None:
        raise ValueError("Unable to construct MCS SMARTS query.")

    return result, query


def symmetry_aware_direct_rmsd(native, docked):
    """
    Calculate lowest direct RMSD over symmetry-equivalent MCS mappings.
    Coordinates remain in the receptor frame: no fitting/alignment is done.
    """
    mcs_result, query = find_mcs_query(native, docked)

    native_matches = native.GetSubstructMatches(
        query,
        uniquify=True,
        maxMatches=MAX_MATCHES_PER_MOL,
    )
    docked_matches = docked.GetSubstructMatches(
        query,
        uniquify=True,
        maxMatches=MAX_MATCHES_PER_MOL,
    )

    if not native_matches or not docked_matches:
        raise ValueError("MCS query did not match native and docked ligand.")

    best_rmsd = math.inf
    pairs_tested = 0

    for native_match, docked_match in itertools.product(
        native_matches, docked_matches
    ):
        native_xyz = coordinates(native, native_match)
        docked_xyz = coordinates(docked, docked_match)

        best_rmsd = min(best_rmsd, raw_rmsd(native_xyz, docked_xyz))
        pairs_tested += 1

    native_heavy = heavy_atom_count(native)
    docked_heavy = heavy_atom_count(docked)
    smaller_heavy_count = min(native_heavy, docked_heavy)
    mcs_fraction = mcs_result.numAtoms / smaller_heavy_count

    return {
        "rmsd": best_rmsd,
        "mcs_atoms": mcs_result.numAtoms,
        "mcs_fraction": mcs_fraction,
        "mapping_pairs_tested": pairs_tested,
        "native_heavy": native_heavy,
        "docked_heavy": docked_heavy,
    }


def parse_vina_affinities(pdbqt_path):
    """
    Read Vina scores in MODEL order from:
      REMARK VINA RESULT:  -10.5  0.000  0.000
    """
    affinities = []

    with pdbqt_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = VINA_AFFINITY_PATTERN.match(line)
            if match:
                affinities.append(float(match.group(1)))

    if not affinities:
        raise ValueError(f"No 'REMARK VINA RESULT' scores in: {pdbqt_path}")

    return affinities


def convert_pdbqt_to_sdf(pdbqt_path, sdf_path):
    """Convert all PDBQT models to a multi-record SDF with Open Babel."""
    obabel = shutil.which("obabel")
    if obabel is None:
        raise RuntimeError(
            "Open Babel executable 'obabel' was not found in PATH. "
            "Run: which obabel"
        )

    command = [
        obabel,
        "-ipdbqt",
        str(pdbqt_path),
        "-osdf",
        "-O",
        str(sdf_path),
    ]

    completed = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "Open Babel conversion failed.\n"
            f"Command: {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )

    if not sdf_path.is_file() or sdf_path.stat().st_size == 0:
        raise RuntimeError(f"Open Babel produced no SDF: {sdf_path}")


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    CONVERTED_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    all_pose_rows = []
    summary_rows = []

    for pdbid in PDBIDS:
        native_path = NATIVE_DIR / f"{pdbid}_ligand_only.pdb"
        docked_pdbqt_path = VINA_DIR / f"{pdbid}_docked.pdbqt"
        converted_sdf_path = CONVERTED_DIR / f"{pdbid}_all_vina_poses.sdf"

        try:
            if not native_path.is_file():
                raise FileNotFoundError(f"Missing native ligand: {native_path}")
            if not docked_pdbqt_path.is_file():
                raise FileNotFoundError(
                    f"Missing Vina output PDBQT: {docked_pdbqt_path}"
                )

            native = load_native_ligand(native_path)
            affinities = parse_vina_affinities(docked_pdbqt_path)
            convert_pdbqt_to_sdf(docked_pdbqt_path, converted_sdf_path)
            poses = load_sdf_poses(converted_sdf_path)

            if len(poses) != len(affinities):
                raise ValueError(
                    f"Pose count mismatch: PDBQT contains {len(affinities)} Vina "
                    f"scores but converted SDF contains {len(poses)} molecules."
                )

            pdb_rows = []

            for rank, (affinity, docked) in enumerate(
                zip(affinities, poses),
                start=1,
            ):
                result = symmetry_aware_direct_rmsd(native, docked)

                # An MCS smaller than the molecule would mean the graphs were
                # not fully comparable; preserve that audit information.
                full_graph_match = result["mcs_fraction"] >= 0.95

                row = {
                    "pdbid": pdbid,
                    "vina_rank": rank,
                    "vina_affinity_kcal_mol": f"{affinity:.3f}",
                    "direct_heavy_atom_rmsd_angstrom": f"{result['rmsd']:.3f}",
                    "native_like_at_2A": (
                        "yes"
                        if result["rmsd"] < NATIVE_LIKE_RMSD_THRESHOLD_A
                        else "no"
                    ),
                    "native_heavy_atom_count": result["native_heavy"],
                    "docked_heavy_atom_count": result["docked_heavy"],
                    "mcs_atom_count": result["mcs_atoms"],
                    "mcs_fraction_of_smaller_ligand": f"{result['mcs_fraction']:.3f}",
                    "full_graph_match": "yes" if full_graph_match else "no",
                    "mapping_pairs_tested": result["mapping_pairs_tested"],
                    "status": "ok",
                }

                pdb_rows.append(row)
                all_pose_rows.append(row)

            best_row = min(
                pdb_rows,
                key=lambda row: float(row["direct_heavy_atom_rmsd_angstrom"]),
            )
            top1_row = next(row for row in pdb_rows if row["vina_rank"] == 1)

            best_rmsd = float(best_row["direct_heavy_atom_rmsd_angstrom"])
            top1_rmsd = float(top1_row["direct_heavy_atom_rmsd_angstrom"])
            best_rank = int(best_row["vina_rank"])

            native_like_rows = [
                row for row in pdb_rows if row["native_like_at_2A"] == "yes"
            ]

            if native_like_rows:
                first_native_like_rank = min(
                    int(row["vina_rank"]) for row in native_like_rows
                )

                if first_native_like_rank == 1:
                    interpretation = "top1_native_like_pose_recovered"
                else:
                    interpretation = "native_like_pose_sampled_but_ranked_below_top1"
            else:
                first_native_like_rank = ""
                interpretation = (
                    "no_native_like_pose_among_reported_vina_poses_"
                    "under_current_docking_setup"
                )

            summary_rows.append(
                {
                    "pdbid": pdbid,
                    "status": "ok",
                    "total_vina_poses_analyzed": len(pdb_rows),
                    "top1_affinity_kcal_mol": top1_row[
                        "vina_affinity_kcal_mol"
                    ],
                    "top1_rmsd_angstrom": f"{top1_rmsd:.3f}",
                    "best_rmsd_among_all_poses_angstrom": f"{best_rmsd:.3f}",
                    "rank_of_best_rmsd_pose": best_rank,
                    "affinity_of_best_rmsd_pose_kcal_mol": best_row[
                        "vina_affinity_kcal_mol"
                    ],
                    "native_like_pose_count_at_2A": len(native_like_rows),
                    "first_native_like_rank_at_2A": first_native_like_rank,
                    "rmsd_improvement_top1_minus_best_angstrom": (
                        f"{top1_rmsd - best_rmsd:.3f}"
                    ),
                    "interpretation": interpretation,
                }
            )

            print(
                f"[OK] {pdbid}: poses={len(pdb_rows)}, "
                f"top1 RMSD={top1_rmsd:.3f} A, "
                f"best RMSD={best_rmsd:.3f} A (rank {best_rank})"
            )

        except Exception as exc:
            print(f"[ERROR] {pdbid}: {exc}", file=sys.stderr)

            summary_rows.append(
                {
                    "pdbid": pdbid,
                    "status": f"error: {exc}",
                }
            )

    pose_fields = [
        "pdbid",
        "vina_rank",
        "vina_affinity_kcal_mol",
        "direct_heavy_atom_rmsd_angstrom",
        "native_like_at_2A",
        "native_heavy_atom_count",
        "docked_heavy_atom_count",
        "mcs_atom_count",
        "mcs_fraction_of_smaller_ligand",
        "full_graph_match",
        "mapping_pairs_tested",
        "status",
    ]

    summary_fields = [
        "pdbid",
        "status",
        "total_vina_poses_analyzed",
        "top1_affinity_kcal_mol",
        "top1_rmsd_angstrom",
        "best_rmsd_among_all_poses_angstrom",
        "rank_of_best_rmsd_pose",
        "affinity_of_best_rmsd_pose_kcal_mol",
        "native_like_pose_count_at_2A",
        "first_native_like_rank_at_2A",
        "rmsd_improvement_top1_minus_best_angstrom",
        "interpretation",
    ]

    write_tsv(POSE_TSV, pose_fields, all_pose_rows)
    write_tsv(SUMMARY_TSV, summary_fields, summary_rows)

    n_ok = sum(row.get("status") == "ok" for row in summary_rows)

    print(f"\nWrote pose-level table: {POSE_TSV}")
    print(f"Wrote summary table:    {SUMMARY_TSV}")
    print(f"Completed successfully: {n_ok}/{len(PDBIDS)} structures")


if __name__ == "__main__":
    main()
