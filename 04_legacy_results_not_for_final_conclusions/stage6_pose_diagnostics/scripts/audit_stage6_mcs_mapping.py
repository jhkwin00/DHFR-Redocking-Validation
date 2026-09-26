#!/usr/bin/env python3
"""
Stage 6.2: native-vs-docked ligand mapping audit.

Purpose
-------
1. Determine whether graph-based MCS covers the whole ligand.
2. Calculate symmetry-aware direct-coordinate RMSD for the best MCS mapping.
3. Calculate an aligned RMSD only as a conformer-shape reference.

Important:
- direct_mcs_rmsd is evaluated in the receptor coordinate frame:
  it retains pose/orientation differences.
- aligned_mcs_rmsd removes rigid-body translation/rotation:
  it is NOT a pose-recovery metric; it only indicates geometric/conformer similarity.
"""

from pathlib import Path
import csv
import itertools
import math
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS

PROJECT = Path(__file__).resolve().parents[3]

NATIVE_DIR = PROJECT / "results" / "stage2_ligand_extraction" / "ligand_only"
DOCKED_DIR = PROJECT / "results" / "stage5_rmsd_qc" / "docked_top1"
OUT_TSV = (
    PROJECT
    / "results"
    / "stage6_pose_diagnostics"
    / "reports"
    / "stage6_mcs_mapping_audit.tsv"
)

PDBIDS = ["1KMV", "3FS6", "3NTZ", "3NXO", "4M6K", "5HSR"]

# Avoid excessive enumeration in highly symmetric molecules.
MAX_MATCHES_PER_MOL = 200


def load_native_pdb(path):
    # proximityBonding=True permits an inferred graph when PDB CONECT records
    # are absent/incomplete. Coordinates remain from the original PDB.
    mol = Chem.MolFromPDBFile(
        str(path),
        sanitize=False,
        removeHs=False,
        proximityBonding=True,
    )
    if mol is None or mol.GetNumConformers() == 0:
        raise ValueError(f"Could not read native ligand PDB: {path}")
    mol = Chem.RemoveHs(mol, sanitize=False)
    # sanitize=False로 PDB를 읽었으므로, MCS의 ring-aware 옵션이
    # 사용할 property cache와 ring information을 명시적으로 초기화한다.
    mol.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(mol)
    return mol


def load_docked_sdf(path):
    supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
    for mol in supplier:
        if mol is not None and mol.GetNumConformers() > 0:
            mol = Chem.RemoveHs(mol, sanitize=False)
            # Docked SDF에도 동일하게 ring information을 준비한다.
            mol.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(mol)
            return mol
    raise ValueError(f"Could not read docked SDF: {path}")


def heavy_count(mol):
    return sum(atom.GetAtomicNum() > 1 for atom in mol.GetAtoms())


def get_coords(mol, atom_indices):
    conf = mol.GetConformer()
    return np.asarray(
        [
            [
                conf.GetAtomPosition(i).x,
                conf.GetAtomPosition(i).y,
                conf.GetAtomPosition(i).z,
            ]
            for i in atom_indices
        ],
        dtype=float,
    )


def raw_rmsd(native_xyz, docked_xyz):
    return float(np.sqrt(np.mean(np.sum((native_xyz - docked_xyz) ** 2, axis=1))))


def kabsch_aligned_rmsd(native_xyz, docked_xyz):
    """
    RMSD after optimal alignment of docked coordinates to native coordinates.
    Used only to show conformer/shape similarity, not receptor-frame pose quality.
    """
    native_centered = native_xyz - native_xyz.mean(axis=0)
    docked_centered = docked_xyz - docked_xyz.mean(axis=0)

    covariance = docked_centered.T @ native_centered
    u, _, vt = np.linalg.svd(covariance)

    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt

    docked_fitted = docked_centered @ rotation
    return raw_rmsd(native_centered, docked_fitted)


def find_mcs(native, docked):
    """
    Attempt exact element/bond-order MCS first.
    If PDB-inferred bond orders prevent matching, relax bond comparison but
    preserve atom-element matching.
    """
    common = dict(
        timeout=30,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
    )

    strict = rdFMCS.FindMCS(
        [native, docked],
        bondCompare=rdFMCS.BondCompare.CompareOrderExact,
        **common,
    )

    if strict.numAtoms >= 3:
        return strict, "exact_bond_order"

    relaxed = rdFMCS.FindMCS(
        [native, docked],
        bondCompare=rdFMCS.BondCompare.CompareAny,
        **common,
    )
    return relaxed, "bond_order_relaxed"


def main():
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

    rows = []

    for pdbid in PDBIDS:
        native_path = NATIVE_DIR / f"{pdbid}_ligand_only.pdb"
        docked_path = DOCKED_DIR / f"{pdbid}_top1_heavy_from_pdbqt.sdf"

        try:
            if not native_path.is_file():
                raise FileNotFoundError(native_path)
            if not docked_path.is_file():
                raise FileNotFoundError(docked_path)

            native = load_native_pdb(native_path)
            docked = load_docked_sdf(docked_path)

            native_heavy = heavy_count(native)
            docked_heavy = heavy_count(docked)

            mcs_result, method = find_mcs(native, docked)
            if mcs_result.numAtoms < 3:
                raise ValueError(
                    f"MCS too small ({mcs_result.numAtoms} atoms); "
                    "native PDB bond inference may be unsuitable."
                )

            query = Chem.MolFromSmarts(mcs_result.smartsString)
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
                raise ValueError("MCS SMARTS did not match both ligand molecules.")

            best_direct = math.inf
            best_aligned = math.inf
            pairs_tested = 0

            for native_match, docked_match in itertools.product(
                native_matches, docked_matches
            ):
                native_xyz = get_coords(native, native_match)
                docked_xyz = get_coords(docked, docked_match)

                direct = raw_rmsd(native_xyz, docked_xyz)
                aligned = kabsch_aligned_rmsd(native_xyz, docked_xyz)

                best_direct = min(best_direct, direct)
                best_aligned = min(best_aligned, aligned)
                pairs_tested += 1

            min_heavy = min(native_heavy, docked_heavy)
            mcs_fraction = mcs_result.numAtoms / min_heavy

            if mcs_fraction >= 0.95:
                interpretation = (
                    "MCS covers essentially the entire ligand; "
                    "MCS is not a tail-excluding core definition."
                )
            else:
                interpretation = (
                    "MCS excludes part of the ligand; candidate for "
                    "tail-excluding core analysis."
                )

            row = {
                "pdbid": pdbid,
                "status": "ok",
                "native_heavy_atom_count": native_heavy,
                "docked_heavy_atom_count": docked_heavy,
                "mcs_atom_count": mcs_result.numAtoms,
                "mcs_fraction_of_smaller_ligand": f"{mcs_fraction:.3f}",
                "mcs_method": method,
                "native_mcs_match_count": len(native_matches),
                "docked_mcs_match_count": len(docked_matches),
                "mapping_pairs_tested": pairs_tested,
                "best_direct_mcs_rmsd_angstrom": f"{best_direct:.3f}",
                "best_aligned_mcs_rmsd_angstrom": f"{best_aligned:.3f}",
                "interpretation": interpretation,
            }
            rows.append(row)

            print(
                f"[OK] {pdbid}: MCS={mcs_result.numAtoms}/{min_heavy} atoms, "
                f"direct={best_direct:.3f} A, aligned={best_aligned:.3f} A"
            )

        except Exception as exc:
            print(f"[ERROR] {pdbid}: {exc}", file=sys.stderr)
            rows.append(
                {
                    "pdbid": pdbid,
                    "status": f"error: {exc}",
                }
            )

    fieldnames = [
        "pdbid",
        "status",
        "native_heavy_atom_count",
        "docked_heavy_atom_count",
        "mcs_atom_count",
        "mcs_fraction_of_smaller_ligand",
        "mcs_method",
        "native_mcs_match_count",
        "docked_mcs_match_count",
        "mapping_pairs_tested",
        "best_direct_mcs_rmsd_angstrom",
        "best_aligned_mcs_rmsd_angstrom",
        "interpretation",
    ]

    with OUT_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    n_ok = sum(row.get("status") == "ok" for row in rows)
    print(f"\nWrote: {OUT_TSV}")
    print(f"Completed successfully: {n_ok}/{len(PDBIDS)} structures")


if __name__ == "__main__":
    main()
