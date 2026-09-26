#!/usr/bin/env python3
"""
Stage 6.3: Bemis-Murcko scaffold (core) pose diagnostics.

For each native/docked ligand pair:
  - Extract the native ligand's Bemis-Murcko scaffold.
  - Match that scaffold in native and docked ligand structures.
  - Enumerate symmetry-equivalent matches.
  - Calculate:
      * direct scaffold RMSD in receptor coordinates (pose metric)
      * aligned scaffold RMSD (shape/conformer reference only)

Important:
  direct_scaffold_rmsd_angstrom retains translation/orientation differences
  in the receptor frame and is the relevant pose-recovery metric.

  aligned_scaffold_rmsd_angstrom removes rigid-body translation/rotation.
  It is NOT a docking-pose score.
"""

from pathlib import Path
import csv
import itertools
import math
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


PROJECT = Path(__file__).resolve().parents[3]

NATIVE_DIR = PROJECT / "results" / "stage2_ligand_extraction" / "ligand_only"
DOCKED_DIR = PROJECT / "results" / "stage5_rmsd_qc" / "docked_top1"

OUT_TSV = (
    PROJECT
    / "results"
    / "stage6_pose_diagnostics"
    / "reports"
    / "stage6_murcko_core_rmsd.tsv"
)

PDBIDS = ["1KMV", "3FS6", "3NTZ", "3NXO", "4M6K", "5HSR"]
MAX_MATCHES_PER_MOL = 200


def prepare_mol(mol):
    """Remove H atoms and initialize caches/ring info required by RDKit."""
    mol = Chem.RemoveHs(mol, sanitize=False)
    if mol is None or mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no usable coordinates after H removal.")

    mol.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(mol)
    return mol


def load_native_pdb(path):
    mol = Chem.MolFromPDBFile(
        str(path),
        sanitize=False,
        removeHs=False,
        proximityBonding=True,
    )
    if mol is None:
        raise ValueError(f"Could not read native ligand PDB: {path}")
    return prepare_mol(mol)


def load_docked_sdf(path):
    supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
    for mol in supplier:
        if mol is not None and mol.GetNumConformers() > 0:
            return prepare_mol(mol)
    raise ValueError(f"Could not read docked ligand SDF: {path}")


def heavy_atom_count(mol):
    return sum(atom.GetAtomicNum() > 1 for atom in mol.GetAtoms())


def coords(mol, atom_indices):
    conformer = mol.GetConformer()
    return np.asarray(
        [
            [
                conformer.GetAtomPosition(idx).x,
                conformer.GetAtomPosition(idx).y,
                conformer.GetAtomPosition(idx).z,
            ]
            for idx in atom_indices
        ],
        dtype=float,
    )


def direct_rmsd(native_xyz, docked_xyz):
    return float(
        np.sqrt(np.mean(np.sum((native_xyz - docked_xyz) ** 2, axis=1)))
    )


def kabsch_aligned_rmsd(native_xyz, docked_xyz):
    """
    RMSD after optimally rotating/translating docked coordinates to native.
    Only used as a conformer/shape reference.
    """
    native_centered = native_xyz - native_xyz.mean(axis=0)
    docked_centered = docked_xyz - docked_xyz.mean(axis=0)

    covariance = docked_centered.T @ native_centered
    u, _, vt = np.linalg.svd(covariance)
    rotation = u @ vt

    # Enforce a proper rotation, not a reflection.
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt

    docked_fitted = docked_centered @ rotation
    return direct_rmsd(native_centered, docked_fitted)


def scaffold_smiles(scaffold):
    """Coordinate-independent, reproducible representation for reporting."""
    return Chem.MolToSmiles(scaffold, canonical=True)


def main():
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    rows = []

    for pdbid in PDBIDS:
        native_path = NATIVE_DIR / f"{pdbid}_ligand_only.pdb"
        docked_path = DOCKED_DIR / f"{pdbid}_top1_heavy_from_pdbqt.sdf"

        try:
            if not native_path.is_file():
                raise FileNotFoundError(f"Missing native ligand: {native_path}")
            if not docked_path.is_file():
                raise FileNotFoundError(f"Missing docked ligand: {docked_path}")

            native = load_native_pdb(native_path)
            docked = load_docked_sdf(docked_path)

            native_heavy = heavy_atom_count(native)
            docked_heavy = heavy_atom_count(docked)

            scaffold = MurckoScaffold.GetScaffoldForMol(native)
            scaffold.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(scaffold)

            scaffold_atoms = heavy_atom_count(scaffold)
            if scaffold_atoms < 3:
                raise ValueError(
                    f"Murcko scaffold has only {scaffold_atoms} heavy atoms; "
                    "a core RMSD is not meaningful."
                )

            native_matches = native.GetSubstructMatches(
                scaffold,
                uniquify=True,
                maxMatches=MAX_MATCHES_PER_MOL,
            )
            docked_matches = docked.GetSubstructMatches(
                scaffold,
                uniquify=True,
                maxMatches=MAX_MATCHES_PER_MOL,
            )

            if not native_matches:
                raise ValueError("Native ligand did not match its own scaffold.")
            if not docked_matches:
                raise ValueError(
                    "Docked ligand did not match the native Murcko scaffold."
                )

            best_direct = math.inf
            best_aligned = math.inf
            pair_count = 0

            for native_match, docked_match in itertools.product(
                native_matches, docked_matches
            ):
                native_xyz = coords(native, native_match)
                docked_xyz = coords(docked, docked_match)

                best_direct = min(
                    best_direct,
                    direct_rmsd(native_xyz, docked_xyz),
                )
                best_aligned = min(
                    best_aligned,
                    kabsch_aligned_rmsd(native_xyz, docked_xyz),
                )
                pair_count += 1

            fraction_native = scaffold_atoms / native_heavy
            fraction_docked = scaffold_atoms / docked_heavy

            if fraction_native >= 0.95:
                interpretation = (
                    "Scaffold covers nearly the entire ligand; "
                    "tail exclusion is limited."
                )
            else:
                interpretation = (
                    "Scaffold excludes terminal substituent(s); "
                    "direct scaffold RMSD is a tail-reduced pose metric."
                )

            row = {
                "pdbid": pdbid,
                "status": "ok",
                "native_heavy_atom_count": native_heavy,
                "docked_heavy_atom_count": docked_heavy,
                "murcko_scaffold_heavy_atom_count": scaffold_atoms,
                "scaffold_fraction_native": f"{fraction_native:.3f}",
                "scaffold_fraction_docked": f"{fraction_docked:.3f}",
                "native_scaffold_match_count": len(native_matches),
                "docked_scaffold_match_count": len(docked_matches),
                "mapping_pairs_tested": pair_count,
                "direct_scaffold_rmsd_angstrom": f"{best_direct:.3f}",
                "aligned_scaffold_rmsd_angstrom": f"{best_aligned:.3f}",
                "murcko_scaffold_smiles": scaffold_smiles(scaffold),
                "interpretation": interpretation,
            }
            rows.append(row)

            print(
                f"[OK] {pdbid}: scaffold={scaffold_atoms}/{native_heavy} heavy atoms, "
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
        "murcko_scaffold_heavy_atom_count",
        "scaffold_fraction_native",
        "scaffold_fraction_docked",
        "native_scaffold_match_count",
        "docked_scaffold_match_count",
        "mapping_pairs_tested",
        "direct_scaffold_rmsd_angstrom",
        "aligned_scaffold_rmsd_angstrom",
        "murcko_scaffold_smiles",
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
