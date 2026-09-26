#!/usr/bin/env python3
"""
Stage 6 pose diagnostics, first pass:
  1) Heavy-atom mass-weighted ligand COM distance
  2) Receptor-residue contact overlap (Jaccard similarity)

Input:
  - prepared receptor PDB
  - native ligand-only PDB
  - docked top1 SDF
  - Stage 5 RMSD / visual-QC TSV files

Output:
  results/stage6_pose_diagnostics/reports/stage6_pose_diagnostics.tsv
"""

from pathlib import Path
import csv
import sys

import numpy as np
from rdkit import Chem

# ---------------------------------------------------------------------
# Project paths and analysis settings
# ---------------------------------------------------------------------
PROJECT = Path(__file__).resolve().parents[3]

RECEPTOR_DIR = PROJECT / "results" / "prepared_structures"
NATIVE_DIR = PROJECT / "results" / "stage2_ligand_extraction" / "ligand_only"
DOCKED_DIR = PROJECT / "results" / "stage5_rmsd_qc" / "docked_top1"

STAGE5_RMSD_TSV = (
    PROJECT / "results" / "stage5_rmsd_qc" / "reports" / "stage5_rmsd_summary.tsv"
)
STAGE5_VISUAL_TSV = (
    PROJECT
    / "results"
    / "stage5_rmsd_qc"
    / "reports"
    / "stage5_visual_qc_template.tsv"
)

OUT_TSV = (
    PROJECT
    / "results"
    / "stage6_pose_diagnostics"
    / "reports"
    / "stage6_pose_diagnostics.tsv"
)

PDBIDS = ["1KMV", "3FS6", "3NTZ", "3NXO", "4M6K", "5HSR"]

# A residue is counted as a contact residue if any ligand heavy atom is
# within this distance of any receptor heavy atom.
CONTACT_CUTOFF_A = 4.0

# Standard protein residues only. This prevents any retained ion/cofactor/
# water from being treated as a protein-residue contact.
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL",
    # Common histidine/cysteine naming variants:
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "CYM",
}


def require_file(path: Path, label: str) -> None:
    """Stop with an explicit message instead of silently using bad input."""
    if not path.is_file():
        raise FileNotFoundError(f"{label} file not found:\n  {path}")


def load_pdb(path: Path):
    """
    Load coordinate-bearing PDB without requiring successful bond sanitization.
    For this analysis we require coordinates, atomic numbers and PDB residue info;
    ligand bond orders are not used.
    """
    mol = Chem.MolFromPDBFile(
        str(path),
        sanitize=False,
        removeHs=False,
        proximityBonding=False,
    )
    if mol is None or mol.GetNumConformers() == 0:
        raise ValueError(f"Could not read coordinates from PDB: {path}")
    return mol


def load_sdf(path: Path):
    """Load the first valid molecule from a docked SDF."""
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    for mol in supplier:
        if mol is not None and mol.GetNumConformers() > 0:
            return mol
    raise ValueError(f"Could not read a coordinate-bearing molecule from SDF: {path}")


def heavy_atom_indices(mol):
    """Return indices of non-hydrogen atoms."""
    indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not indices:
        raise ValueError("No heavy atoms were found.")
    return indices


def coordinates_for_indices(mol, indices):
    """Return an (N, 3) numpy coordinate array."""
    conformer = mol.GetConformer()
    return np.array(
        [
            [
                conformer.GetAtomPosition(idx).x,
                conformer.GetAtomPosition(idx).y,
                conformer.GetAtomPosition(idx).z,
            ]
            for idx in indices
        ],
        dtype=float,
    )


def heavy_atom_mass_weighted_com(mol):
    """
    Mass-weighted COM calculated with heavy atoms only.
    Hydrogen positions/absence in PDB/PDBQT-derived files therefore cannot
    distort the comparison.
    """
    indices = heavy_atom_indices(mol)
    coords = coordinates_for_indices(mol, indices)
    masses = np.array(
        [mol.GetAtomWithIdx(idx).GetMass() for idx in indices],
        dtype=float,
    )
    return np.average(coords, axis=0, weights=masses)


def residue_label(atom):
    """
    Make a stable residue identifier:
      CHAIN:RESNUM[:INSERTION_CODE]:RESNAME
    e.g. A:31:ILE or A:31:A:ILE
    """
    info = atom.GetPDBResidueInfo()
    if info is None:
        return None

    residue_name = info.GetResidueName().strip().upper()
    if residue_name not in STANDARD_AA:
        return None

    chain = info.GetChainId().strip() or "_"
    residue_number = info.GetResidueNumber()
    insertion_code = info.GetInsertionCode().strip()

    if insertion_code:
        return f"{chain}:{residue_number}:{insertion_code}:{residue_name}"
    return f"{chain}:{residue_number}:{residue_name}"


def receptor_heavy_atoms_by_residue(receptor_mol):
    """
    Build:
      {residue_label: Nx3 coordinates of receptor heavy atoms}

    Only standard amino-acid residues are retained.
    """
    grouped = {}

    for atom in receptor_mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue

        label = residue_label(atom)
        if label is None:
            continue

        pos = receptor_mol.GetConformer().GetAtomPosition(atom.GetIdx())
        grouped.setdefault(label, []).append([pos.x, pos.y, pos.z])

    if not grouped:
        raise ValueError(
            "No standard amino-acid receptor residues were parsed. "
            "Check receptor PDB formatting / residue names."
        )

    return {
        label: np.asarray(coords, dtype=float)
        for label, coords in grouped.items()
    }


def contacting_residues(receptor_coords_by_residue, ligand_mol, cutoff_a):
    """
    Return receptor residues with at least one receptor-heavy-atom to
    ligand-heavy-atom distance <= cutoff_a.
    """
    ligand_coords = coordinates_for_indices(ligand_mol, heavy_atom_indices(ligand_mol))
    cutoff_sq = cutoff_a * cutoff_a
    contacts = set()

    for residue, receptor_coords in receptor_coords_by_residue.items():
        displacement = receptor_coords[:, None, :] - ligand_coords[None, :, :]
        squared_distances = np.sum(displacement * displacement, axis=2)

        if np.any(squared_distances <= cutoff_sq):
            contacts.add(residue)

    return contacts


def read_tsv_by_pdbid(path: Path):
    """Read a TSV into {pdbid: row_dictionary}."""
    require_file(path, "TSV input")
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["pdbid"].strip(): row for row in reader}


def fmt_set(values):
    """Semicolon-separated, reproducible residue-list formatting."""
    return ";".join(sorted(values))


def main():
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

    stage5_rmsd = read_tsv_by_pdbid(STAGE5_RMSD_TSV)
    stage5_visual = read_tsv_by_pdbid(STAGE5_VISUAL_TSV)

    output_rows = []

    for pdbid in PDBIDS:
        receptor_path = RECEPTOR_DIR / f"{pdbid}_receptor.pdb"
        native_path = NATIVE_DIR / f"{pdbid}_ligand_only.pdb"
        docked_path = DOCKED_DIR / f"{pdbid}_top1_heavy_from_pdbqt.sdf"

        try:
            require_file(receptor_path, f"{pdbid} receptor")
            require_file(native_path, f"{pdbid} native ligand")
            require_file(docked_path, f"{pdbid} docked top1 ligand")

            receptor = load_pdb(receptor_path)
            native = load_pdb(native_path)
            docked = load_sdf(docked_path)

            native_com = heavy_atom_mass_weighted_com(native)
            docked_com = heavy_atom_mass_weighted_com(docked)
            com_distance = float(np.linalg.norm(native_com - docked_com))

            receptor_by_residue = receptor_heavy_atoms_by_residue(receptor)
            native_contacts = contacting_residues(
                receptor_by_residue, native, CONTACT_CUTOFF_A
            )
            docked_contacts = contacting_residues(
                receptor_by_residue, docked, CONTACT_CUTOFF_A
            )

            intersection = native_contacts & docked_contacts
            union = native_contacts | docked_contacts
            jaccard = (len(intersection) / len(union)) if union else float("nan")

            rmsd_row = stage5_rmsd.get(pdbid, {})
            visual_row = stage5_visual.get(pdbid, {})

            output_rows.append(
                {
                    "pdbid": pdbid,
                    "status": "ok",
                    "contact_cutoff_angstrom": f"{CONTACT_CUTOFF_A:.1f}",
                    "full_heavy_atom_rmsd_angstrom": rmsd_row.get(
                        "rmsd_heavy_atom_angstrom", ""
                    ),
                    "vina_affinity_kcal_mol": rmsd_row.get(
                        "vina_affinity_kcal_mol", ""
                    ),
                    "same_binding_pocket_visual": visual_row.get(
                        "same_binding_pocket", ""
                    ),
                    "visual_class": visual_row.get("visual_class", ""),
                    "ligand_com_distance_angstrom": f"{com_distance:.3f}",
                    "native_heavy_atom_count": len(heavy_atom_indices(native)),
                    "docked_heavy_atom_count": len(heavy_atom_indices(docked)),
                    "native_contact_residue_count": len(native_contacts),
                    "docked_contact_residue_count": len(docked_contacts),
                    "contact_intersection_count": len(intersection),
                    "contact_union_count": len(union),
                    "contact_jaccard": f"{jaccard:.3f}",
                    "native_contact_residues": fmt_set(native_contacts),
                    "docked_contact_residues": fmt_set(docked_contacts),
                    "shared_contact_residues": fmt_set(intersection),
                    "native_only_contact_residues": fmt_set(
                        native_contacts - docked_contacts
                    ),
                    "docked_only_contact_residues": fmt_set(
                        docked_contacts - native_contacts
                    ),
                    "notes": (
                        "Heavy-atom mass-weighted COM; standard amino-acid "
                        f"residue contacts at <= {CONTACT_CUTOFF_A:.1f} A."
                    ),
                }
            )

            print(
                f"[OK] {pdbid}: COM={com_distance:.3f} A, "
                f"Jaccard={jaccard:.3f}, "
                f"contacts native/docked={len(native_contacts)}/{len(docked_contacts)}"
            )

        except Exception as exc:
            print(f"[ERROR] {pdbid}: {exc}", file=sys.stderr)
            output_rows.append(
                {
                    "pdbid": pdbid,
                    "status": f"error: {exc}",
                    "contact_cutoff_angstrom": f"{CONTACT_CUTOFF_A:.1f}",
                    "full_heavy_atom_rmsd_angstrom": stage5_rmsd.get(
                        pdbid, {}
                    ).get("rmsd_heavy_atom_angstrom", ""),
                    "vina_affinity_kcal_mol": stage5_rmsd.get(
                        pdbid, {}
                    ).get("vina_affinity_kcal_mol", ""),
                }
            )

    fieldnames = [
        "pdbid",
        "status",
        "contact_cutoff_angstrom",
        "full_heavy_atom_rmsd_angstrom",
        "vina_affinity_kcal_mol",
        "same_binding_pocket_visual",
        "visual_class",
        "ligand_com_distance_angstrom",
        "native_heavy_atom_count",
        "docked_heavy_atom_count",
        "native_contact_residue_count",
        "docked_contact_residue_count",
        "contact_intersection_count",
        "contact_union_count",
        "contact_jaccard",
        "native_contact_residues",
        "docked_contact_residues",
        "shared_contact_residues",
        "native_only_contact_residues",
        "docked_only_contact_residues",
        "notes",
    ]

    with OUT_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    n_ok = sum(row["status"] == "ok" for row in output_rows)
    print(f"\nWrote: {OUT_TSV}")
    print(f"Completed successfully: {n_ok}/{len(PDBIDS)} structures")


if __name__ == "__main__":
    main()
