#!/usr/bin/env python3
"""
Validate Vina pose order and RMSD by processing each PDBQT MODEL independently.

This avoids assuming that record order from a multi-model PDBQT -> multi-record
SDF conversion necessarily equals Vina MODEL order.

For each MODEL:
  - Extract its original PDBQT block.
  - Retain and parse its own REMARK VINA RESULT affinity.
  - Convert that individual model to SDF.
  - Compute native-vs-pose direct, symmetry-aware heavy-atom RMSD.
  - Compare with the legacy Stage 5 "top1" SDF to identify which MODEL it is.

Outputs:
  stage7_individual_model_pose_rmsd.tsv
  stage7_individual_model_summary.tsv
  stage5_legacy_pose_model_identity.tsv
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
STAGE5_DIR = PROJECT / "results" / "stage5_rmsd_qc" / "docked_top1"

WORK_DIR = PROJECT / "results" / "stage7_topN_pose_analysis" / "individual_models"
REPORT_DIR = PROJECT / "results" / "stage7_topN_pose_analysis" / "reports"

POSE_TSV = REPORT_DIR / "stage7_individual_model_pose_rmsd.tsv"
SUMMARY_TSV = REPORT_DIR / "stage7_individual_model_summary.tsv"
LEGACY_TSV = REPORT_DIR / "stage5_legacy_pose_model_identity.tsv"

PDBIDS = ["1KMV", "3FS6", "3NTZ", "3NXO", "4M6K", "5HSR"]

NATIVE_LIKE_THRESHOLD_A = 2.0
MAX_MATCHES_PER_MOL = 200

AFFINITY_RE = re.compile(
    r"^REMARK VINA RESULT:\s*(-?\d+(?:\.\d+)?)"
)


def prepare_mol(mol, label):
    mol = Chem.RemoveHs(mol, sanitize=False)

    if mol is None or mol.GetNumConformers() == 0:
        raise ValueError(f"{label}: no usable coordinates after H removal.")

    mol.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(mol)
    return mol


def load_native(path):
    mol = Chem.MolFromPDBFile(
        str(path),
        sanitize=False,
        removeHs=False,
        proximityBonding=True,
    )
    if mol is None:
        raise ValueError(f"Cannot read native ligand: {path}")
    return prepare_mol(mol, "Native ligand")


def load_first_sdf(path, label):
    supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)

    for mol in supplier:
        if mol is not None and mol.GetNumConformers() > 0:
            return prepare_mol(mol, label)

    raise ValueError(f"Cannot read coordinate-bearing SDF: {path}")


def heavy_atom_count(mol):
    return sum(atom.GetAtomicNum() > 1 for atom in mol.GetAtoms())


def coordinates(mol, atom_indices):
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


def raw_rmsd(xyz_a, xyz_b):
    return float(np.sqrt(np.mean(np.sum((xyz_a - xyz_b) ** 2, axis=1))))


def symmetry_aware_direct_rmsd(reference, mobile):
    """
    RMSD in the receptor coordinate frame. No spatial fitting is performed.
    Symmetry-equivalent MCS matches are enumerated.
    """
    mcs = rdFMCS.FindMCS(
        [reference, mobile],
        timeout=30,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrderExact,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
    )

    if mcs.numAtoms < 3:
        raise ValueError(f"MCS has only {mcs.numAtoms} atoms.")

    query = Chem.MolFromSmarts(mcs.smartsString)
    if query is None:
        raise ValueError("Could not build MCS SMARTS query.")

    ref_matches = reference.GetSubstructMatches(
        query, uniquify=True, maxMatches=MAX_MATCHES_PER_MOL
    )
    mob_matches = mobile.GetSubstructMatches(
        query, uniquify=True, maxMatches=MAX_MATCHES_PER_MOL
    )

    if not ref_matches or not mob_matches:
        raise ValueError("MCS query did not match both molecules.")

    best = math.inf
    pairs_tested = 0

    for ref_match, mob_match in itertools.product(ref_matches, mob_matches):
        rmsd = raw_rmsd(
            coordinates(reference, ref_match),
            coordinates(mobile, mob_match),
        )
        best = min(best, rmsd)
        pairs_tested += 1

    fraction = mcs.numAtoms / min(
        heavy_atom_count(reference),
        heavy_atom_count(mobile),
    )

    return best, mcs.numAtoms, fraction, pairs_tested


def extract_model_blocks(pdbqt_path):
    """
    Return a list of (model_number, affinity, complete_pdbqt_text).

    Any header text before MODEL is included in every extracted file. Vina
    model-specific REMARK VINA RESULT lines stay associated with that model.
    """
    lines = pdbqt_path.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines(keepends=True)

    header = []
    models = []

    in_model = False
    current = []
    model_number = None

    for line in lines:
        if line.startswith("MODEL"):
            if in_model:
                raise ValueError(
                    f"Encountered nested MODEL in {pdbqt_path}: {line.strip()}"
                )

            in_model = True
            parts = line.split()
            model_number = int(parts[1]) if len(parts) > 1 else len(models) + 1
            current = list(header)
            current.append(line)

        elif line.startswith("ENDMDL"):
            if not in_model:
                raise ValueError(
                    f"ENDMDL found before MODEL in {pdbqt_path}"
                )

            current.append(line)
            model_text = "".join(current)

            affinity = None
            for model_line in current:
                match = AFFINITY_RE.match(model_line)
                if match:
                    affinity = float(match.group(1))
                    break

            if affinity is None:
                raise ValueError(
                    f"MODEL {model_number} has no REMARK VINA RESULT affinity."
                )

            models.append((model_number, affinity, model_text))
            in_model = False
            current = []
            model_number = None

        elif in_model:
            current.append(line)
        else:
            header.append(line)

    if in_model:
        raise ValueError(f"Unclosed final MODEL in {pdbqt_path}")

    if not models:
        raise ValueError(f"No MODEL blocks found in {pdbqt_path}")

    return models


def convert_one_pdbqt_to_sdf(pdbqt_path, sdf_path):
    obabel = shutil.which("obabel")
    if obabel is None:
        raise RuntimeError("obabel not found in PATH.")

    command = [
        obabel,
        "-ipdbqt",
        str(pdbqt_path),
        "-osdf",
        "-O",
        str(sdf_path),
    ]

    result = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"Open Babel failed for {pdbqt_path}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )

    if not sdf_path.is_file() or sdf_path.stat().st_size == 0:
        raise RuntimeError(f"No SDF generated: {sdf_path}")


def write_tsv(path, fieldnames, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    if shutil.which("obabel") is None:
        raise SystemExit("ERROR: obabel is not available in this environment.")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    pose_rows = []
    summary_rows = []
    legacy_rows = []

    for pdbid in PDBIDS:
        try:
            native_path = NATIVE_DIR / f"{pdbid}_ligand_only.pdb"
            vina_path = VINA_DIR / f"{pdbid}_docked.pdbqt"
            legacy_path = STAGE5_DIR / f"{pdbid}_top1_heavy_from_pdbqt.sdf"

            native = load_native(native_path)
            legacy = load_first_sdf(legacy_path, "Legacy Stage 5 pose")
            model_blocks = extract_model_blocks(vina_path)

            pdb_work_dir = WORK_DIR / pdbid
            pdb_work_dir.mkdir(parents=True, exist_ok=True)

            pdb_pose_rows = []
            legacy_comparisons = []

            for sequential_rank, (model_number, affinity, text) in enumerate(
                model_blocks, start=1
            ):
                model_pdbqt = pdb_work_dir / f"{pdbid}_model_{model_number:02d}.pdbqt"
                model_sdf = pdb_work_dir / f"{pdbid}_model_{model_number:02d}.sdf"

                model_pdbqt.write_text(text, encoding="utf-8")
                convert_one_pdbqt_to_sdf(model_pdbqt, model_sdf)
                pose = load_first_sdf(
                    model_sdf,
                    f"{pdbid} MODEL {model_number}",
                )

                native_rmsd, mcs_atoms, mcs_fraction, mapping_pairs = (
                    symmetry_aware_direct_rmsd(native, pose)
                )
                legacy_rmsd, _, _, _ = symmetry_aware_direct_rmsd(legacy, pose)

                row = {
                    "pdbid": pdbid,
                    "vina_rank": sequential_rank,
                    "pdbqt_model_number": model_number,
                    "vina_affinity_kcal_mol": f"{affinity:.3f}",
                    "direct_native_rmsd_angstrom": f"{native_rmsd:.3f}",
                    "native_like_at_2A": (
                        "yes" if native_rmsd < NATIVE_LIKE_THRESHOLD_A else "no"
                    ),
                    "mcs_atom_count": mcs_atoms,
                    "mcs_fraction_of_smaller_ligand": f"{mcs_fraction:.3f}",
                    "mapping_pairs_tested": mapping_pairs,
                    "legacy_stage5_coordinate_rmsd_angstrom": f"{legacy_rmsd:.6f}",
                    "status": "ok",
                }

                pdb_pose_rows.append(row)
                pose_rows.append(row)
                legacy_comparisons.append((legacy_rmsd, sequential_rank, model_number))

            top1 = pdb_pose_rows[0]
            best_native = min(
                pdb_pose_rows,
                key=lambda row: float(row["direct_native_rmsd_angstrom"]),
            )
            native_like = [
                row for row in pdb_pose_rows if row["native_like_at_2A"] == "yes"
            ]

            best_legacy_rmsd, legacy_rank, legacy_model = min(legacy_comparisons)

            if native_like:
                first_native_rank = min(int(row["vina_rank"]) for row in native_like)
                if first_native_rank == 1:
                    interpretation = "top1_native_like_pose_recovered"
                else:
                    interpretation = "native_like_pose_sampled_but_ranked_below_top1"
            else:
                first_native_rank = ""
                interpretation = (
                    "no_native_like_pose_among_reported_vina_poses_"
                    "under_current_docking_setup"
                )

            summary_rows.append(
                {
                    "pdbid": pdbid,
                    "status": "ok",
                    "total_models_analyzed": len(pdb_pose_rows),
                    "top1_model_number": top1["pdbqt_model_number"],
                    "top1_affinity_kcal_mol": top1["vina_affinity_kcal_mol"],
                    "top1_native_rmsd_angstrom": top1[
                        "direct_native_rmsd_angstrom"
                    ],
                    "best_native_rmsd_among_models_angstrom": best_native[
                        "direct_native_rmsd_angstrom"
                    ],
                    "rank_of_best_native_rmsd": best_native["vina_rank"],
                    "model_of_best_native_rmsd": best_native["pdbqt_model_number"],
                    "affinity_of_best_native_rmsd_kcal_mol": best_native[
                        "vina_affinity_kcal_mol"
                    ],
                    "native_like_model_count_at_2A": len(native_like),
                    "first_native_like_rank_at_2A": first_native_rank,
                    "interpretation": interpretation,
                }
            )

            legacy_identity = (
                "legacy_stage5_file_matches_this_model"
                if best_legacy_rmsd < 0.001
                else "legacy_stage5_file_does_not_exactly_match_any_model"
            )

            legacy_rows.append(
                {
                    "pdbid": pdbid,
                    "status": "ok",
                    "closest_vina_rank_to_legacy_stage5_file": legacy_rank,
                    "closest_pdbqt_model_number": legacy_model,
                    "legacy_to_closest_model_rmsd_angstrom": f"{best_legacy_rmsd:.6f}",
                    "legacy_identity_conclusion": legacy_identity,
                }
            )

            print(
                f"[OK] {pdbid}: top1={top1['direct_native_rmsd_angstrom']} A; "
                f"best={best_native['direct_native_rmsd_angstrom']} A "
                f"(rank {best_native['vina_rank']}); "
                f"legacy closest=model {legacy_model}, "
                f"RMSD={best_legacy_rmsd:.6f} A"
            )

        except Exception as exc:
            print(f"[ERROR] {pdbid}: {exc}", file=sys.stderr)
            summary_rows.append({"pdbid": pdbid, "status": f"error: {exc}"})
            legacy_rows.append({"pdbid": pdbid, "status": f"error: {exc}"})

    pose_fields = [
        "pdbid",
        "vina_rank",
        "pdbqt_model_number",
        "vina_affinity_kcal_mol",
        "direct_native_rmsd_angstrom",
        "native_like_at_2A",
        "mcs_atom_count",
        "mcs_fraction_of_smaller_ligand",
        "mapping_pairs_tested",
        "legacy_stage5_coordinate_rmsd_angstrom",
        "status",
    ]

    summary_fields = [
        "pdbid",
        "status",
        "total_models_analyzed",
        "top1_model_number",
        "top1_affinity_kcal_mol",
        "top1_native_rmsd_angstrom",
        "best_native_rmsd_among_models_angstrom",
        "rank_of_best_native_rmsd",
        "model_of_best_native_rmsd",
        "affinity_of_best_native_rmsd_kcal_mol",
        "native_like_model_count_at_2A",
        "first_native_like_rank_at_2A",
        "interpretation",
    ]

    legacy_fields = [
        "pdbid",
        "status",
        "closest_vina_rank_to_legacy_stage5_file",
        "closest_pdbqt_model_number",
        "legacy_to_closest_model_rmsd_angstrom",
        "legacy_identity_conclusion",
    ]

    write_tsv(POSE_TSV, pose_fields, pose_rows)
    write_tsv(SUMMARY_TSV, summary_fields, summary_rows)
    write_tsv(LEGACY_TSV, legacy_fields, legacy_rows)

    n_ok = sum(row.get("status") == "ok" for row in summary_rows)
    print(f"\nWrote: {POSE_TSV}")
    print(f"Wrote: {SUMMARY_TSV}")
    print(f"Wrote: {LEGACY_TSV}")
    print(f"Completed successfully: {n_ok}/{len(PDBIDS)} structures")


if __name__ == "__main__":
    main()
