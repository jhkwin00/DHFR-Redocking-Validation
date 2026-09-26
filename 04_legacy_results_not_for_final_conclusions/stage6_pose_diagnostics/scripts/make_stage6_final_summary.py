#!/usr/bin/env python3
"""
Merge Stage 5 and Stage 6 diagnostics into one final, report-ready TSV.

Interpretation:
- All structures in this project show high Murcko-core RMSD.
- Therefore classification is intentionally descriptive, not a universal
  docking-success threshold.
"""

from pathlib import Path
import csv

PROJECT = Path(__file__).resolve().parents[3]

RMSD_TSV = PROJECT / "results/stage5_rmsd_qc/reports/stage5_rmsd_summary.tsv"
VISUAL_TSV = PROJECT / "results/stage5_rmsd_qc/reports/stage5_visual_qc_template.tsv"
COM_TSV = PROJECT / "results/stage6_pose_diagnostics/reports/stage6_pose_diagnostics.tsv"
CORE_TSV = PROJECT / "results/stage6_pose_diagnostics/reports/stage6_murcko_core_rmsd.tsv"

OUT_TSV = PROJECT / "results/stage6_pose_diagnostics/reports/stage6_final_pose_diagnostic_summary.tsv"

PDBIDS = ["1KMV", "3FS6", "3NTZ", "3NXO", "4M6K", "5HSR"]


def read_by_pdbid(path):
    with path.open("r", newline="", encoding="utf-8") as handle:
        return {
            row["pdbid"].strip(): row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def as_float(row, key):
    return float(row[key])


def main():
    rmsd = read_by_pdbid(RMSD_TSV)
    visual = read_by_pdbid(VISUAL_TSV)
    com = read_by_pdbid(COM_TSV)
    core = read_by_pdbid(CORE_TSV)

    rows = []

    for pdbid in PDBIDS:
        full_rmsd = as_float(rmsd[pdbid], "rmsd_heavy_atom_angstrom")
        com_distance = as_float(com[pdbid], "ligand_com_distance_angstrom")
        jaccard = as_float(com[pdbid], "contact_jaccard")
        core_rmsd = as_float(core[pdbid], "direct_scaffold_rmsd_angstrom")
        aligned_core_rmsd = as_float(
            core[pdbid], "aligned_scaffold_rmsd_angstrom"
        )
        core_delta = full_rmsd - core_rmsd

        # This is a descriptive project-specific label, not a universal cutoff.
        if core_rmsd < 2.0 and core_delta >= 1.0:
            diagnostic_class = "tail_dominated_mismatch"
        elif core_rmsd >= 3.0:
            diagnostic_class = "core_dominated_pose_mismatch"
        else:
            diagnostic_class = "mixed_or_intermediate_mismatch"

        if com_distance < 0.5 and jaccard >= 0.4:
            site_interpretation = "same_site_supported_by_COM_and_contacts"
        elif com_distance < 1.0:
            site_interpretation = "site_proximity_supported_but_contact_difference"
        else:
            site_interpretation = "site_assignment_requires_visual_recheck"

        rows.append(
            {
                "pdbid": pdbid,
                "vina_affinity_kcal_mol": rmsd[pdbid].get(
                    "vina_affinity_kcal_mol", ""
                ),
                "full_heavy_atom_rmsd_angstrom": f"{full_rmsd:.3f}",
                "ligand_com_distance_angstrom": f"{com_distance:.3f}",
                "contact_jaccard_4A": f"{jaccard:.3f}",
                "murcko_scaffold_heavy_atom_count": core[pdbid].get(
                    "murcko_scaffold_heavy_atom_count", ""
                ),
                "scaffold_fraction_native": core[pdbid].get(
                    "scaffold_fraction_native", ""
                ),
                "direct_murcko_core_rmsd_angstrom": f"{core_rmsd:.3f}",
                "aligned_murcko_core_rmsd_angstrom": f"{aligned_core_rmsd:.3f}",
                "full_minus_core_rmsd_angstrom": f"{core_delta:.3f}",
                "visual_class": visual[pdbid].get("visual_class", ""),
                "site_interpretation": site_interpretation,
                "diagnostic_class": diagnostic_class,
                "final_interpretation": (
                    "Binding-site region recovered, but native pose is not "
                    "recovered: high receptor-frame core RMSD indicates a "
                    "core-dominated orientation/conformation mismatch."
                ),
            }
        )

    fields = [
        "pdbid",
        "vina_affinity_kcal_mol",
        "full_heavy_atom_rmsd_angstrom",
        "ligand_com_distance_angstrom",
        "contact_jaccard_4A",
        "murcko_scaffold_heavy_atom_count",
        "scaffold_fraction_native",
        "direct_murcko_core_rmsd_angstrom",
        "aligned_murcko_core_rmsd_angstrom",
        "full_minus_core_rmsd_angstrom",
        "visual_class",
        "site_interpretation",
        "diagnostic_class",
        "final_interpretation",
    ]

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote: {OUT_TSV}")
    print(f"Completed: {len(rows)}/{len(PDBIDS)} structures")


if __name__ == "__main__":
    main()
