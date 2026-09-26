from pathlib import Path
import csv

visual_tsv = Path("results/stage5_rmsd_qc/reports/stage5_visual_qc_template.tsv")

updates = {
    "3FS6": {
        "same_binding_pocket": "yes",
        "core_overlap": "moderate",
        "orientation_match": "moderate",
        "tail_match": "poor",
        "visual_class": "core_overlap_tail_shifted",
        "overall_judgment": "partial_pose_recovery",
        "notes": "core overlap present, substituent/tail mismatch",
    },
    "4M6K": {
        "same_binding_pocket": "yes",
        "core_overlap": "moderate",
        "orientation_match": "poor",
        "tail_match": "moderate",
        "visual_class": "same_pocket_rotated",
        "overall_judgment": "site_recovered_pose_poor",
        "notes": "same pocket, noticeable orientation mismatch",
    },
    "5HSR": {
        "same_binding_pocket": "yes",
        "core_overlap": "moderate",
        "orientation_match": "moderate",
        "tail_match": "moderate",
        "visual_class": "core_overlap_tail_shifted",
        "overall_judgment": "partial_pose_recovery",
        "notes": "same pocket, partial core overlap with modest pose mismatch",
    },
}

with open(visual_tsv, newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
    fieldnames = rows[0].keys()

for row in rows:
    pdbid = row["pdbid"]
    if pdbid in updates:
        for k, v in updates[pdbid].items():
            row[k] = v

with open(visual_tsv, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Updated: {visual_tsv}")
for row in rows:
    print(
        row["pdbid"],
        row["vina_affinity_kcal_mol"],
        row["rmsd_heavy_atom_angstrom"],
        row["same_binding_pocket"],
        row["core_overlap"],
        row["orientation_match"],
        row["tail_match"],
        row["visual_class"],
        row["overall_judgment"],
        row["notes"],
        sep="\t"
    )
