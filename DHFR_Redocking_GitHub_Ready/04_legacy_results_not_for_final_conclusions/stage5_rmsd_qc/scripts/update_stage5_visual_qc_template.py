from pathlib import Path
import csv

root = Path("results/stage5_rmsd_qc/reports")
rmsd_tsv = root / "stage5_rmsd_summary.tsv"
visual_tsv = root / "stage5_visual_qc_template.tsv"

output_columns = [
    "pdbid",
    "vina_affinity_kcal_mol",
    "rmsd_heavy_atom_angstrom",
    "same_binding_pocket",
    "core_overlap",
    "orientation_match",
    "tail_match",
    "visual_class",
    "overall_judgment",
    "notes",
]

default_manual = {
    "same_binding_pocket": "NA",
    "core_overlap": "NA",
    "orientation_match": "NA",
    "tail_match": "NA",
    "visual_class": "NA",
    "overall_judgment": "NA",
    "notes": "NA",
}

def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

if not rmsd_tsv.exists():
    raise FileNotFoundError(f"Missing RMSD summary: {rmsd_tsv}")

rmsd_rows = read_tsv(rmsd_tsv)
rmsd_map = {row["pdbid"]: row for row in rmsd_rows}

existing_map = {}
if visual_tsv.exists():
    existing_rows = read_tsv(visual_tsv)
    existing_map = {row["pdbid"]: row for row in existing_rows}

final_rows = []
for pdbid in [row["pdbid"] for row in rmsd_rows]:
    rmsd_row = rmsd_map[pdbid]
    old = existing_map.get(pdbid, {})

    row = {
        "pdbid": pdbid,
        "vina_affinity_kcal_mol": rmsd_row.get("vina_affinity_kcal_mol", "NA"),
        "rmsd_heavy_atom_angstrom": rmsd_row.get("rmsd_heavy_atom_angstrom", "NA"),
    }

    for col in default_manual:
        row[col] = old.get(col, default_manual[col]) if old.get(col, "") != "" else default_manual[col]

    final_rows.append(row)

with open(visual_tsv, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=output_columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(final_rows)

print(f"Updated: {visual_tsv}")
for row in final_rows:
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
