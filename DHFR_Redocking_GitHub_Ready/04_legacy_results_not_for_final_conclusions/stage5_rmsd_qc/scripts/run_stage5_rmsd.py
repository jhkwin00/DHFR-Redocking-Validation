from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolAlign

ROOT = Path("results")
STAGE3_SDF = ROOT / "stage3_rdkit_ligand_prep" / "sdf"
DOCK_DIR = ROOT / "stage4_vina_redocking" / "docking_out"
FINAL_SUMMARY = ROOT / "stage4_vina_redocking" / "reports" / "vina_redocking_summary_v2_final.tsv"

OUTDIR = ROOT / "stage5_rmsd_qc"
DOCKED_TOP1_DIR = OUTDIR / "docked_top1"
REPORT_DIR = OUTDIR / "reports"

DOCKED_TOP1_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

PDBIDS = ["1KMV", "3FS6", "3NTZ", "3NXO", "4M6K", "5HSR"]

def read_affinity_map(tsv_path):
    d = {}
    lines = tsv_path.read_text().splitlines()
    header = lines[0].split("\t")
    idx_pdbid = header.index("pdbid")
    idx_aff = header.index("best_affinity_kcal_mol")
    idx_status = header.index("status")
    for line in lines[1:]:
        parts = line.split("\t")
        d[parts[idx_pdbid]] = {
            "affinity": parts[idx_aff],
            "stage4_status": parts[idx_status],
        }
    return d

def load_native_noh_template(pdbid):
    sdf = STAGE3_SDF / f"{pdbid}_ligand_only.sdf"
    if not sdf.exists():
        raise FileNotFoundError(f"Missing template SDF: {sdf}")
    mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Failed to read template SDF: {sdf}")
    mol = Chem.RemoveHs(mol)
    return mol, sdf

def is_hydrogen_pdbqt_line(line):
    parts = line.split()
    if not parts:
        return False
    atom_type = parts[-1].upper()
    atom_name = parts[2].upper() if len(parts) > 2 else ""
    # PDBQT atom type starts with H for hydrogens: H, HD, HS, etc.
    if atom_type.startswith("H"):
        return True
    # backup rule
    if atom_name.startswith("H"):
        return True
    return False

def extract_model1_heavy_coords_from_pdbqt(pdbqt_path):
    coords = []
    in_model1 = False
    saw_model = False

    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                model_no = line.split()[1]
                if model_no == "1":
                    in_model1 = True
                    saw_model = True
                    continue
                elif saw_model:
                    break

            if line.startswith("ENDMDL") and in_model1:
                break

            if line.startswith(("ATOM", "HETATM")):
                if saw_model and not in_model1:
                    continue
                if is_hydrogen_pdbqt_line(line):
                    continue
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))

    return coords

def build_mol_with_coords(template_mol, coords):
    mol = Chem.Mol(template_mol)
    if mol.GetNumAtoms() != len(coords):
        raise ValueError(
            f"Heavy-atom count mismatch: template={mol.GetNumAtoms()} pdbqt_heavy={len(coords)}"
        )
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, (x, y, z))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol

def write_sdf(mol, path):
    w = Chem.SDWriter(str(path))
    w.write(mol)
    w.close()

affinity_map = read_affinity_map(FINAL_SUMMARY)
rows = []

for pdbid in PDBIDS:
    pdbqt = DOCK_DIR / f"{pdbid}_docked.pdbqt"
    aff = affinity_map.get(pdbid, {}).get("affinity", "NA")
    st4 = affinity_map.get(pdbid, {}).get("stage4_status", "NA")

    try:
        native_mol, native_path = load_native_noh_template(pdbid)
        coords = extract_model1_heavy_coords_from_pdbqt(pdbqt)

        if not coords:
            raise ValueError("No heavy-atom coordinates found in MODEL 1")

        docked_mol = build_mol_with_coords(native_mol, coords)
        rmsd = rdMolAlign.GetBestRMS(docked_mol, native_mol)

        docked_sdf = DOCKED_TOP1_DIR / f"{pdbid}_top1_heavy_from_pdbqt.sdf"
        native_sdf_copy = DOCKED_TOP1_DIR / f"{pdbid}_native_noh_template_used.sdf"

        write_sdf(docked_mol, docked_sdf)
        write_sdf(native_mol, native_sdf_copy)

        rows.append([
            pdbid,
            aff,
            f"{rmsd:.3f}",
            str(native_mol.GetNumAtoms()),
            st4
        ])

    except Exception as e:
        rows.append([
            pdbid,
            aff,
            "NA",
            "NA",
            f"failed:{e}"
        ])

out_tsv = REPORT_DIR / "stage5_rmsd_summary.tsv"
with open(out_tsv, "w") as f:
    f.write("pdbid\tvina_affinity_kcal_mol\trmsd_heavy_atom_angstrom\theavy_atom_count\tstage4_status_or_note\n")
    for r in rows:
        f.write("\t".join(r) + "\n")

print(f"Wrote: {out_tsv}")
for r in rows:
    print("\t".join(r))
