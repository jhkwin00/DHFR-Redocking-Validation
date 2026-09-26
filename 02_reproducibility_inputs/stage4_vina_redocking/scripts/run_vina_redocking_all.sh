#!/usr/bin/env bash
set -euo pipefail

BASE="results/stage4_vina_redocking"
PREP="results/prepared_structures"
LIGDIR="results/stage3_rdkit_ligand_prep/sdf_h"

mkdir -p \
  "$BASE/receptor_pdbqt" \
  "$BASE/ligand_pdbqt" \
  "$BASE/config" \
  "$BASE/docking_out" \
  "$BASE/logs" \
  "$BASE/reports" \
  "$BASE/tmp"

# PDBID -> native ligand residue name in reference_complex.pdb
declare -A RESNAME
RESNAME["1KMV"]="MTX"
RESNAME["3FS6"]="DH1"
RESNAME["3NTZ"]="726"
RESNAME["3NXO"]="741"
RESNAME["4M6K"]="18Q"
RESNAME["5HSR"]="0HG"

PDBIDS=(1KMV 3FS6 3NTZ 3NXO 4M6K 5HSR)

SUMMARY="$BASE/reports/vina_redocking_summary.tsv"
echo -e "pdbid\treceptor_pdbqt\tligand_pdbqt\tcenter_x\tcenter_y\tcenter_z\tbest_affinity_kcal_mol\tstatus" > "$SUMMARY"

for pdbid in "${PDBIDS[@]}"; do
  echo "=================================================="
  echo "[INFO] Processing $pdbid"
  echo "=================================================="

  receptor_pdb="$PREP/${pdbid}_receptor.pdb"
  complex_pdb="$PREP/${pdbid}_reference_complex.pdb"
  ligand_sdf="$LIGDIR/${pdbid}_ligand_only_H.sdf"

  receptor_prefix="$BASE/receptor_pdbqt/${pdbid}_receptor"
  ligand_pdbqt="$BASE/ligand_pdbqt/${pdbid}_ligand.pdbqt"
  config_file="$BASE/config/${pdbid}_vina_config.txt"
  out_pdbqt="$BASE/docking_out/${pdbid}_docked.pdbqt"
  log_file="$BASE/logs/${pdbid}_vina.log"
  receptor_log="$BASE/logs/${pdbid}_receptor_prep.log"
  ligand_log="$BASE/logs/${pdbid}_ligand_prep.log"

  resname="${RESNAME[$pdbid]}"

  if [[ ! -f "$receptor_pdb" ]]; then
    echo "[ERROR] Missing receptor PDB: $receptor_pdb"
    echo -e "${pdbid}\tNA\tNA\tNA\tNA\tNA\tNA\tmissing_receptor_pdb" >> "$SUMMARY"
    continue
  fi

  if [[ ! -f "$complex_pdb" ]]; then
    echo "[ERROR] Missing reference complex PDB: $complex_pdb"
    echo -e "${pdbid}\tNA\tNA\tNA\tNA\tNA\tNA\tmissing_reference_complex" >> "$SUMMARY"
    continue
  fi

  if [[ ! -f "$ligand_sdf" ]]; then
    echo "[ERROR] Missing ligand SDF: $ligand_sdf"
    echo -e "${pdbid}\tNA\tNA\tNA\tNA\tNA\tNA\tmissing_ligand_sdf" >> "$SUMMARY"
    continue
  fi

  echo "[INFO] Preparing receptor..."
  if ! mk_prepare_receptor.py \
      --read_pdb "$receptor_pdb" \
      -o "$receptor_prefix" \
      -p > "$receptor_log" 2>&1; then
    echo "[ERROR] Receptor preparation failed for $pdbid"
    echo -e "${pdbid}\tNA\tNA\tNA\tNA\tNA\tNA\treceptor_prep_failed" >> "$SUMMARY"
    continue
  fi

  receptor_pdbqt=$(find "$BASE/receptor_pdbqt" -maxdepth 1 -type f \
    $$ -name "${pdbid}_receptor.pdbqt" -o -name "${pdbid}_receptor_rigid.pdbqt" $$ | head -n 1 || true)

  if [[ -z "${receptor_pdbqt:-}" || ! -f "$receptor_pdbqt" ]]; then
    echo "[ERROR] Receptor PDBQT not found for $pdbid"
    echo -e "${pdbid}\tNA\tNA\tNA\tNA\tNA\tNA\treceptor_pdbqt_missing" >> "$SUMMARY"
    continue
  fi

  echo "[INFO] Preparing ligand..."
  if ! mk_prepare_ligand.py \
      -i "$ligand_sdf" \
      -o "$ligand_pdbqt" > "$ligand_log" 2>&1; then
    echo "[ERROR] Ligand preparation failed for $pdbid"
    echo -e "${pdbid}\t${receptor_pdbqt}\tNA\tNA\tNA\tNA\tNA\tligand_prep_failed" >> "$SUMMARY"
    continue
  fi

  if [[ ! -f "$ligand_pdbqt" ]]; then
    echo "[ERROR] Ligand PDBQT not found for $pdbid"
    echo -e "${pdbid}\t${receptor_pdbqt}\tNA\tNA\tNA\tNA\tNA\tligand_pdbqt_missing" >> "$SUMMARY"
    continue
  fi

  echo "[INFO] Calculating box center from native ligand $resname ..."
  center_line=$(python - <<PY
pdb_file = "$complex_pdb"
resname = "$resname"

coords = []
with open(pdb_file) as f:
    for line in f:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        rn = line[17:20].strip()
        if rn != resname:
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coords.append((x, y, z))

if not coords:
    raise SystemExit("FAILED: native ligand not found")

cx = sum(x for x, y, z in coords) / len(coords)
cy = sum(y for x, y, z in coords) / len(coords)
cz = sum(z for x, y, z in coords) / len(coords)
print(f"{cx:.3f}\t{cy:.3f}\t{cz:.3f}")
PY
)

  center_x=$(echo "$center_line" | cut -f1)
  center_y=$(echo "$center_line" | cut -f2)
  center_z=$(echo "$center_line" | cut -f3)

  cat > "$config_file" <<CFG
receptor = $receptor_pdbqt
ligand = $ligand_pdbqt

center_x = $center_x
center_y = $center_y
center_z = $center_z

size_x = 22
size_y = 22
size_z = 22

exhaustiveness = 16
num_modes = 20
energy_range = 4
CFG

  echo "[INFO] Running Vina..."
  if ! vina \
      --config "$config_file" \
      --out "$out_pdbqt" \
      --log "$log_file"; then
    echo "[ERROR] Vina failed for $pdbid"
    echo -e "${pdbid}\t${receptor_pdbqt}\t${ligand_pdbqt}\t${center_x}\t${center_y}\t${center_z}\tNA\tvina_failed" >> "$SUMMARY"
    continue
  fi

  best_affinity=$(awk '/^[[:space:]]*1[[:space:]]+-/ {print $2; exit}' "$log_file")
  if [[ -z "${best_affinity:-}" ]]; then
    best_affinity="NA"
    status="vina_done_but_score_not_parsed"
  else
    status="ok"
  fi

  echo -e "${pdbid}\t${receptor_pdbqt}\t${ligand_pdbqt}\t${center_x}\t${center_y}\t${center_z}\t${best_affinity}\t${status}" >> "$SUMMARY"
  echo "[INFO] Done: $pdbid  best_affinity=${best_affinity}"
done

echo
echo "[INFO] Summary written to: $SUMMARY"
