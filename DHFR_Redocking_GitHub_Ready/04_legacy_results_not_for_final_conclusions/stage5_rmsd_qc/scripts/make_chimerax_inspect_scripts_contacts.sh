#!/usr/bin/env bash
set -euo pipefail

root="$HOME/lab_prep/DHFR_project"
outdir="$root/results/stage5_rmsd_qc/scripts/chimerax_contacts"
mkdir -p "$outdir"

pdbids=(1KMV 3FS6 3NTZ 3NXO 4M6K 5HSR)

for pdbid in "${pdbids[@]}"; do
    receptor="$root/results/prepared_structures/${pdbid}_receptor.pdb"
    native="$root/results/stage2_ligand_extraction/ligand_only/${pdbid}_ligand_only.pdb"
    docked="$root/results/stage5_rmsd_qc/docked_top1/${pdbid}_top1_heavy_from_pdbqt.sdf"

    receptor_win=$(wslpath -w "$receptor")
    native_win=$(wslpath -w "$native")
    docked_win=$(wslpath -w "$docked")

    cat > "$outdir/${pdbid}_inspect_contacts.cxc" <<EOF
close all
open "$receptor_win"
open "$native_win"
open "$docked_win"

hide all models
show #1 cartoons
hide #1 atoms
show #2 models
show #3 models

style #2 stick
style #3 stick

color #1 gray
color #2 cyan
color #3 yellow

show protein & within 5 of #2 models
show protein & within 5 of #3 models
style protein & within 5 of #2,#3 stick
color protein & within 5 of #2,#3 wheat

label off
set bgColor white
view
EOF
done

echo "Generated ChimeraX contact scripts in: $outdir"
ls -1 "$outdir"
