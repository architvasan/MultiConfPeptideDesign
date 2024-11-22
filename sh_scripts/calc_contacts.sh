#!/bin/bash
module load conda
conda activate chai1
cd /eagle/datascience/avasan/Simulations/MultiConfPeptideDesign
cutoff=12
for i in {0..3}
do
    python mcbind/analyze/contact_pairs.py \
        -p inputs/NMNAT-2/conformations/meta${i}.pdb \
        -sA "protein and resid 119-192 and name CA" \
        -sB "protein and name CA" \
        -CA A \
        -CB A \
        -c ${cutoff} \
        -of common/contact_pairs/const_meta${i}_freqs.csv \
        -or common/contact_pairs/contact_meta${i}.restraints
done