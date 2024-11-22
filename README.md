# MultiConfPeptideDesign

#### getting contacts to input into chai1:
1. obtain constraints.restraint file
    over metastates 0-3 
``` bash
python src/analyze/contact_pairs.py \
        -p inputs/NMNAT-2/conformations/meta1.pdb \
        -sA "protein and resid 119-192 and name CA" \
        -sB "protein and name CA" \
        -CA A \
        -CB B \
        -c 8 \
        -of common/contact_pairs/const_meta1_freqs.csv \
        -or common/contact_pairs/contact_meta1.restraints
 ```