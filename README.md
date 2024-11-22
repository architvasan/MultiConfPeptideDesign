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
2. Loading constraints file into chai1
```bash
    python src/fold_ai/chai_pred_context.py \
     -T file -s1 inputs/NMNAT-2/sequences/nmnat2.fasta \
      -cf common/contact_pairs/contact_meta3.restraints \
      -od common/folds/folds_nmnat2/meta3 \
      -f common/folds/folds_nmnat2/meta3/meta3.fasta \
      -d 0
```