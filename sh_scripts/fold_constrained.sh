#!/bin/bash
module load conda
conda activate chai1
cd /eagle/datascience/avasan/Simulations/MultiConfPeptideDesign
for i in {0..3}
do
    python mcbind/fold_ai/chai_pred_context.py \
     -T file -s1 inputs/NMNAT-2/sequences/nmnat2.fasta \
      -cf common/contact_pairs/contact_meta${i}.restraints \
      -od common/folds/folds_nmnat2/meta${i} \
      -f common/folds/folds_nmnat2/meta${i}/fold.fasta \
      -d 0
done