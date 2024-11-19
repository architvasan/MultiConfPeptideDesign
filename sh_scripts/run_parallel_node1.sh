module use /soft/modulefiles
module load conda

cd /eagle/datascience/avasan/Simulations/AntiBodyDesign
conda activate /lus/eagle/projects/datascience/avasan/envs/chai1

data_dir=iedb_tables_pos_and_neg
output_dir_general=trials/T5_ant3HFM_body4NCO
dbfil=fab_lyso.csv
chaintarget=heavy_chain
cdrid=heavy_cdr3
ndesigns=1
log_name=logs/run_t5_a3HFM_b4NCO
mkdir $output_dir_general
for i in $(seq 0 2);
do
   let j=i+4
#    python run_pipeline_full_newcdrs.py \
#     --dbfil $data_dir/$dbfil \
#     --cdrfil $data_dir/heavy-cdr3_rank${j}.json \
#     --chaintarget $chaintarget \
#     --cdrid $cdrid \
#     --ndesigns $ndesigns \
#     --device $i &
   echo $j
   python run_pipeline_full_newcdrs.py \
    --dbfil $data_dir/$dbfil \
    --cdrfil $data_dir/heavy-cdr3_rank${j}.json \
    --chaintarget $chaintarget \
    --cdrid $cdrid \
    --outdir $output_dir_general/rank$j \
    --ndesigns $ndesigns \
    --device $i \
    > ${log_name}${j}.log 2> ${log_name}${j}.err &

done
j=7
i=3
echo $j
python run_pipeline_full_newcdrs.py \
 --dbfil $data_dir/$dbfil \
 --cdrfil $data_dir/heavy-cdr3_rank${j}.json \
 --chaintarget $chaintarget \
 --cdrid $cdrid \
 --outdir $output_dir_general/rank$j \
 --ndesigns $ndesigns \
 --device $i \
 > ${log_name}${j}.log 2> ${log_name}${j}.err