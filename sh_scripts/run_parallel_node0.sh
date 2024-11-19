module use /soft/modulefiles
module load conda

cd /eagle/datascience/avasan/Simulations/AntiBodyDesign
conda activate /lus/eagle/projects/datascience/avasan/envs/chai1

module load PrgEnv-gnu
export CRAY_ACCEL_TARGET="nvidia80"
export CRAY_TCMALLOC_MEMFS_FORCE="1"
export CRAYPE_LINK_TYPE="dynamic"
export CRAY_ACCEL_VENDOR="nvidia"
export CRAY_CPU_TARGET="x86-64"

data_dir=iedb_tables_pos_and_neg
output_dir_general=trials/T5_ant3HFM_body4NCO
log_name=logs/run_t5_a3HFM_b4NCO
dbfil=fab_lyso.csv
chaintarget=heavy_chain
cdrid=heavy_cdr3
ndesigns=1

mkdir $output_dir_general

for i in $(seq 0 2);
do
   echo $i
   python run_pipeline_full_newcdrs.py \
    --dbfil $data_dir/$dbfil \
    --cdrfil $data_dir/heavy-cdr3_rank${i}.json \
    --chaintarget $chaintarget \
    --cdrid $cdrid \
    --outdir $output_dir_general/rank$i \
    --ndesigns $ndesigns \
    --device $i \
    > ${log_name}${i}.log 2> ${log_name}${i}.err &
done

echo 3

python run_pipeline_full_newcdrs.py \
 --dbfil $data_dir/$dbfil \
 --cdrfil $data_dir/heavy-cdr3_rank3.json \
 --chaintarget $chaintarget \
 --cdrid $cdrid \
 --outdir $output_dir_general/rank3 \
 --ndesigns $ndesigns \
 --device 3 \
> ${log_name}3.log 2> ${log_name}3.err