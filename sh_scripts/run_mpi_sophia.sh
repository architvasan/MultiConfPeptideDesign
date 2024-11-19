#!/bin/bash
#PBS -N rfdiff_14
#PBS -l select=10
#PBS -l walltime=1:00:00
#PBS -q debug-scaling
#PBS -A datascience
#PBS -l filesystems=eagle
#PBS -m abe
#PBS -M avasan@anl.gov

module use /soft/modulefiles
module load conda

conda activate /lus/eagle/projects/datascience/avasan/envs/chai1
conda env list
#module load mpiwrappers/cray-mpich-llvm 
#module load PrgEnv-gnu 
#export CRAY_ACCEL_TARGET="nvidia80" 
#export CRAY_TCMALLOC_MEMFS_FORCE="1" 
#export CRAYPE_LINK_TYPE="dynamic" 
#export CRAY_ACCEL_VENDOR="nvidia"
#export CRAY_CPU_TARGET="x86-64"

trial=1
cd /eagle/datascience/avasan/Simulations/AntiBodyDesign_Constrained
PWD=/eagle/datascience/avasan/Simulations/AntiBodyDesign_Constrained
stdout=logs/run_t${trial}.mpi.a3HFM_b4NCO.log
stderr=logs/run_t${trial}.mpi.a3HFM_b4NCO.err
NDEPTH=16

data_dir=iedb_tables_pos_and_neg
output_dir_general=trials/T${trial}_ant3HFM_body4NCO_const
log_name=logs/run_t${trial}_a3HFM_b4NCO
dbfil=fab_lyso.csv
chaintarget=heavy_chain
cdrid=heavy_cdr3
ndesigns=20

totalrank=8
pernode=8

mkdir $output_dir_general

# Run each rank separately on sophia since there are issues with
# mpi4py

rank=0
for i in $(seq 0 6);
do
    rank=${i}
    echo $i
    /lus/eagle/projects/datascience/avasan/envs/chai1/bin/python src/run_pipeline_full_newcdrs.py \
        --dbfil $PWD/$data_dir/$dbfil \
        --cdrfil $PWD/$data_dir/heavy-cdr3_rank${rank}.json \
        --chaintarget $chaintarget \
        --cdrid $cdrid \
        --outdir $PWD/$output_dir_general/${rank} \
        --ndesigns $ndesigns \
        -G $pernode \
        -D ${rank} > $PWD/$log_name.${rank}.log 2> $PWD/$log_name.$rank.err &
done    

rank=7
echo $i
/lus/eagle/projects/datascience/avasan/envs/chai1/bin/python src/run_pipeline_full_newcdrs.py \
    --dbfil $PWD/$data_dir/$dbfil \
    --cdrfil $PWD/$data_dir/heavy-cdr3_rank${rank}.json \
    --chaintarget $chaintarget \
    --cdrid $cdrid \
    --outdir $PWD/$output_dir_general/${rank} \
    --ndesigns $ndesigns \
    -G $pernode \
    -D ${rank} > $PWD/$log_name.${rank}.log 2> $PWD/$log_name.$rank.err



#mpiexec -n $totalrank -ppn $pernode \
#    --depth=${NDEPTH} --cpu-bind depth \
#    ./set_affinity_gpu_polaris.sh \
#    /lus/eagle/projects/datascience/avasan/envs/chai1/bin/python run_pipeline_full_newcdrs.py \
#    --dbfil $data_dir/$dbfil \
#    --cdrfil $data_dir/heavy-cdr3_all.json \
#    --chaintarget $chaintarget \
#    --usempi \
#    --cdrid $cdrid \
#    --outdir $output_dir_general \
#    --ndesigns $ndesigns \
#    -G $pernode \
#    > $stdout 2> $stderr
