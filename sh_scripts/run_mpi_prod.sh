#!/bin/bash
#PBS -N rfdiff
#PBS -l select=10
#PBS -l walltime=6:00:00
#PBS -q preemptable
#PBS -A datascience
#PBS -l filesystems=eagle
#PBS -m abe
#PBS -M avasan@anl.gov

module use /soft/modulefiles
module load conda

conda activate /lus/eagle/projects/datascience/avasan/envs/chai1
conda env list
#module load mpiwrappers/cray-mpich-llvm 
module load PrgEnv-gnu 
export CRAY_ACCEL_TARGET="nvidia80" 
export CRAY_TCMALLOC_MEMFS_FORCE="1" 
export CRAYPE_LINK_TYPE="dynamic" 
export CRAY_ACCEL_VENDOR="nvidia"
export CRAY_CPU_TARGET="x86-64"

trial=12

cd /eagle/datascience/avasan/Simulations/AntiBodyDesign
stdout=logs/run_t${trial}.mpi.a3HFM_b4NCO.log
stderr=logs/run_t${trial}.mpi.a3HFM_b4NCO.err
NDEPTH=16

data_dir=iedb_tables_pos_and_neg
output_dir_general=trials/T${trial}_ant3HFM_body4NCO
log_name=logs/run_t${trial}_a3HFM_b4NCO
dbfil=fab_lyso.csv
chaintarget=heavy_chain
cdrid=heavy_cdr3
ndesigns=40

totalrank=40
pernode=4

mkdir $output_dir_general

mpiexec -n $totalrank -ppn $pernode \
    --depth=${NDEPTH} --cpu-bind depth \
    ./set_affinity_gpu_polaris.sh \
    /lus/eagle/projects/datascience/avasan/envs/chai1/bin/python run_pipeline_full_newcdrs.py \
    --dbfil $data_dir/$dbfil \
    --cdrfil $data_dir/heavy-cdr3_all.json \
    --chaintarget $chaintarget \
    --usempi \
    --cdrid $cdrid \
    --outdir $output_dir_general \
    --ndesigns $ndesigns \
    > $stdout 2> $stderr