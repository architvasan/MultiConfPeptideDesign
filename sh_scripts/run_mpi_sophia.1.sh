#!/bin/bash

module use /soft/modulefiles
module load conda

cd /eagle/datascience/avasan/Simulations/AntiBodyDesign
conda activate /lus/eagle/projects/datascience/avasan/envs/chai1
conda env list
#module load mpiwrappers/cray-mpich-llvm 

inputfil=iedb_tables_pos_and_neg/hiv_ab.csv
chaiout=trials/T1/chaiout
rfout=trials/T2/rfout
logout=trials/T2/logout
resmaploc=trials/T1/resmaps
numdesign=10
foldinit=False
totalrank=1
pernode=1
stdout=logs/rf_t2.mpi.sophia.log
stderr=logs/rf_t2.mpi.sophia.err


CUDA_VISIBLE_DEVICES=0 python run_pipeline_full_mpi4py.py \
    -i $inputfil \
    -C $chaiout \
    -R $rfout/0_0 \
    -L $logout \
    -M $resmaploc \
    -N $numdesign \
    -G $pernode \
    -F $foldinit \
        > $stdout.0.log 2> $stderr.0.log &

CUDA_VISIBLE_DEVICES=1 python run_pipeline_full_mpi4py.py \
    -i $inputfil \
    -C $chaiout \
    -R $rfout/1_0 \
    -L $logout \
    -M $resmaploc \
    -N $numdesign \
    -G $pernode \
    -F $foldinit \
        > $stdout.1.log 2> $stderr.1.log &

CUDA_VISIBLE_DEVICES=2 python run_pipeline_full_mpi4py.py \
    -i $inputfil \
    -C $chaiout \
    -R $rfout/2_0 \
    -L $logout \
    -M $resmaploc \
    -N $numdesign \
    -G $pernode \
    -F $foldinit \
        > $stdout.2.log 2> $stderr.2.log &

CUDA_VISIBLE_DEVICES=3 python run_pipeline_full_mpi4py.py \
    -i $inputfil \
    -C $chaiout \
    -R $rfout/3_0 \
    -L $logout \
    -M $resmaploc \
    -N $numdesign \
    -G $pernode \
    -F $foldinit \
        > $stdout.3.log 2> $stderr.3.log 

#mpiexec -n $totalrank -ppn $pernode \
#    python run_pipeline_full_mpi4py.py \
#        -i $inputfil \
#        -C $chaiout \
#        -R $rfout \
#        -L $logout \
#        -M $resmaploc \
#        -N $numdesign \
#        -G $pernode \
#        -F $foldinit \
#            > $stdout 2> $stderr