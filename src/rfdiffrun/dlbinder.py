import os

def protein_mpnn(dlbind_env,
                 mpnn_path,
                 rf_dir,
                 device_ind,
                 ):
    command = f"""\
              module use /soft/modulefiles &&\
              module load conda &&\
              conda activate {dlbind_env} &&\
              {mpnn_path}/dl_interface_design.py\
              -silent {rf_dir}/rfout.silent\
              -outsilent {rf_dir}/mpnnout.silent\
              -checkpoint_name {rf_dir}/checkpoint_mpnn.dat\
              -relax_cycles 0\
             """
    os.system(command)

