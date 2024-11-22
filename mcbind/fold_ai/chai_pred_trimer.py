from pathlib import Path

import numpy as np
import torch

import os
from chai_lab.chai1 import run_inference

def fold_chai(s1,
              s2,
              s3,
              fastapath,
              outputdir,
              device=0,
              ):

    fasta = f"""
    >protein|name=prot1
{s1}
>protein|name=prot2
{s2}
>protein|name=prot3
{s3}
    """.strip()
    
    fasta_path = Path(fastapath)
    fasta_path.write_text(fasta)
    output_dir = Path(outputdir)

    candidates = run_inference(
        fasta_file=fasta_path,
        output_dir=output_dir,
        # 'default' setup
        num_trunk_recycles=3,
        num_diffn_timesteps=80,
        seed=42,
        device=torch.device(f"cuda"),
        use_esm_embeddings=True,
    )

    cif_paths = candidates.cif_paths
    scores = [rd.aggregate_score for rd in candidates.ranking_data]
    scores = np.load(output_dir.joinpath("scores.model_idx_2.npz"))
    return scores


if __name__ == "__main__": 
    import argparse
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-T',
                        '--typeseq',
                        type=str,
                        help='how are sequences formatted? (file/str)')

    parser.add_argument('-s1',
                        '--seq1',
                        type=str,
                        help='Sequence 1')
    
    parser.add_argument('-s2',
                        '--seq2',
                        type=str,
                        required=False,
                        default="none",
                        help='Sequence 2')

    parser.add_argument('-s3',
                        '--seq3',
                        type=str,
                        required=False,
                        default="none",
                        help='Sequence 3')

    parser.add_argument('-od',
                        '--outdir',
                        type=str,
                        required=True,
                        help='output directory')

    parser.add_argument('-f',
                        '--fastapath',
                        type=str,
                        required=True,
                        help='fastapath')


    parser.add_argument('-d',
                        '--device',
                        type=str,
                        required=False,
                        default=0,
                        help='Device to place job')
    
    args = parser.parse_args()
    
    #scores = fold_chai(
    #          s1,
    #          s2,
    #          args.fastapath,
    #          args.outputdir,
    #          )

    os.environ['CHAI_DOWNLOADS_DIR'] = '/lus/eagle/projects/datascience/avasan/Software/CHAI1_downloads'
    try:
        os.mkdir(f'{args.outdir}')
    except:
        print("couldn't make dir")
        pass

    if args.typeseq == 'str':
        scores = fold_chai(
                args.seq1,
                args.seq2,
                args.seq3,
                args.fastapath,
                args.outdir,
                )

    elif args.typeseq == 'file':
        with open(args.seq1, 'r') as file:
            seq1_str = file.read().replace('\n', '')

        if args.seq2 != "none":
            with open(args.seq2, 'r') as file:
                seq2_str = file.read().replace('\n', '')
        else:
            seq2_str = args.seq2

        if args.seq3 != "none":
            with open(args.seq3, 'r') as file:
                seq3_str = file.read().replace('\n', '')
        else:
            seq3_str = args.seq3

        scores = fold_chai(
                seq1_str,
                seq2_str,
                seq3_str,
                args.fastapath,
                args.outdir,
                )


    print(scores)

