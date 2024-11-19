import numpy as np
import json
import pandas as pd

def split_savelocs(
                    central_loc_file,
                    outpattern,
                    number_ranks=8,
                    ):

    df = pd.read_csv(central_loc_file)#'./all_pdbs_save_loc.t16.csv')
    df_split = np.array_split(df, number_ranks)
    for i in range(0,len(df_split)):
        df_split[i].to_csv(f'{outpattern}_{i}.csv', index=False)
    return df_split

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-fc',
                        '--file_central',
                        type=str,
                        help='input csv file with paths of pdbs')

    parser.add_argument('-op',
                        '--outpattern',
                        type=str,
                        help='output_pattern')
    
    parser.add_argument('-R',
                        '--ranks',
                        type=int,
                        required=False,
                        default=8,
                        help='number of ranks')

    args = parser.parse_args()

    df_split = split_savelocs(
                    args.file_central,
                    args.outpattern,
                    number_ranks=args.ranks,
                    )

# if False:
#     with open('./all_pdbs_save_loc.t16.csv','r') as file:
#         cdrlist = json.load(file)
    
#     cdrlist_set = set(cdrlist['all_train_sequences'])
#     cdrlist_uniq = list(cdrlist_set)
#     with open(f'iedb_tables_pos_and_neg/heavy-cdr3_all.json', 'w') as file:
#         json.dump({'all_train_sequences': list(cdrlist_uniq)}, file)
    
#     cdrlist_split = np.array_split(cdrlist_uniq, 8)
#     for i in range(8):
#         with open(f'iedb_tables_pos_and_neg/heavy-cdr3_rank{i}.json', 'w') as file:
#             json.dump({'all_train_sequences': list(cdrlist_split[i])}, file)
