import pandas as pd
import numpy as np
import glob
import json
from tqdm import tqdm
# import module sars_functions.py
import sars_functions
import functions.AmyloComp_functions as amp

### Open Score matrices
matr_path = 'matrices/'
MatrixOuter =      amp.OpenScores(matr_path + 'ScoreOuter_2.txt')
MatrixArc =        amp.OpenScores(matr_path + 'ScoreArc_5.txt')
MatrixInner =      amp.OpenScores(matr_path + 'result_matrix_b4.tsv')
MatrixInnerExtra = amp.OpenScores(matr_path + 'result_matrix_extra_b4.tsv')

arc_score_path = 'matrices/ScoreArc_5.txt'
compatible_arches = sars_functions.eval_compatible_arches(arc_score_path)

human_path = 'path/to/dir/with/human_beta_arches.csv'
iupred_res_path = 'path/to/iupred/res'
sars_path = 'path/to/SARS_beta_arches.csv'
res_path = 'path/where/we/want/to/save/results/'


# write results to out
sars_functions.concat_iupred_res(iupred_res_path, 'out.txt')
human_df = sars_functions.concat_human_csv(human_path)

full_proteome_dict = sars_functions.read_web_iupred('out.txt')
human_df = sars_functions.make_iupred_df(human_df, full_proteome_dict)
sars_df = sars_functions.read_sars_csv(sars_path)

###

## if you would like to save dict and then be able to read it use json
# with open('proteome_unstructured_dict.txt', 'w') as outfile:
#       outfile.write(json.dumps(full_proteome_dict))

# with open('proteome_unstructured_dict.txt') as json_file:
#     full_proteome_dict = json.load(json_file)
# for key in full_proteome_dict.keys():
#     full_proteome_dict[key] = {int(k):v for k,v in full_proteome_dict[key].items()}

###

iupred_df = human_df[human_df.IUpred == 1]
iupred_df_filt = iupred_df[iupred_df['Score_total'] >= 0.55]
sars_df_filt = sars_df[sars_df["Score_total"] > 0.55]

sars_functions.find_amyloid_connection(sars_df_filt, iupred_df_filt, compatible_arches, res_path)
sars_functions.detect_pathogenic_amyloids(res_path,amyloidase_base,additional_base)