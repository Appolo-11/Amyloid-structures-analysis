import pandas as pd
import numpy as np
import glob
import json
from tqdm import tqdm
from joblib import Parallel, delayed
import pickle
import functions.AmyloComp_functions as amp
import functions.BArchClass as ba


def concat_human_csv(path):
    """ These functions concat all .csv files contained in path directory
    """
    all_files = glob.glob(path + "/*.csv")

    # sort file names from _0 to _20000
    all_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))

    df_list = []

    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0, quotechar="'")
        df['ID'] = list(map(lambda x: x.split('|')[1], df['ID']))
        if df.isnull().values.any():
          print(f'{filename[-15:]} has missing values')
        # df.Start = int(df.Start)
        df_list.append(df)

    frame = pd.concat(df_list, axis=0, ignore_index=True)
    clean_df = frame.dropna()
    clean_df.Start=clean_df.Start.astype(np.int64)
    clean_df.Stop=clean_df.Stop.astype(np.int64)
    clean_df['IUpred'] = pd.Series(dtype=np.int64)
    clean_df.replace('"', '')
    clean_df.columns.str.replace(r"[' ']", "_", regex=True)

    new_cols = ['ID', 'Sequence', 'Arch', 'Start', 'Score_Bstrand_length',
            'Score_total', 'IUpred']
    for col in clean_df.columns:
      if col not in new_cols:
        clean_df.drop(col, inplace=True, axis=1)

    return clean_df

def read_sars_csv(sars_path):
  """ These functions read .csv files defined with path
  """
  sars_df = pd.read_csv(sars_path, index_col=None, header=0, quotechar="'")
  sars_df['ID'] = list(map(lambda x: x.split('|')[1], sars_df['ID']))
  if sars_df.isnull().values.any():
      print(f'{filename[50:]} has missing values')

  sars_df.columns = sars_df.columns.str.replace(r"[' ']", "_", regex=True)

  new_cols = ['ID', 'Sequence', 'Arch', 'Start', 'Score_Bstrand_length', 'Score_total',
              'CompArch_sars_seq', 'CompArch_sars_start', 'HumanID',
              'CompArch_human_seq', 'CompArch_human_start', 'CompScore', 'OtherScores']

  # delete extra columns
  for col in sars_df.columns:
    if col not in new_cols:
      sars_df.drop(col, inplace=True, axis=1)

  return sars_df


# df should have ID column with proteins id's
def write_unique_proteins(df):
	with open('unique_proteins.txt', 'w') as f:
    		f.writelines(df.ID.unique()+'\n')


def read_web_iupred(filepath):
    """ This function read IUpred web-interface result *filepath* which contain info about proteins position
    	and their IUpred score and return information about unstructured regions.
    		Format of return is a dictionary with keys = protein-IDs and
    		values = dictionaries as like {start_1: stop_1, ... , start_n: stop_n}
    		start_i is an index of start i-th unstructured region, stop_i is the end of that region
    	Protein region considered as unctructured in case of:
    			- IUpred score of positions it contains are greater than 0.3
    			- IUpred score of positions it contains are less than 0.3
    					but total region length are less than 30 aminoacids (AA)
    			- protein length is less than 45 AA: in this case whole protein considered as unstructured

    Important! At the end of the file please add extra string "# IU"
    """
    with open(filepath, 'r') as iupred:
        proteome_dict = {}
        read_unstructured_region = False
        current_pos = 1000
        # count_structured = 0
        structured_start = 1
        current_unstructured_start=1000
        protein=''

        for line in iupred.readlines():

          # all proteins shorter than 45 AA we consider as unstructured
          if line.startswith('# IU') and current_pos <= 45:
                proteome_dict[protein] = {1: current_pos}

          # new protein = initializing region as a structured by default
          if line.startswith('# IU') and read_unstructured_region:
              read_unstructured_region = False

          # if we end to read protein at structured region it's necessary to check
          # was that length more than 30. If not - that region became unstructured
          elif line.startswith('# IU') and not read_unstructured_region:
              if current_pos - structured_start < 30:
                  proteome_dict[protein][structured_start] = current_pos

          # read new protein ID
          elif line[0]=='>':
              # count_structured = 0
              structured_start = 1
              protein = line[1:-1]
              if protein not in proteome_dict:
                proteome_dict[protein] = {}

          # check if it is information about protein position
          # (string starts with protein position number)
          elif line[:1].isdigit():
              split_line = line.split()
              current_pos = int(split_line[0])
              iupred_score = float(split_line[2])
              # anchor_score = float(split_line[3])

              if iupred_score >= 0.3:
                  if not read_unstructured_region:
                      structured_len = current_pos - structured_start + 1
                      current_unstructured_start = current_pos
                      # check if previous structured length is enough to be considered as structured
                      if structured_len < 30:
                      	  # if it's short we note it as unstructured
                          proteome_dict[protein][structured_start] = current_pos
                      else:
                          proteome_dict[protein][current_unstructured_start] = current_pos
                      # turn on unstructured flag
                      read_unstructured_region = True
                  # if we already read unstructured region just change its end to the current_pos
                  else:
                      proteome_dict[protein][current_unstructured_start] = current_pos

              # if it is chenge from unstructured to sctructured we note new structured start and
              # turn off unstructured flag
              elif iupred_score < 0.3 and read_unstructured_region:
                  read_unstructured_region = False
                  structured_start = current_pos
              # else:
              #     count_structured += 1

    # there we unite region-neighbours for one protein
    # region-neighbours are such one: 
    # 		- one's end is another's start
    # 		- one's end is another's start-1
    for protein in proteome_dict.keys():
        protein_dict = proteome_dict[protein]
        # starts are protein_dict keys
        protein_starts = sorted(protein_dict.keys())
        #  if there is only one region there is nothing to do
        if len(protein_starts) > 1:
            cur_start = protein_starts[0]
            cur_end = protein_dict[cur_start]

            for i in range(1, len(protein_starts)):
                prev_end = cur_end
                prev_start = cur_start
                cur_start = protein_starts[i]
                cur_end = protein_dict[cur_start]

                # if prev_start not in protein_dict:
                #     continue

                if cur_start == prev_end or cur_start == prev_end+1:
                    protein_dict[prev_start] = cur_end
                    del protein_dict[cur_start]
                    cur_start = prev_start



    return proteome_dict


def concat_iupred_res(path, path_to_out_file):
  """ These functions concat files from IUPred and made one file from them
  """
  all_files = glob.glob(path + "/*.result")

  # sort file names from _0 to _20000
  all_files.sort(key=lambda x: int(x.split('/')[-1].split('.')[0][5:]))
  with open(path_to_out_file, 'w') as outfile:
      for fname in all_files:
        with open(fname) as infile:
          for line in infile:
            outfile.write(line)
  outfile.write('# IUPred')


def make_iupred_df(human_df, full_proteome_dict):
  """ These functions add information about structured/unstructured region belonging 
  to the human_df according to the full_proteome_dict data, contained IUPred scores for each protein in human_df
  """
  df_proteins_list = []
  for protein in tqdm(full_proteome_dict.keys()):
        df_protein = clean_df[clean_df.ID==protein]
        for ind in df_protein.index:
            start = df_protein.loc[ind].Start
            stop = df_protein.loc[ind].Stop
            is_unstructured = False
            for possible_start in sorted(full_proteome_dict[protein].keys()):
                possible_stop = full_proteome_dict[protein][possible_start]
                if start >= possible_start and start < possible_stop:
                    if stop <= possible_stop:
                        is_unstructured = True
                    break
                else:
                    continue
            df_protein.at[ind, 'IUpred'] = 1 if is_unstructured else 0
        df_proteins_list.append(df_protein)

  human_iupred_df = pd.concat(df_proteins_list, axis=0, ignore_index=True)
  human_iupred_df.IUpred = human_iupred_df.IUpred.astype(np.int64)

  return human_iupred_df


def eval_compatible_arches(arc_score_path):
	with open(arc_score_path) as f:
	  arches = f.readline().split()[1:]
	  compatible_arches = {}
	  for line in f:
	    split_line = line.split()
	    i_arch = split_line[0]
	    compatible_arches[i_arch] = []
	    for i, score in enumerate(split_line[1:]):
	      if float(score) > 0.01:
	        compatible_arches[i_arch].append(arches[i])

	return compatible_arches


def find_amyloid_bound_parallel(humandf, sarsdf, compatible_arches, arch_type):
    new_cols = ['CompArch_sars_seq', 'CompArch_sars_start', 'HumanID',
                'CompArch_human_seq', 'CompArch_human_start',
                'CompScore', 'OtherScores']

    human_proteins = set()

    h_compatible = humandf[humandf.Arch.isin(compatible_arches[arch_type])]
    sars_arch_type = sarsdf[sarsdf.Arch == arch_type]
    sars_compscores_res = pd.DataFrame(columns=new_cols, index=sars_arch_type.index)

    # iters over the SARS arches of one type - arch_type
    for sars_arch_row in tqdm(sars_arch_type.itertuples()):
        sars_seqs = []
        sars_starts = []
        h_IDs = []
        h_seqs = []
        h_starts = []
        comp_scores = []
        other_scores = []
        sars_ba = ba.BetaArch(arch_type, start=1,
                              sequence=sars_arch_row.Sequence,
                              score=sars_arch_row.Score_total) #._16 means total score value (cause of space in the column name it replaced with number)
        # iters over Human arches, compatible with arch_type
        for h_arch_row in h_compatible.itertuples():
            h_protein = h_arch_row.ID
            h_ba = ba.BetaArch(h_arch_row.Arch, start=1,
                        sequence=h_arch_row.Sequence,
                        score=h_arch_row.Score_total)
            output = amp.CompArches(sars_ba, h_ba, MatrixInner, MatrixOuter, MatrixArc, MatrixInnerExtra)
            comp_score = output['Scores']['CompScore']
            if  comp_score >= 0.5: #0.15
                sars_comp_ba = output['ArchR']
                sars_seqs.append(sars_comp_ba.sequence)
                sars_starts.append(sars_comp_ba.start)
                h_comp_ba = output['ArchC']
                h_seqs.append(h_comp_ba.sequence)
                h_starts.append(h_comp_ba.start)
                h_IDs.append(h_protein)
                comp_scores.append(comp_score)
                other_scores.append(output['Scores'])
                human_proteins.add(h_protein)

        sars_compscores_res.loc[sars_arch_row.Index] = [sars_seqs, sars_starts,
                                                        h_IDs, h_seqs, h_starts,
                                                        comp_scores, other_scores]
    return sars_compscores_res, human_proteins

def find_amyloid_connection(sars_df_filt, iupred_df_filt, compatible_arches, res_path):
	for sars_protein in tqdm(sars_df_filt.ID.unique()):
	    sars_df_filt_protein = sars_df_filt[sars_df_filt.ID == sars_protein]
	    results = Parallel(n_jobs=8)(delayed(find_amyloid_bound_parallel)(iupred_df_filt, sars_df_filt_protein, compatible_arches, arch_type) \
	                                  for arch_type in sars_df_filt_protein.Arch.unique())
	    npresults = np.array(results, dtype=object)
	    pd_result = pd.concat(npresults[:,0], axis=0).sort_index()
	    sars_res_df = sars_df_filt_protein.merge(pd_result, left_index=True, right_index=True)

	    human_proteins = set()
	    for h_proteins in npresults[:, 1][0:]:
	       human_proteins = human_proteins.union(h_proteins)

	    # write results to files
	    sars_res_df.to_csv(res_path + sars_protein +'_compscore.csv')
	    with open(res_path+'human_proteins_'+sars_protein+'.txt', 'w') as outfile:
	        for protein in human_proteins:
	            outfile.write(protein+'\n')
	    with open(res_path + sars_protein +'_compscore.pickle', 'wb') as f:
	        # Pickle the 'data' dictionary using the highest protocol available.
	        pickle.dump(sars_res_df, f, pickle.HIGHEST_PROTOCOL)


def detect_pathogenic_amyloids(path_to_proteins_txts, amyloidase_base, additional_base):

    all_files = glob.glob(path_to_proteins_txts + "/*.txt")
    with open('amyloidase_res.txt', 'w') as amyl:
      for filename in all_files:
          sars_prot = filename[-10:-3]
          with open(filename) as f:
            for line in f:
              # print(line, sars_prot)
              if line[:-1] in amyloidase_base:
                amyl.write(f'{line[:-1]} is a {amyloidase_base[line[:-1]]} (main) and could be bound with {sars_prot} \n')
              if line[:-1] in additional_base:
                amyl.write(f'{line[:-1]} is a {additional_base[line[:-1]]} (additional) and could be bound with {sars_prot} \n')


