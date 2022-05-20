import pandas as pd
import pandas as pd
import numpy as np
import glob
import json
from tqdm import tqdm


def concat_csv(path):
""" These functions concat all .csv files containde in path directory
"""	
	all_files = glob.glob(path + "/*.csv")

	# sort file names from _0 to _20000
	all_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))

	li = []

	for filename in all_files:
	    df = pd.read_csv(filename, index_col=None, header=0, quotechar="'")
	    df['ID'] = list(map(lambda x: x.split('|')[1], df['ID']))
	    if df.isnull().values.any():
	      print(f'{filename[-15:]} has missing values')
	    # df.Start = int(df.Start)
	    li.append(df)

	frame = pd.concat(li, axis=0, ignore_index=True)
	clean_df = frame.dropna()
	clean_df.Start=clean_df.Start.astype(np.int64)
	clean_df.Stop=clean_df.Stop.astype(np.int64)
	clean_df['IUpred'] = pd.Series(dtype=np.int64)

	return clean_df


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

    