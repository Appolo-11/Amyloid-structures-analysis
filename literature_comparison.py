import pandas as pd
import pickle
from bases import amyloidase_base, additional_base

def find_avg_arch_len(df):
    """ Find avg SARS-CoV-2 beta_arches length
    """
    for ind in df.index:
      count = 0
      sum_len = 0
      for arch in df.loc[ind].CompArch_sars_seq:
          sum_len += len(arch)
          count += 1

    avg_len = round(sum_len/count)
    return avg_len

def find_arches_in_spikes(sars_df):
    """ Take SARS-CoV-2 protein's dataframe with beta-arches information and compare 
    found beta_arches with known amyloid regions.
    Then choose SARS-CoV-2 beta-arches intersect that region at least on a half. 
    Return human proteins able to co-aggregate with filtered amount of beta-arches.
        
    """
    spikes191 = set() # amyloid 20AA peptide 191-210
    spikes532 = set() # amyloid 20AA peptide 532-551
    spikes599 = set() # amyloid 20AA peptide 599-618
    spikes1165 = set() # amyloid 20AA peptide 1165-1184

    AA_peptide_len = 20

    all_spikes = {191: spikes191, 532: spikes532, 599: spikes599, 1165: spikes1165}
    # amyloid_regions = {191: [198, 203], 532: [537, 544], 599: [606, 611], 1165: [1174, 1178]}
    for id in sars_df.index:
      for spike_start in all_spikes:
        for (sars_start, sars_seq) in zip(sars_df.loc[id].CompArch_sars_start, sars_df.loc[id].CompArch_sars_seq):
            region_start = sars_df.loc[id].Start
            loc_start = sars_start+region_start - 1
            loc_end = loc_start + len(sars_seq) 
            if loc_start < spike_start-avg_len/2 or loc_end > spike_start+AA_peptide_len-1 + avg_len/2:
                continue
            all_spikes[spike_start].add(id)
    return all_spikes

with open('res_path/P0DTC2_compscore.pickle', 'rb') as f:
    sars_TC2 = pickle.load(f)

arches_in_spikes_found = find_arches_in_spikes(sars_TC2)

article_human_protein_arches_in_spikes = set()
for s in arches_in_spikes_found:
  for id in arches_in_spikes_found[s]:
    article_human_protein_arches_in_spikes.update(sars_TC2.loc[id].HumanID)

with open('article_intersects_amyloids.txt') as f:
    for i in article_human_protein_arches_in_spikes:
      if i in amyloidase_base:
        f.write(f'{i} is a {amyloidase_base[i]} - main amyloidose\n')
      elif i in additional_base:
        f.write(f'{i} is a {additional_base[i]} - additional amyloidose\n')