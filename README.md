# Amyloid-structures-analysis
The set of programms for revealing possible SARS-CoV-2 proteins able to co-aggregate with human proteins.
_Remark:_ all amyloids of our interests are pathogenic amyloids. But we would name them just as amyloid.

## Aims and goals of project
**Aims:** To assess the probability of developing human amyloidosis as a result of coaggregation with SARS-CoV-2 proteins. 
To solve this problem, the AmyloComp program should be used, developed in SPbU laboratory together with the A.V. Kayava group (University of Montpellier, France)

**Goals:** 
- Database (db) preparation:
  - Creating a db of β-arches of the human proteome from existing data
  - Add to the db info about belonging to a structured site (by IUPred)
- Make a program predicting the probability of human proteins and SARS-CoV-2 proteins coaggregation based on 
available data with preliminary filtering by type of β-arch and the possibility to take into account (or not) belonging to an unstructured site 
- Enrichment analysis of GO terms in the obtained protein set
- Comparison of the obtained results with the literature data on the amyloid properties of the identified human proteins


## Results

As a result of work, we obtained the program, which are able to predict human proteins apotentially able to co-aggregate (i.e. construct amyloid structure) 
with SARS-CoV-2 proteins. The research was conducted based on the infromation about possible $\beta$-arches among human proteome and 8 known SARS-CoV-2 proteins.
That initial data are still not printed and could not be placed in a public repository. Since that we save structure of files and directories, 
but replace confidential files with empty or randomly filled ones.
We found several thousands potential candidates (human proteins) for each of the SARS-CoV-2 protein. But don't be such scared: the only information we use was information about β-arches. And our results should be experimentally approved.
In Table 1 presented the results of our search: two numbers in each column means different filtering bound - minimal score (probability) of $\beta$-arches existing. Both bounds were chosen according to the experiment testing. Upper bound is a median of scores for really existing β-arches and lower one as a 10%-quantile.

Table 1. Statstics of found human proteins for each SARS-CoV-2 protein
      ![image](https://user-images.githubusercontent.com/36989596/169540724-ee3a86ff-24a4-4b85-9e18-7ddaaf4064d7.png)

However, among 30 pathogenic variants we found some "interesting" candidates, such as $\alpha$-synuclein which amyloids provided Parkinson's desease, 
huntingtin (initiator of Huntington chorea) or proteins which amyloids initiate Alzheimer desease. By the way, there is an experimental work that approved _in vitro_ that SARS-CoV-2 S-protein increase intesivity of $\alpha$-synuclein aggregation [1]. And since that our results could help to interpret mechanics of this process since we predict not only proteins, but the regions of possible co-aggregation. Identification of the pathogenic amyloids was done based on the published information [2].

## Literature comparison

We found recently published work [3], that experimetally approved several amyloid regions (20AA-peptides) of SARS-CoV-2 S-protein. We compare our foundations with these regions and observe nearly 500 human proteins, which were identified as possibly c0-aggregating with these regions. Among them were identified one known pathogenic amyloid - Corneodesmosin, which initiate some skin deseases.


## References

[1] Slav A. Semerdzhiev et al. Interactions between SARS-CoV-2 N-Protein and α-Synuclein Accelerate Amyloid Formation. CS Chem. Neurosci. 2022, 13, 1, 143–150
Publication Date:December 3, 2021 https://doi.org/10.1021/acschemneuro.1c00666

[2] Matiiv A.B. et al. AMYLOID AND AMYLOID-LIKE AGGREGATES: DIVERSITY AND THE TERM CRISIS. Biochemistry (Moscow). 2020. Т. 85. № 9. С. 1011-1034.  doi:10.31857/S0320972520090043

[3] Sofie Nyström, Per Hammarström. Amyloidogenesis of SARS-CoV-2 Spike Protein. bioRxiv 2021.12.16.472920; doi: https://doi.org/10.1101/2021.12.16.472920

# Data and programm code description

_Remark 1:_ Unless otherwise specified, all functions are placed in sars_functions.py
_Remark 2_: Everything from `matrices` and `functions` directories are private data, since that we placed just "construction' and need organization, but the real program code or numerical values are removed or replaced with random data.

## 1. Input data
### Human data
Human data assumed to be a set of .csv files with the structure described below:

| ID  | Sequence | Arch  | Start | Score_Bstrand_length | Score_total | 
| ------------- | ------------- | ------------- | ------------- |------------- | ------------- |
| Q969G3  | NNYRLGG  | PPL  | 31 | 0.999  | 0.644  |
where:
- `ID` -- Protein UniProtID 
- `Sequence` -- β-arch sequence
- `Start` -- coordinate of $\beta$-arch start in the protein
- `Arch` -- $\beta$-arch type
- `Score_Bstrand_length` - β-arch length score
-  `Score_total` - score of β-arch from ArchCandy (probability of $\beta$-arch existence)

Function `concat_human_csv(path_to_csvs)` take path to the directory as an input and concat all .csv to the big one - we would name it `human_df`.

Then we need to embed that with the IUPred data (data taken from the [web-service](https://iupred2a.elte.hu/)) - that soft predict which regions of protein are structured/unstructured.
If your data as big as ours, you get from IUPred several big files. To concat them to the huge one, use `concat_iupred_res(path_to_iupred_results, path_to_out_file)` function.
Then we should transform it and integrate it to the `human_df`: this could be done with assistance of `read_web_iupred(path_to_cincat_iupred_results)` (returns full_proteome_dict) and `make_iupred_df(human_df, full_proteome_dict)`. Dataframe with IUPred results integrated would be named as `iupred_df`.

### SARS-CoV-2 data

 SARS-CoV-2 data assumed to be stored in one table, so we need no concatenation (if you do just use `concat_human_df` with columns you need). 
 The format of data should be same as initial human data. To read it use `read_sars_csv(path_to_sars_csv)`. The result database would be named as `sars_df`.
 
 ## 2. Other preparations
 
### Get compatible arches
 We need to create database where place compatible arches type. We do that based on the table of experimental compatibility scores. To do that and get the dictionary of compatible arches `compatible_arches` just run `eval_compatible_arches(path_to_Score_arc_5)`.
 
### Filter our data

We remember, that we have the scores-analogues of β-arch existence probability. We would like to use that and take the only data with probability higher than som bound. We use `Score_total > 0.55` bound and also we consider the only human β-arches, who placed in the unstructured region (since that fact increase the real probability to constrain β-arches structure). 
That is done by the next three command:
```
iupred_df = human_df[human_df.IUpred == 1]
iupred_df_filt = iupred_df[iupred_df['Score_total'] > 0.55]
sars_df_filt = sars_df[sars_df["Score_total"] > 0.55]
```

## 3. Searching for possible co-aggregations

The quintessence of that work is a program, predicted the probable co-aggregations of human proteins and SARS-CoV-2 ones. That could be done with `find_amyloid_connection(sars_df_filt, iupred_df_filt, compatible_arches, path_where_to_write_the_results)`. That functions create three files **for each SARS-Cov-2 protein**:

- <ProteinID>_compscore.pickle - binary file, but it allows to save lists and dictionaries in the table and then read them in the same format. 
- <ProteinID>_compscore.csv - format that is convenient to open in many tools and comfortable to look and analyze analyze by the gaze :)
- human_proteins_<ProteinID>.txt - list of human proteins, potentially possible to co-aggregate with SARs-CoV-2 protein

  The `.pickle` and `.csv` files took new columns, are presented below:
  
  					
| CompArch_sars_seq | CompArch_sars_start | HumanID  | CompArch_human_seq | CompArch_human_start | CompScore | OtherScores |
| ------------- | ------------- | ------------- | ------------- |------------- | ------------- | ------------- |
| NNYRLGG | 2 | [QN25W9, P034FT, ... ]  | [NNYR, LGG, ...] | [2, 5, ...]  | [0.55, 0.77, ... ]  | [other scores] |
  
  where:
  - CompArch_sars_seq --  SARS β-arch sequence taking active part in co-aggregation 

  - CompArch_sars_start -- CompArch_sars_seq start coordinate relatively ot the initial SARS β-arch

  - HumanID -- list of the human UniProtID's

  - CompArch_human_seq -- β-arch sequence taking active part in co-aggregation

  - CompArch_human_start -- CompArch_sars_seq start coordinate relatively ot the initial human β-arch

  - CompScore -- CompScore from CompArches

  - OtherScores --- other scores from CompArches
  
### Found known pathogenic amyloids 
  
  We go through the list of found proteins and detec amyloid ones among them with assistance of `detect_pathogenic_amyloids(res_path,amyloidase_base,additional_base)` - it write results in the `amyloidase_res.txt`.
  

## 4. Comparison with literature 
  
  Comparison with literature data are conducted with `literature_comparison.py`. Just run it after all previous steps and it write the pathogenic amyloids consists of β-arches from approved amyloidogenic regions into `article_intersects_amyloids.txt` file. 
  
## 5. Have a cup of your favorite drinks - you are awesome person, finished the pipeline! 
