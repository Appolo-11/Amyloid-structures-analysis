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
huntingtin (initiator of Huntington chorea) or proteins which amyloids initiate Alzheimer desease. By the way, there is an experimental work that approved _in vitro_ 
that SARS-CoV-2 S-protein increase intesivity of $\alpha$-synuclein aggregation [1]. And since that our results could help to interpret mechanics of this process 
since we predict not only proteins, but the regions of possible co-aggregation. Identification of the pathogenic amyloids was done based on the published information [2].

## Literature comarison

We found recently published work [3], that experimetally approved several amyloid regions (20AA-peptides) of SARS-CoV-2 S-protein. We compare our foundations with these regions and observe nearly 500 human proteins, which were identified as possibly c0-aggregating with these regions. Among them were identified one known pathogenic amyloid - Corneodesmosin, which initiate some skin deseases.


## References

[1] Slav A. Semerdzhiev et al. Interactions between SARS-CoV-2 N-Protein and α-Synuclein Accelerate Amyloid Formation. CS Chem. Neurosci. 2022, 13, 1, 143–150
Publication Date:December 3, 2021 https://doi.org/10.1021/acschemneuro.1c00666

[2] Matiiv A.B. et al. AMYLOID AND AMYLOID-LIKE AGGREGATES: DIVERSITY AND THE TERM CRISIS. Biochemistry (Moscow). 2020. Т. 85. № 9. С. 1011-1034.  doi:10.31857/S0320972520090043

[3] Sofie Nyström, Per Hammarström. Amyloidogenesis of SARS-CoV-2 Spike Protein. bioRxiv 2021.12.16.472920; doi: https://doi.org/10.1101/2021.12.16.472920
