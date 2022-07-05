This code will run the bulk of the pipeline for identifying protein-protein interactions between a virus and host, filtering for host proteins that are involved in a given biological process, and analysis of variants in the genes encoding these host proteins.
In addition, we've included some additional scripts from the Scoring-Multiple-Variants project, which have been updated to be compatible with newer utility functions.

# SCRIPTS
## analyze_all_variants.ipynb
This script will run most of the computations by calling the other scripts. Aside from processing intermediate data.

## get_dbSNP.py
This script will query variants from dbSNP in the associated host genes, and then analyze them with a number of tools. Includes functions for querying and scoring missense and synonymous variants from a desired gene.

## parse_gwas.py
This script will filter out GWAS data to exclude irrelevant variants, then analyze all variants. Currently, only splicing scores and data from dbSNP and ClinVar are processed for variants, as most variants found are usually intronic or in UTR's.

## tools.py
Contains a variety of tools used in other scripts, including parsing variant specifics from a string, finding the associated AA variant for an NT variant, aligning multiple sequences, etc. This script is widely used, as it contains a number of utility functions.

## bio_tools.py
Contains a series of utility functions, specifically for parsing and calling Entrez resources. Widely used, across multiple projects.

## compute_features.py
Contains many utility functions and wrappers for dependencies. Primarily used for calculation of variant features, such as conservation, mRNA scores, splicing, miRNA binding, etc.

## scores_to_csv.py
This script will compile all Rosetta score files in a specified directory and sort them by the specified column, return only the lowest five entries. This is useful for finding "optimal" protein structures.

## compute_RSCU.py
Computes relative synonymous codon usage (RSCU) and RSCPU for codon and codon pair usage. This is not frequently used, as the code is easy to reproduce.

## stat_tools.py
Contains utility functions for statistical analyses. Many useful functions for sampling, cross-validation, and others. This is primarily used for analysis of data and for projects involving maching learning predictors.

## rate4site2.py
An attempt to recreate the overall function of rate4site, extended to other types of sequences (NTs, codons, etc). This script is not complete, and is included only for availability to other projects.

## score_variants.py
Input is spreadsheet containing the variant combinations as rows and variants as columns, where 1 indicates variant contained in combination.
Utilizes gene spreadsheet, VEP, and other tools to compute variant scores, then combines multiple variants based on model. This script is used primarily in our multiple mutation scoring project.

# DEPENDENCIES 
Clustal Omega
NetSurfP2
netphos
netNglyc
netOglyc
netPhos
NUPACK
mFold
KineFold
remuRNA
ViennaFold (RNAfold)
PLMC
Python libraries: pandas, Biopython, scikit-learn, numpy, scipy, numba, networkx, prody, BeautifulSoup (see requirements.txt)
Full dependencies are included in requirements.txt.

## RARE CODON ENRICHMENT CALCULATION (rc_enrichment)
Originally from Jacobs and Shakhnovich, 2017 (https://faculty.chemistry.harvard.edu/shakhnovich/software/coarse-grained-co-translational-folding-analysis), but optimized here to improve Poisson binomial calculation. This calculated rare codon enrichment and conservation of rare codon enriched regions, as well as changes in protein folding energy that occur during translation.

## EVMUTATION (EVmutation)
Originally from Marks lab (Hopf et. al., 2017) (https://github.com/debbiemarkslab/EVmutation). Given here with a wrapper to realign with the focus sequence.
Requires PLMC installation.
This is a very powerful conservation calculation that incorporates not only position-specific amino acids but also position pairing (which may impact protein structure contacts.

# HOW TO USE
We strongly recommend having a disk mounted in RAM, as many intermediate programs use disk writes, which can reduce the lifespan of SSDs. This variable is specified in each .py file as "ram_disk".
Set the appropriate Entrez email and API key at the top of each Python script. 
Most computations are run from the analyze_all_variants Jupyter notebook. Set the variables appropriately at the top of this notebook before running. Before the last cell of the notebook, you will need to manually generate the set of structural homologues by generating structures and querying them in Dali.
