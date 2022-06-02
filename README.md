This code will run the bulk of the pipeline for identifying protein-protein interactions between a virus and host, filtering for host proteins that are involved in a given biological process, and analysis of variants in the genes encoding these host proteins.

# SCRIPTS
## analyze_all_variants.ipynb
This script will run most of the computations by calling the other scripts.

## tools.py
Contains a variety of tools used in other scripts, including parsing variant specifics from a string, finding the associated AA variant for an NT variant, aligning multiple sequences, etc.

## bio_tools.py
Contains a series of utility functions, specifically for parsing and calling Entrez resources.

## compute_features.py
Contains many utility functions and wrappers for dependencies.

## get_dbSNP.py
This script will query variants from dbSNP in the associated host genes, and then analyze them with a number of tools.

## parse_gwas.py
This script will filter out GWAS data to exclude irrelevant variants, then analyze all variants.

## compute_RSCU.py
Computed RSCU and RSCPU for codon and codon pair usage.

## stat_tools.py
Contains utility functions for statistical analyses.

## rate4site2.py
An attempt to recreate the overall function of rate4site, extended to other types of sequences (NTs, codons, etc).

## parse_lovd2.py
Parses LOVD for combinations of disease-associated variants from patients, given a disease/gene. Then, queries dbSNP for benign (or at least non-pathogenic) variants to generate comparable control combinations. 

## score_variants.py
Input is spreadsheet containing the variant combinations as rows and variants as columns, where 1 indicates variant contained in combination.
Utilizes gene spreadsheet, VEP, and other tools to compute variant scores, then combines multiple variants based on model.

## combine_documents.py
Used to simplify analyses by reducing dimensions.

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
Python libraries: pandas, Biopython, scikit-learn, numpy, scipy, numba, networkx, prody, BeautifulSoup

## RARE CODON ENRICHMENT CALCULATION (rc_enrichment)
Originally from Jacobs and Shakhnovich, 2017 (https://faculty.chemistry.harvard.edu/shakhnovich/software/coarse-grained-co-translational-folding-analysis), but optimized here to improve Poisson binomial calculation. 

## EVMUTATION (EVmutation)
Originally from Marks lab (Hopf et. al., 2017) (https://github.com/debbiemarkslab/EVmutation). Given here with a wrapper to realign with the focus sequence.
Requires PLMC installation.

# HOW TO USE
We strongly recommend having a disk mounted in RAM, as many intermediate programs use disk writes, which can reduce the lifespan of SSDs. This variable is specified in each .py file as "ram_disk".
Set the appropriate Entrez email and API key at the top of each Python script. 
Most computations are run from the analyze_all_variants Jupyter notebook. Set the variables appropriately at the top of this notebook before running. Before the last cell of the notebook, you will need to manually generate the set of structural homologues by generating structures and querying them in Dali.
