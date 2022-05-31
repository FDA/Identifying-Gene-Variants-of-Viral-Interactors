This software package performs the calculations described in the publication

WM Jacobs and EI Shakhnovich (2017),
``Evidence of evolutionary selection for cotranslational folding,''
Proceedings of the National Academy of Sciences 114 (43), pp. 11434-11439.

Package contents
----------------
Executable scripts for free-energy calculations:
  calc_consensus_contacts.py
  calc_absolute_energies.py
  calc_elongation.py

Executable scripts for codon-usage calculations:
  calc_codon_usage.py
  calc_rare_enrichment.py
  calc_profile_comparison.py

Supporting files:
  codons.py
  fe.py
  polymer.py
  substructures.py

Dependencies
------------
In addition to packages included in the standard Python distribution,
this software package depends on the following Python packages:
  numpy     [www.numpy.org]
  scipy	    [www.scipy.org]
  networkx  [networkx.github.io]
  Bio	    [biopython.org]
  prody	    [prody.csb.pitt.edu]

The calc_consensus_contacts.py script also uses Clustal Omega [www.clustal.org/omega].
It is assumed that this tool is installed in the usual place, i.e., /usr/bin/clustalo,
although this path can be changed in the code.

Usage
-----
These scripts were designed to be run with Python version 3+.
The usage and options for all primary scripts can be found using the --help option, e.g.

  python3 calc_consensus_contacts.py --help

The default values for all optional arguments, shown in brackets in the help text,
match those used in the main text of the publication.

Example free-energy calculation
-------------------------------
Assuming that the files 1DDR.pdb.gz (downloaded from www.rcsb.org)
and folA.fasta (downloaded, for example, from www.ecogene.org) are located in the
current working directory, run:

  python3 calc_consensus_contacts.py folA.fasta 1DDR
  	  # Outputs polymer.dat
  python3 calc_absolute_energies.py polymer.dat
  	  # Outputs polymer_abs_energies.dat
  python3 calc_elongation.py folA polymer_abs_energies.dat
  	  # Outputs folA_elongation_profile.dat

Free-energy and elongation calculations
---------------------------------------
- These scripts require crystal structures (in PDB format) for the proteins of interest, which can be found at the Protein Data Bank. They also require a fasta file containing the protein sequence.
- The calc_consensus_contacts script reads the crystal structures to compute amino acid pair bond energies, based on native contacts. It will only include such contacts if a minimum fraction of crystal structures contain the contacts.
- The calc_absolute_energies.py script computes energies based on the model given in the paper.
- The calc_elongation.py script computes free-energy estimates as the nascent amino acid chain elongates. It may take a while.


Codon-usage and rare-codon calculations
---------------------------------------
- These scripts require a nucleotide sequence fasta files (not necessarily aligned), which can be 
generated from a nucleotide BLAST search. They can be aligned with align_codons.py. Protein abundances for 
the genes of interest may be used as input, but are not necessary. They are available from PaxDB.
- The calc_codon_usage script should be run first to generate the genome-wide codon statistics.  
In general, it should be run again to get decile codon statistics.
- The calc_enrichment_script should then be run once for each gene. It may take a while.
- Finally, calc_profile_comparison.py should be run. It takes the elongation profile and rare codon profiles as input, assuming they have names like "gene_elongation_profile.dat" and "gene_rc_profile.dat".

Potential issues: 
- Frame shifts in the alignment. 
  - Solution: This should be solved by using align_codons.py to align instead of Clustal Omega. If not, use a different version of the sequence.
- Stop codons early on. 
  - Solution: Use a different version of the sequence.
- Very few sequences are used to compute rare codons.
  - Solution: Use a different version of the sequence.
- Some NT sequences contain amiguity characters (not A, T, U, G, C).
  - Solution: Remove those sequences.

Utility scripts
---------------------------------------
gene.sh will run all necessary scripts sequentially, given a comma separated list of the PDB IDs, the name of the gene of interest, and the directory containing everything. This is asssuming you have the PDB files, the nucleotide fasta file containing the gene of interest (gene_nt.fasta), and a nucleotide BLAST fasta file containing all similar sequences (gene_blast.fasta) in the directory.

find_stop_codons.py will search a nucleotide sequence in a fasta file for any stop codons.

clean_duplicates.py will remove sequences from a fasta file with the same name.

clean_fasta.py will remove the descriptions from the names of sequences in a fasta file, leaving only the accession numbers.

align_codons.py will convert the nucleotide sequences into amino acid sequences, align them, and then convert the aligned sequences back into nucleotides. It also removes sequences with the same name, and removes the description part of the sequence names, leaving only the NCBI accession id. It requires a path to the fasta file containing the sequences to align.

plot_polymer.py will plot the polymer model of the protein, which gives the amino-acid pair bond energies.

plot_energy.py will plot the cotranslational folding energy profile (elongation_profile.dat) output from calc_elongation.py. It also identifies regions where the energy drops. It requires a path to the elongation_profile.dat file, and an optional path to the output file.

plot_rare_enrichment.py will plot the rare codon enrichment (rc_profile.dat) over the length of the sequence. It uses both the fraction of sequence enriched with rare codons, and the probability that a sequence would be enriched on the same region in the null model.

plot_rc_and_energy.py combines the plots from plot_energy.py and plot_rare_enrichment.py into a single graph, aligning their x-axes.

compare_pdbs.py will plot the region where each pdb matches the sequence of interest. It requires a path to the fasta file for the sequence of interest, a comma-separated list of pdb ids, a path to the folder where the pdb.gz files are contained, and an optional path to the output file.

Note on documentation
---------------------------------------------
The code uses multiple graph models to represent proteins:
 - In one model, vertices represent amino acids, and edges represent native contacts.
 - In another model, vertices represent native contacts, and edges occur between native contacts involving neighboring amino acids. Connected components in this model can represent substructures, which frequently correspond to secondary structures.

Keep this in mind, as "vertex" may represent either an amino acid or a native contact, and "edge" may represent either a native contact or two native contacts involving neighboring vertices.
