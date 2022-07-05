#!/usr/bin/python
#EVmutation uses approximately (1/486)^2*len(ntseq)^2 GB RAM for amino acid sequences, and (1/806)^2*len(ntseq)^2 for nt sequences

import random, subprocess, re, os, sys, os.path, warnings, math, time, argparse, mygene, warnings, inspect, string, traceback, requests, json, shutil, psutil, pickle, statistics, bio_tools, hgvs, copy, tools, datetime
import numpy as np
import pandas as pd
import align_codons
import hmmlearn.hmm
from collections import Counter
from Bio import SeqIO, pairwise2, Entrez, SearchIO
from Bio.Blast import NCBIWWW
from urllib.request import urlopen
import urllib.parse
from multiprocessing import Queue, Process, Manager, Value
try:
	from BeautifulSoup import BeautifulSoup
except ImportError:
	from bs4 import BeautifulSoup

curpath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, curpath + '/EVmutation')
sys.path.insert(0, curpath + '/rc_enrichment')
import create_analysis
import calc_codon_usage, calc_rare_enrichment

#update the following variables to match the locations in your path
pos = "ORF nt#        (A of ATG=1)" #column name of location of variant in coding sequence
gene = "GENE" #column name of gene 
nm = "NM"
nt_wt = "nt (WT)" #column name of wild type nucleotide
nt_mut = "nt (mut)" #column name of mutant nucleotide
remuRNA_cmd = "/media/home/software/remuRNA/remuRNA" #path to remuRNA executable
uniclust = "/media/home/software/hh-suite/data/uniclust30_2018_08" #path to uniclust for NetSurfP2
hhsuite_model = "/media/home/software/NetSurfP/models/hhsuite.pb"
entrez_email = "insert here"
Entrez.email = entrez_email
entrez_apikey = "insert here"
Entrez.api_key = entrez_apikey
ram_disk = "insert here"


#dict containing BLOSUM62 values for amino acid sequences
blosum = {'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': -2, 'Z': -1, 'X': 0, '*': -4}, \
	'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'B': -1, 'Z': 0, 'X': -1, '*': -4}, \
	'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3, 'B': 3, 'Z': 0, 'X': -1, '*': -4}, \
	'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1, '*': -4}, \
	'C': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1, 'B': -3, 'Z': -3, 'X': -2, '*': -4}, \
	'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2, 'B': 0, 'Z': 3, 'X': -1, '*': -4}, \
	'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1, '*': -4}, \
	'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3, 'B': -1, 'Z': -2, 'X': -1, '*': -4}, \
	'H': {'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3, 'B': 0, 'Z': 0, 'X': -1, '*': -4}, \
	'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3, 'B': -3, 'Z': -3, 'X': -1, '*': -4}, \
	'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1, 'B': -4, 'Z': -3, 'X': -1, '*': -4}, \
	'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 1, 'X': -1, '*': -4}, \
	'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1, 'B': -3, 'Z': -1, 'X': -1, '*': -4}, \
	'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1, 'B': -3, 'Z': -3, 'X': -1, '*': -4}, \
	'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2, 'B': -2, 'Z': -1, 'X': -2, '*': -4}, \
	'S': {'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 0, 'X': 0, '*': -4}, \
	'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0, 'B': -1, 'Z': -1, 'X': 0, '*': -4}, \
	'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3, 'B': -4, 'Z': -3, 'X': -2, '*': -4}, \
	'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1, 'B': -3, 'Z': -2, 'X': -1, '*': -4}, \
	'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4, 'B': -3, 'Z': -2, 'X': -1, '*': -4}, \
	'B': {'A': -2, 'R': -1, 'N': 3, 'D': 4, 'C': -3, 'Q': 0, 'E': 1, 'G': -1, 'H': 0, 'I': -3, 'L': -4, 'K': 0, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1, '*': -4}, \
	'Z': {'A': -1, 'R': 0, 'N': 0, 'D': 1, 'C': -3, 'Q': 3, 'E': 4, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1, '*': -4}, \
	'X': {'A': 0, 'R': -1, 'N': -1, 'D': -1, 'C': -2, 'Q': -1, 'E': -1, 'G': -1, 'H': -1, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -1, 'P': -2, 'S': 0, 'T': 0, 'W': -2, 'Y': -1, 'V': -1, 'B': -1, 'Z': -1, 'X': -1, '*': -4}, \
	'*': {'A': -4, 'R': -4, 'N': -4, 'D': -4, 'C': -4, 'Q': -4, 'E': -4, 'G': -4, 'H': -4, 'I': -4, 'L': -4, 'K': -4, 'M': -4, 'F': -4, 'P': -4, 'S': -4, 'T': -4, 'W': -4, 'Y': -4, 'V': -4, 'B': -4, 'Z': -4, 'X': -4, '*': 1}}

#codon to amino acid translation for standard genetic code
rev_translate = {'I' : ['ATT', 'ATC', 'ATA'],
               'L' : ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
               'V' : ['GTT', 'GTC', 'GTA', 'GTG'],
               'F' : ['TTT', 'TTC'],
               'M' : ['ATG'],
               'C' : ['TGT', 'TGC'],
               'A' : ['GCT', 'GCC', 'GCA', 'GCG'],
               'G' : ['GGT', 'GGC', 'GGA', 'GGG'],
               'P' : ['CCT', 'CCC', 'CCA', 'CCG'],
               'T' : ['ACT', 'ACC', 'ACA', 'ACG'],
               'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
               'Y' : ['TAT', 'TAC'],
               'W' : ['TGG'],
               'Q' : ['CAA', 'CAG'],
               'N' : ['AAT', 'AAC'],
               'H' : ['CAT', 'CAC'],
               'E' : ['GAA', 'GAG'],
               'D' : ['GAT', 'GAC'],
               'K' : ['AAA', 'AAG'],
               'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
               '*' : ['TAA', 'TAG', 'TGA']}
translate = {codon : aa for aa,v in rev_translate.items() for codon in v}

#three-letter and one-letter code correspondence for amino acids
aa_to_3aa = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'B': 'ASX', 'C': 'CYS', 'E': 'GLU', 'Q': 'GLN', 'Z': 'GLX', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}
triple_to_aa = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'ASX': 'B', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLX': 'Z', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'CSD': 'C', 'HYP': 'P', 'BMT': 'T', '5HP': 'E', 'ABA': 'A', 'AIB': 'A', 'CSW': 'C', 'OCS': 'C', 'DAL': 'A', 'DAR': 'R', 'DSG': 'N', 'DSP': 'D', 'DCY': 'C', 'DGL': 'E', 'DGN': 'Q', 'DHI': 'H', 'DIL': 'I', 'DIV': 'V', 'DLE': 'L', 'DLY': 'K', 'DPN': 'F', 'DPR': 'P', 'DSN': 'S', 'DTH': 'T', 'DTY': 'Y', 'DVA': 'V', 'CGU': 'E', 'KCX': 'K', 'LLP': 'K', 'CXM': 'M', 'FME': 'M', 'MLE': 'L', 'MVA': 'V', 'NLE': 'L', 'PTR': 'Y', 'ORN': 'A', 'SEP': 'S', 'TPO': 'T', 'PCA': 'E', 'SAR': 'G', 'CEA': 'C', 'CSO': 'C', 'CSS': 'C', 'CSX': 'C', 'CME': 'C', 'TYS': 'Y', 'TPQ': 'F', 'STY': 'Y'}

#amino acid data
aa_data = {
	'A' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 1, "Hydrophobicity": 1.8},
	'C' : {"Charge": 0, "Polar Neutral": 1, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": 2.5},
	'D' : {"Charge": -1, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 1, "Hydrophobicity": -3.5},
	'E' : {"Charge": -1, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": -3.5},
	'F' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 1, "Large": 1, "Small": 0, "Hydrophobicity": 2.8},
	'G' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 1, "Hydrophobicity": -0.4},
	'H' : {"Charge": 1, "Polar Neutral": 0, "Aromatic": 1, "Large": 1, "Small": 0, "Hydrophobicity": -3.2},
	'I' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": 4.5},
	'K' : {"Charge": 1, "Polar Neutral": 0, "Aromatic": 0, "Large": 1, "Small": 0, "Hydrophobicity": -3.9},
    'L' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": 3.8},
	'M' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 1, "Small": 0, "Hydrophobicity": 1.9},
	'N' : {"Charge": 0, "Polar Neutral": 1, "Aromatic": 0, "Large": 0, "Small": 1, "Hydrophobicity": -3.5},
	'P' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 1, "Hydrophobicity": -1.6},
	'Q' : {"Charge": 0, "Polar Neutral": 1, "Aromatic": 0, "Large": 1, "Small": 0, "Hydrophobicity": -3.5},
	'R' : {"Charge": 1, "Polar Neutral": 0, "Aromatic": 0, "Large": 1, "Small": 0, "Hydrophobicity": -4.5},
	'S' : {"Charge": 0, "Polar Neutral": 1, "Aromatic": 0, "Large": 0, "Small": 1, "Hydrophobicity": -0.8},
	'T' : {"Charge": 0, "Polar Neutral": 1, "Aromatic": 0, "Large": 0, "Small": 1, "Hydrophobicity": -0.7},
	'V' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": 4.2},
	'W' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 1, "Large": 1, "Small": 0, "Hydrophobicity": -0.9},
	'Y' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 1, "Large": 1, "Small": 0, "Hydrophobicity": -1.3}	
}

#Input: Location of spreadsheet containing SNP geneview output
def parse_SNP_geneview(df, sheet_name="Coding Seq", columns = {('Chr.', 'position'): 'Chr. pos', ('mRNA', 'pos'): 'mRNA pos', ('dbSNP rs#', 'cluster id'): 'dbSNP rsid', ('Hetero-', 'zygosity'): 'Heterozygosity', ('Clinical', 'Significance'): 'Clin sig', ('Protein', 'residue'): 'Amino acid', ('Codon', 'pos'): 'Codon pos', ('Amino acid', 'pos'): 'AA position'}, header=[0,1], skiprows=[0]):
	df = pd.read_excel(df, sheet_name=sheet_name, header=header, skiprows=skiprows)
	data = {}
	for i, row in df.iterrows():
		try:
			data[i] = {v:row[k] for k,v in columns.items()}
			data[i].update({"MAF":row[[x for x in df.columns if x[0] == "MAF"][0]], "Function":row[[x for x in df.columns if x[0] == "Function"][0]], "nt (mut)":row[[x for x in df.columns if x[0] == "dbSNP" and "allele" in x[1]][0]], "ORF nt#":3 * (int(row[('Amino acid', 'pos')])-1)+int(row[('Codon', 'pos')]), "PubMed":row[[x for x in df.columns if x[0] == "PubMed"][0]]})
		except Exception as e:
			print(e)
	for k,v in data.items():
		for k2, v2 in data.items():
			if v["ORF nt#"] == v2["ORF nt#"]:
				if not np.isnan(v2["Chr. pos"]):
					data[k]["Chr. pos"] = int(v2["Chr. pos"])
				if not np.isnan(v2["mRNA pos"]):
					data[k]["mRNA pos"] = int(v2["mRNA pos"])
	return(data)

#Input: Name of gene of interest, RefSeq accession ID, numeric gene ID, taxonomy ID of organism of interest (human by default)
#Output: SNP geneview of all SNVs in gene
def get_SNP_geneview(genename, nm="", geneid="", taxid=9606, funcfilt={}):
	if len(nm) < 1:
		accids = bio_tools.get_accids(genename)["mRNA"]
		nm = min(accids, key=lambda kv: int(re.findall("\_(\d+)\.?\d*", kv[0])[0]))[0]
	if not tools.is_numeric(geneid):
		#get the gene ID (works much better than gene name)
		geneid = str(0)
		compiled_data = {}
		mg = mygene.MyGeneInfo()
		hits = tools.retry_func(mg.query, [genename], {'size':10})["hits"]
		for hit in sorted(hits, key=lambda hit: hit["_score"], reverse=True):
			if hit["taxid"] == taxid:
				geneid = str(hit["_id"])
				break
	r = tools.retry_request(requests.get, ["https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?locusId=" + str(geneid)])
	trs = tools.get_element(r.text, k="TR")
	trs_match = []
	for tr in trs:
		nm_possible = re.findall("[A-Z]{2}\_\d+\.?\d*", tr)
		if any([True for i in range(len(nm_possible)) if bio_tools.acc_wo_version(nm) == bio_tools.acc_wo_version(nm_possible[i])]) and len(nm_possible) >= 3:
			try:
				nm_possible = [list({x for x in nm_possible if x.startswith("NM") or x.startswith("XM")})[0], list({x for x in nm_possible if x.startswith("XP") or x.startswith("NP")})[0], list({x for x in nm_possible if x.startswith("NT") or x.startswith("NW")})[0]]
				trs_match.append(nm_possible)
			except IndexError as e:
				pass
	def score_nm_possible(nm_pos_i):
		total = int(re.findall("\.(\d+)", nm_pos_i[0])[0])
		if nm_pos_i[0] == nm:
			total += 20
		if nm_pos_i[0][:2] == nm[:2]:
			total += 10
		if nm_pos_i[1][:2] == "NP":
			total += 5
		if nm_pos_i[2][:2] == "NT":
			total += 1
		return(total)
	trs_match = sorted(trs_match, key=lambda kv: score_nm_possible(kv), reverse=True)[0]
	nm = trs_match[0]
	r = tools.retry_request(requests.get, ["https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?showRare=on&chooseRs=coding&locusId=" + str(geneid) + "&mrna=" + trs_match[0] + "&ctg=" + trs_match[2] + "&prot=" + trs_match[1] + "&orien=forward"])
	table = [x for x in tools.get_element(r.text, k="TABLE") if "<BR>position </TH>\n<TH>mRNA<BR>pos  </TH>\n<TH> dbSNP rs#<BR>cluster id </TH>\n<TH> Hetero-<BR>zygosity </TH>\n<TH><a href=\"javascript:openPopup('/projects/SNP/Images/snp_legend_validation.jpg')\"> Validation</a></TH>\n<TH>MAF  </TH>\n<TH>Allele <br/>origin</TH>\n<TH>3D  </TH>\n<TH> <a href=\"\" STYLE=\"text-decoration:none; color:black\" onclick=\"this.href='javascript:void(0)';this.disabled=1\" title=\"With links to OMIM,LSDB or Clinical Lab\"> Clinically<br>Associated</a></TH>\n<TH>Clinical<br/>Significance</TH>\n<TH>Function  </TH>\n<TH> 	dbSNP<BR>allele </TH>\n<TH> Protein<BR>residue </TH>\n<TH> Codon<BR>pos </TH>\n<TH>Amino acid<BR>pos  </TH>\n<TH> PubMed </TH>" in x][0]
	headers = tools.get_element(tools.get_element(table, "THEAD")[0], "TH")
	headers = [re.sub("\<.+\>", " ", re.findall("\<TH\>(.+)\<\/TH\>", x)[0]).strip() for x in headers]
	rows = [tools.get_element(x, "td") for x in tools.get_element(tools.get_element(table, "TBODY")[0], "TR")]
	rows = [[re.sub("\<.+?\>", " ", y).strip() for y in row] for row in rows]
	data = {}
	for i, row in enumerate(rows):
		data[i+1] = {headers[j]:rows[i][j] for j in range(len(headers))}
		data[i+1][gene] = genename
		for k in ['Chr. position', 'mRNA pos']:
			v = data[i+1][k]
			if i in data.keys() and not any([True for k2 in ['Codon pos', 'Amino acid pos'] if data[i][k2] != data[i+1][k2]]):
				data[i+1][k] = data[i][k]
	data = verify_SNP_geneview(data, nm=nm, funcfilt=funcfilt)
	return(data)
	
#Input: data dictionary of SNP geneview, RefSeq accession ID
#Output: data dictionary with warnings for inconsistent/mismatching nucleotide and amino acid information and inconsistant clinical significance information
def verify_SNP_geneview(data, nm="", funcfilt={}):
	if len(nm) < 1:
		accids = bio_tools.get_accids(data[list(data.keys())[0]]["GENE"])["mRNA"]
		nm = min(accids, key=lambda kv: int(re.findall("\_(\d+)\.?\d*", kv[0])[0]))[0]
	transcript, ORF = bio_tools.nm_seqs(nm)
	offset = transcript.index(ORF)
	for k,v in list(data.items()):
		try:
			data[k].update({k2:int(data[k][k2]) for k2 in ['Chr. position', 'mRNA pos', 'Codon pos', 'Amino acid pos']})
			data[k][pos] = (v['Amino acid pos']-1)*3 + v['Codon pos']
			data[k][nt_wt] = ORF[data[k]["ORF nt#"]-1]
			data[k][nt_mut] = v["dbSNP allele"].strip()
			data[k]["index"] = str(v["GENE"]) + ":c." + str(data[k][pos]) + data[k][nt_wt] + data[k][nt_wt] + ">" + data[k][nt_mut]
			data[k]["warnings"] = []
			if offset + v["ORF nt#"] != v["mRNA pos"]:
				data[k]["warnings"].append("Position mismatch")
			if v["Function"] != "contig reference" and (v['dbSNP allele'] == transcript[v['mRNA pos']-1] or v['dbSNP allele'] == ORF[v["ORF nt#"]-1]):
				data[k]["warnings"].append("Allele mismatch")
			aa = re.findall("\[([A-Z])\]", v["Protein residue"])
			if re.match("^[A-Z]$", v["dbSNP allele"].strip()) and len(aa) > 0:
				mutseq = tools.update_str(ORF, v['dbSNP allele'].strip(), (v["ORF nt#"]-1))
				if translate_seq(mutseq)[v["Amino acid pos"]-1] != aa[0]:
					data[k]["warnings"].append("Wrong amino acid")
				if v["Function"] == "synonymous" and translate_seq(ORF)[v["Amino acid pos"]-1] != aa[0]:
					data[k]["warnings"].append("Not synonymous variant")
				if v["Function"] == "missense" and translate_seq(ORF)[v["Amino acid pos"]-1] == aa[0]:
					data[k]["warnings"].append("Not missense variant")
		except Exception as e:
			del data[k]
	data = {k:v for k,v in data.items() if all([v[k2].lower() in [x.lower() for x in v2] for k2, v2 in funcfilt.items()])}
	init_time = time.time()
	for i, (k,v) in enumerate(data.items()):
		try:
			if isinstance(v['dbSNP rs# cluster id'], str) and v['dbSNP rs# cluster id'].startswith("rs"):
				freqd, clinsig, nm, rsid = bio_tools.query_dbSNP(v['dbSNP rs# cluster id'])
				global_freqs = {k2:v2 for k2, v2 in freqd.items() if "Global" in k2 or "Total" in k2} #NOT population specific
				if len(global_freqs) > 0:
					data[k]["MAF"] = statistics.mean(global_freqs.values()) #may have multiple sources
				elif len(freqd) > 0:
					data[k]["MAF"] = statistics.mean(freqd.values())
				clinsig = list(set([x[1] if len(x) == 2 else x for x in clinsig]))
				if not any([all([y.lower() in x.lower() or y.lower() in x.replace(" ", "-").lower() for x in clinsig]) for y in funcfilt["Clinical Significance"]]): #clinical significances for entry match and are included in allowed clinical significances
					data[k]["warnings"].append("inconsistent clinical significance " + "/".join(clinsig))
		except Exception as e:
			pass
		finally:
			data[k]["warnings"] = ", ".join(data[k]["warnings"])
			tools.update_time(i, len(data), init_time)
	
	'''
	try:
		groups = {index:[k for k in data.keys() if data[k]["index"] == index] for index in set([data[k]["index"] for k in data])}
		for index, rows in groups.items():
			init_index = rows[0]
			if len(rows) > 1:
				for k in rows[1:]:
					for k2, v2 in data[k].items():
						if not (tools.is_numeric(v2) or tools.is_numeric(data[init_index][k2])):
							data[init_index][k2] = ", ".join(set([v2, data[init_index][k2]]))
					del data[k]
	except Exception as e:
		print(nm + "\t" + str(e))
	'''
	
	data = pd.DataFrame(data).T
	data = data.set_index("index")
	data = data.sort_values(["warnings", "Function", "Clinical Significance"])
	return(data)

#Input: multiple sequence alignment, name of focus sequence
#Output: distribution of characters at each position of focus sequence in MSA
def compute_distribution(seqs, genename, space="-", alphabet="", freq=True):
	joined_seqs = "".join(seqs.values()).replace(space, "")
	alphabet = sorted(set([alphabet[i] for i in range(len(alphabet))] + [joined_seqs[i] for i in range(len(joined_seqs))]))
	alphabetsize = len(alphabet)
	#print(seqs)
	inverted_seqs = {i:"".join([seqs[name][i] for name in seqs]) for i in range(len(seqs[genename]))}
	
	location = 1
	distribution = {}
	for i in range(len(seqs[genename])):
		if seqs[genename][i] != space:
			distribution[location] = {alpha:0 for alpha in alphabet}
			for j in range(len(inverted_seqs[i])):
				if inverted_seqs[i][j] != space:
					distribution[location][inverted_seqs[i][j]] += 1
			total = sum(distribution[location].values())
			if freq:
				distribution[location] = {k:v/float(total) for k,v in distribution[location].items()}
			location += 1
	return(distribution)

#Input: Command to run as a subprocess (externally)
#Output: Shell output of command
#Designed to avoid hanging processes
def getoutput_timer(cmd, t=60, interval=0.1):
	def getoutput(cmd, q):
		output = subprocess.getoutput(cmd)
		q.put(output)
	def waiter(t, q):
		for i in range(math.ceil(t/interval)):
			if not q.empty():
				return
			else:
				time.sleep(interval)
		return
	q = Queue()
	t1 = Process(target=getoutput, args=(cmd, q))
	t2 = Process(target=waiter, args=(t, q))
	t1.start()
	t2.start()
	t2.join()
	if not q.empty():
		output = q.get()
		t1.terminate()
		t2.terminate()
		return(output)
	else:
		t1.terminate()
		t2.terminate()
		return(None)

#Input: ORF nucleotide sequence (start codon to stop codon)
#Output: translated amino acid sequence
def translate_seq(seq):
	aaseq = ""
	for i in range(len(seq)//3):
		try:
			aaseq += translate[seq[3*i:3*(i+1)]]
		except KeyError as e:
			if seq[3*i:3*(i+1)] == "---":
				aaseq += "-"
			else:
				aaseq += "X"
	return(aaseq)
	#return("".join([translate[seq[3*i:3*(i+1)]] for i in range(len(seq)//3)]))

#Intput: queue, dictionary of genename: data
#Output: None
#Runs processes across multiple processes if possible
def work(q, data, ids={}, seqdir="./gene_data", path=ram_disk, outdir="./gene_data/", skip_done=False, pdb=None, space="-", quiet=False):
	while not q.empty():
		genename = q.get()
		seqs = tools.read_fasta(os.path.join(seqdir, genename + ".fasta"))
		try:
			geneids = ids[genename]
		except KeyError as e:
			print("Accession IDs not found for " + genename)
			geneids = {}
		if not skip_done or genename + ".tsv" not in os.listdir(seqdir):
			try:
				data[genename], geneids = precompute(genename, seqs, ids=geneids, path=path, outdir=outdir, space=space, quiet=quiet)
			except Exception as e:
				tb = traceback.format_exc()
				print("\033[93m" + str(tb) + "\033[0m")
				print("\033[93m" + str(e) + "\033[0m")
			finally:
				#remove temp files for current gene
				for f in os.listdir(path):
					if genename in f:
						if os.path.isfile(path + f):
							os.remove(path + f)
						elif os.path.isdir(path + f):
							shutil.rmtree(path + f)
			ids[genename] = geneids
	pd.DataFrame(ids).to_csv(outdir + "gene_info_updated.tsv", sep="\t")

#Input: nucleotide sequence
#Output: predicted mRNA folding energy
def run_mfold(seq, genename="pholder", t=60, path=ram_disk):
	pwd = os.getcwd()
	os.chdir(path)
	tools.write_fasta({genename:seq}, path + "tmp.fasta")	
	output = getoutput_timer("mfold SEQ='" + path + "tmp.fasta'", t=t)
	if "Segmentation fault (core dumped)" in output:
		raise RuntimeError("Segmentation fault while running mFold.")
	energy = re.findall("Minimum folding energy is (\-?\d+\.?\d*) kcal\/mol", output)
	os.chdir(pwd)
	if len(energy) > 0:
		return(float(energy[0]))
	else:
		return(None)

#Input: nucleotide sequence
#Output: predicted mRNA folding energy
def run_nupack(seq, genename="", t=60, path=ram_disk):
	with open(os.path.join(path, "tmp.in"), "w") as outf:
		outf.write(seq)
	output = getoutput_timer("mfe -pseudo " + os.path.join(path, "tmp") + " > " + os.path.join(path, "tmp.mfe"), t=t)
	if output != None:
		with open(os.path.join(path, "tmp.mfe"), "r") as inf:
			startflag = False
			lines = [line.strip() for line in inf.readlines()]
			for line in lines:
				if re.match("^\-?\d+\.?\d*$", line):
					if not startflag and int(float(line)) == len(seq):
						startflag = True
					else:
						return(float(line))

#Input: nucleotide sequence
#Output: predicted mRNA folding energy
#memory folders are preferred as path, as it writes multiple temporary files
def run_kinefold(seq, path=ram_disk, t=60, seed=None):
	if not isinstance(seed, int):
		seed = random.randint(1, 5000)
	#write sequence file
	with open(path + "tmp.dat", "w") as outf:
		outf.write("< tmp\n")
		outf.write(seq)
	#write configuration file
	with open(path + "tmp.req", "w") as outf:
		outf.write(str(seed) + "\n")
		outf.write(path + "tmp.p\n" + \
			path + "tmp.e\n" + \
			path + "tmp.rnm\n" + \
			path + "tmp.rnms\n" + \
			path + "tmp.rnml\n" + \
			path + "tmp.rnm2\n" + \
			path + "tmp.dat\n" + \
			"0		# 0=RNA ; 1=DN A\n" + \
			"6.3460741	# helix minimum free energy in kcal/mol: 6.3460741=10kT \n" + \
			"10000000	# NA \n" + \
			"1000		# folding time requested in msec \n" + \
			"1		# pseudoknots   1=yes 0=no \n" + \
			"0		# entanglements	1=yes 0=no \n" + \
			"2 3		# simulation type: 1=renaturation; 2 20=cotrans. @ 20msec/nt \n" + \
			"		# add T i j k or F i j k options here \n" + \
			"\n" + \
			"tmp\n" + \
			"tmp.zip\n" + \
			"<SEQNAME>tmp" + "_" + str(seed) + " \n" + \
			"<BASE>tmp\n" + \
			"<SEQUENCE>" + seq + "\n" + \
			"<ZIPFILE>tmp.zip")
	output = getoutput_timer("kinefold_long_static " + path + "tmp.req", t=t)
	if output == None or "Segmentation fault (core dumped)" in output:
		raise RuntimeError("Segmentation fault while running Kinefold.")

	#parse output
	m = re.findall("(-?\d+\.?\d*) kcal/mol", output)
	try:
		energy = float(m[-1])
	except IndexError as e:
		print(output)
		raise(e)
	return(energy)

#Input: nucleotide sequence
#Output: predicted mRNA folding energy
def get_RNAfold(seq, path=ram_disk, t=60):
	#write sequence file
	with open(path + "tmp.in", "w") as outf:
		outf.write(">tmp\n")
		outf.write(seq)
	output = getoutput_timer("RNAfold " + path + "tmp.in")
	m = re.findall("\(\s*(\-?\d+\.?\d*)\s*\)", output.split("\n")[-1])
	try:
		energy = float(m[-1])
	except IndexError as e:
		print(output)
		raise(e)
	return(energy)

#unfinished
'''
def get_seblastian(seq, genename, path=ram_disk, wait_time=60, give_up=12):
	r = tools.retry_request(requests.post, ["https://seblastian.crg.es/cgi-bin/seblastian_cgi.py"], {'data':{"routine": "Seblastian", "seblastian_mode": "knownSP", "seblastian_upstream_length": "5000", "seblastian_blast_evalue": "1e-3", "seblastian_max_distance": "3000", "search_complement": "on", "secis_filter": "secis_filter", "make_images": "make_images", "image_resolution": "150", "use_infernal": "infernal", "infernal_threshold": "10", "use_covels": "covels", "covels_threshold": 5, "use_secisearch1": "secisearch", "secisearch_pattern": "default", "sequence_filename": "(binary)", "sequence_fasta": seq}})
	html = re.findall("https\:\/\/seblastian\.crg\.es\/results\/\d+\/index\.html", r.text)
	if html != None and len(html) > 0:
		html = max(html, key=len)
	else:
		return(None)

	#wait until success or time out (defaults to 12 hours)
	def wait_for_seblastian(seq, genename, html, path=ram_disk, wait_time=60, give_up=12, quiet=True):
		retry_counter = 1
		success = False
		while retry_counter < give_up*60*(60/wait_time) and not success:
			results = tools.retry_request(requests.get, [html])
			if results != None and "FAILED" in results.text.upper(): #failure
				retry_counter += give_up*60*(60/wait_time)
				warnings.warn("\033[93m" + "Seblastian failed on " + genename + "\033[0m")
				return(False)
			elif results != None and results.status_code == 200 and "Results will be displayed here when they are ready." not in results.text: #success!
				results = results.text.split("<candidate_supertitle> SECIS prediction</candidate_supertitle>")[-1]
				scores = re.findall("Infernal\s\(score\s(\d+\.?\d*)\)\sCovels\s\(score\s(\d+\.?\d*)\)\sSECISearch1\.0\s\((.+)\)\<br\>", results)
				location = re.findall("\<candidate\_titles\>Positions\son\starget\:\<\/candidate\_titles\>\s(\d+\-\d+)", results)
				if location != None and len(location) > 0:
					location = [list(re.findall("\d+", location[i])) for i in range(len(location))]
				if scores == None or len(scores) < 1 or location == None or len(location) < 1:
					success = False
				if len(scores) != len(location) or any([True for i in range(len(scores)) if len(scores[i]) < 3]):
					warnings.warn("Unequal numbers of scores and locations read from Seblastian output")
				if len(scores[0]) >= 3 and len(location) >= 2:
					with open(path + genename + "_seblastian.out", "w") as outf:
						for i in range(len(scores)):
							outf.write("Infernal\t" + scores[i][0] + "\tCovels\t" + scores[i][1] + "\tSECISearch1.0\t" + scores[i][2] + "\n")
							outf.write("Location\t" + location[i][0] + "\t" + location[i][1])
					return(True)
			else:
				retry_counter += 1
				time.sleep(wait_time)
				elapsed = retry_counter*wait_time
				if not quiet:
					print("Waiting on Seblastian for " + str(elapsed/60) + " minutes, " + str(elapsed % 60) + " seconds" + " " * 20, end="\r")

	#start the process and return (we can do other stuff while it waits)
	p = Process(target=wait_for_seblastian, args=(seq, genename, html, path, wait_time, give_up, True))
	p.start()
	return(p)
'''
	
#unfinished
'''
def run_genesplicer(seqs, mut, directory="/media/home/software/GeneSplicer/human/", path=ram_disk):
	pos, nts = tools.get_mutation_data(mut)
	genomicpos, error = convert_position(seqs["ORF_aligned"], seqs["genomic_aligned"], pos)
	tools.write_fasta({"WT":seqs["genomic"]}, path + "tmp.fasta")
	wtoutput = getoutput_timer("genesplicer " + path + "tmp.fasta " + directory).split("\n")
	tools.check_mut((genomicpos, nts), seqs["genomic"], index=1)
	mutseq = tools.update_str(seqs["genomic"], nts[1], genomicpos+1)
	tools.write_fasta({"mut":mutseq}, path + "tmp.fasta")
	mutoutput = getoutput_timer("genesplicer " + path + "tmp.fasta " + directory).split("\n")
	wtscores = {}
	for line in wtoutput:
		pieces = line.strip().split()
		if len(pieces) >= 4:
			wtscores[(pieces[0], pieces[1])] = {"score":pieces[2], "confidence":pieces[3], "type":pieces[4]}
	mutscores = {}
	for line in mutoutput:
		pieces = line.strip().split()
		if len(pieces) >= 4:
			mutscores[(pieces[0], pieces[1])] = {"score":pieces[2], "confidence":pieces[3], "type":pieces[4]}
	
	dif = {"Gained":{k:v for k,v in mutscores.items() if k not in wtscores.keys()}, "Lost":{k:v for k,v in wtscores.items() if k not in mutscores.keys()}, "Changes":{k:(wtscores[k], mutscores[k]) for k in wtscores.keys() if k in mutscores and mutscores[k] != wtscores[k]}}
	return(dif)
'''

#Input: nucleotide ORF sequence, name of focus sequence
#Output: %minmax of sequence at each position, and %minmax of a random reverse-translated control sequence
#Note: this computes the %minmax on a network server and is very slow
def get_minmax(ORF, genename):
	headers = {"Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9", "Accept-Encoding": "gzip, deflate", "Accept-Language": "en-US,en;q=0.9", "Cache-Control": "max-age=0", "Connection": "keep-alive", "Content-Length": "2518", "Content-Type": "application/x-www-form-urlencoded", "Host": "www.codons.org", "Origin": "http://www.codons.org", "Referer": "http://www.codons.org/index.html", "Upgrade-Insecure-Requests": "1", "User-Agent": "Mozilla/5.0 (Windows NT 6.2; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.163 Safari/537.36"}
	data = {"Gene_Sequence": ">" + genename + "\n" + ORF, "Get_Codon_Usage": "CI-horfeome31.fa", "User_Codon_Usage_Text": "", "RRT_Check_Box": "On"}
	r = tools.retry_request(requests.post, positional_arguments=["http://www.codons.org/RCC.cgi"], keyword_arguments={'headers':headers, 'data':data})
	link = max(re.findall("\<a href\=\"(.+)\"\>DOWNLOAD RESULTS\<a\>", r.text), key=len)
	headers = {"Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9", "Accept-Encoding": "gzip, deflate", "Accept-Language": "en-US,en;q=0.9", "Connection": "keep-alive", "Host": "www.codons.org", "If-Modified-Since": "Wed, 08 Apr 2020 18:15:44 GMT", "If-None-Match": "\"3ebb-5a2cb7d12d400\"", "Referer": "http://www.codons.org/RCC.cgi", "Upgrade-Insecure-Requests": "1", "User-Agent": "Mozilla/5.0 (Windows NT 6.2; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.163 Safari/537.36"}
	r = tools.retry_request(requests.get, ["http://www.codons.org/" + link], {'headers':headers})
	if ">" + genename in r.text:
		lines = r.text.split("\n")
		lines = [l for l in lines if len(re.findall("\d+\.?\d*", l)) >= 3]
		data = {}
		for l in lines:
			nums = [float(i) for i in l.split(", ")]
			data[int(nums[0])+8] = {"%minmax": nums[1], "control %minmax": nums[2]}
	else: 
		raise Exception("Error when parsing output")
	newdata = {"%minmax":[], "%minmax control":[]}
	for k in range(len(ORF)//3):
		if k in data.keys():
			newdata["%minmax"].append(data[k]["%minmax"])
			newdata["%minmax control"].append(data[k]["control %minmax"])
		else:
			newdata["%minmax"].append("")
			newdata["%minmax control"].append("")
	return(newdata)

#Input: nucleotide ORF sequence, name of focus sequence, dictionary of codon:count data, dictionary of codon pair:codon pair data
#Output: %minmax of sequence at each position, and %minmax of a random reverse-translated control sequence
#Note: this computes the %minmax locally (on your computer) and doesn't use network calls
def minmax_local(ORF, genename, codonusage, cpusage=None, control=100):
	avgcodon = {aa:statistics.mean([codonusage[codon] for codon in rev_translate[aa]]) for aa in rev_translate.keys()}
	mincodon = {aa:min([codonusage[codon] for codon in rev_translate[aa]]) for aa in rev_translate.keys()}
	maxcodon = {aa:max([codonusage[codon] for codon in rev_translate[aa]]) for aa in rev_translate.keys()}
	if sum([codonusage[ORF[3*i:3*(i+1)]] for i in range(len(ORF)//3)]) > sum([avgcodon[translate[ORF[3*i:3*(i+1)]]] for i in range(len(ORF)//3)]):
		minmax = sum([codonusage[ORF[3*i:3*(i+1)]]-avgcodon[translate[ORF[3*i:3*(i+1)]]] for i in range(len(ORF)//3)])/sum([maxcodon[translate[ORF[3*i:3*(i+1)]]]-avgcodon[translate[ORF[3*i:3*(i+1)]]] for i in range(len(ORF)//3)])
	else:
		minmax = -sum([avgcodon[translate[ORF[3*i:3*(i+1)]]]-codonusage[ORF[3*i:3*(i+1)]] for i in range(len(ORF)//3)])/sum([avgcodon[translate[ORF[3*i:3*(i+1)]]]-mincodon[translate[ORF[3*i:3*(i+1)]]] for i in range(len(ORF)//3)])
	if control > 0:
		aaseq = translate_seq(ORF)
		dist_c = {aa:{c:codonusage[c] for c in codonusage if translate[c] == aa} for aa in rev_translate.keys()}
		dist_cp = {(c, aa):{cp[3:]:cpusage[cp] for cp in cpusage if cp[:3] == c and translate[cp[3:]] == aa} for c in translate.keys() for aa in rev_translate.keys()}
		for aa in dist_c.keys():
			total = sum(dist_c[aa].values())
			dist_c[aa] = {k:v/total for k,v in dist_c[aa].items()}
		control_minmax = {}
		for i in range(control):
			tmpseq = tools.custom_sample(dist_c[aaseq[0]])
			for j in range(1, len(aaseq)):
				tmpseq += tools.custom_sample(dist_cp[(tmpseq[3*(j-1):3*j], aaseq[j])])
			control_minmax[i] = minmax_local(tmpseq, genename + "_control", codonusage, control=0)
		return((minmax, statistics.mean(control_minmax.values())))
	else:
		return(minmax)

#Unfinished
'''
def minmax_optimized(ORF, codonusage, control=100, window=18):
	avgcodon = {aa:statistics.mean([codonusage[codon] for codon in rev_translate[aa]]) for aa in rev_translate.keys()}
	mincodon = {aa:min([codonusage[codon] for codon in rev_translate[aa]]) for aa in rev_translate.keys()}
	maxcodon = {aa:max([codonusage[codon] for codon in rev_translate[aa]]) for aa in rev_translate.keys()}

	
	translation = translate_seq(ORF)
	data = {0:{"seq":[ORF[3*i:3*(i+1)] for i in range(len(translation))], "%minmax":{}}}
	data.update({j:{"seq":[random.sample(rev_translate[translation[i]],1)[0] for i in range(len(translation))], "%minmax":{}} for j in range(1, 1+control)})
	averagesum = [avgcodon[translation[i]] for i in range(len(translation))]
	minsum = [mincodon[translation[i]] for i in range(len(translation))]
	maxsum = [maxcodon[translation[i]] for i in range(len(translation))]

	for i in range(window//2, len(translation)-(window//2)):
		averagesum_sub = sum([averagesum[j] for j in range(i-(window//2), i+(window//2))])
		minsum_sub = sum([minsum[j] for j in range(i-(window//2), i+(window//2))])
		maxnsum_sub = sum([maxsum[j] for j in range(i-(window//2), i+(window//2))])
		for j in data.keys():
			actualsum = sum([codonusage[data[j]["seq"][k]] for k in range(i-(window//2), i+(window//2))])
			if actualsum > averagesum_sub:
				data[j]["%minmax"][i] = (actualsum - averagesum_sub)/(maxnsum_sub - averagesum_sub)
			if actualsum <= averagesum_sub:
				data[j]["%minmax"][i] = (averagesum_sub - actualsum)/(averagesum_sub - minsum_sub)
	minmax = data[0]["%minmax"]
	minmax_control = {i:statistics.mean([data[j]["%minmax"][i] for j in data.keys() if j != 0]) for i in data[0]["%minmax"].keys()}
	return(minmax, minmax_control)
'''

#Input: nucleotide sequence
#Optional inputs: mutation string
#Output: predicted mRNA folding energy of WT sequence, mutant sequence, difference of mRNA folding energies, and H value
#memory folders are preferred as path, as it writes multiple temporary files
def run_remuRNA(seq, mut=None, path=ram_disk):
	if isinstance(mut, str):
		mut = tools.get_mutation_data(mut)
		mutstr = mut[1][0] + str(mut[0]) + str(mut[1][1])
	elif isinstance(mut, (list, tuple)) and tools.is_numeric(mut[0]) and len(mut[1]) == 2:
		mutstr = mut[1][0] + str(mut[0]) + str(mut[1][1])
	elif mut == None:
		mustr = seq[len(seq)//2] + str(len(seq)//2+1)
		if mutstr[0] == "A":
			mutstr += "C"
		else:
			mutstr += "A"
	with open(path + "tmp.fa", "w") as outf:
		outf.write(">tmp\n")
		outf.write(seq + "\n")
		outf.write("*" + mutstr)
	output = getoutput_timer(remuRNA_cmd + " " + path + "tmp.fa")
	retvals = (output.split("\n")[-1]).split()
	if len(retvals) < 6:
		raise ValueError("remuRNA encountered errors, or output incorrectly parsed for " + mut)
		print(output)
	return(float(retvals[1]), float(retvals[2]), float(retvals[4]), float(retvals[5]))

#Input: aligned ORF sequence (includes "-" gaps in introns and UTRs), aligned transcript (includes "-" gaps in introns)
#Output: mRNA MFE across ORF sequence (scaled between 0-1 values)
#computes mRNA sequences along sequence, scaling it to in [0,1]
def relative_mRNA_MFE(ORF_aligned, transcript_aligned, l=75, space="-"):
	ORF = ORF_aligned.replace(space, "")
	transcript = transcript_aligned.replace(space, "")
	energies = []
	tail = math.floor(l/2.0)
	init_time = time.time()
	for i in range(1, len(ORF)+1):
		j, error = convert_position(ORF_aligned, transcript_aligned, i)
		if not error:
			subseq = transcript[max(0,j-tail):min(len(transcript),j+tail+1)]
			energies.append(run_kinefold(subseq, seed=0))
		else:
			energies.append("")
		tools.update_time(i, len(ORF), init_time, func_name="mRNA MFE")
	maxenergy = max(energies, key=abs)
	minenergy = min(energies, key=abs)
	energies = [(energies[i] - minenergy)/(maxenergy - minenergy) for i in range(len(energies))]
	return(energies)

#Input: Any sequence, name of focus sequence, mode of operation ("amino acid" or "nucleotide", database to BLAST ("nt" or "nr"), program to use ("blastp" or "blastn"), whether to write to disk, any keywords for sequences to exclude
#Output: sequences returned from BLAST
#Note: if mode, db, and program are not specified, these are inferred from the sequence
def run_blast(seq, genename, mode=None, db=None, program=None, write=False, forbidden=None, path=ram_disk, space="-", size=2000, quiet=False, args={}):
	if mode == None:
		mode = tools.infer_alphabet(seq)
	elif mode in ["p", "protein", "amino acid", "aa"] and tools.infer_alphabet(seq) == "nt":
		try:
			seq = translate_seq(seq)
		except KeyError as e:
			pass
	if db == None:
		if mode in ["p", "protein", "amino acid", "aa"]:
			db = "nr"
		else:
			db = "nt"
	if program == None:
		if mode in ["p", "protein", "amino acid", "aa"]:
			program = "blastp"
		else:
			program = "blastn"
	instr = ">" + genename + "\n" + seq
	def try_blast(mode, db, instr, size):
		handle = NCBIWWW.qblast(mode, db, instr, hitlist_size=size)
		outstr = handle.read()
		handle.close()
		return(outstr)
	outstr = tools.retry_func(try_blast, [program, db, instr, size])
	if outstr != None:
		soup = BeautifulSoup(outstr, "lxml")
		hits = soup.find_all("hit")
		seqs = {}
		for hit in hits:
			accession = hit.find("hit_accession").get_text()
			name = hit.find("hit_def").get_text()
			if forbidden != None and ((isinstance(forbidden, (list, tuple,)) and not any([True for i in range(len(forbidden)) if forbidden[i] in name])) or (isinstance(forbidden, str) and forbidden in name)):
				#print(name)
				continue
			hsps = hit.find_all("hsp")
			for hsp in hsps:
				start = hsp.find("hsp_query-from").get_text()
				end = hsp.find("hsp_query-to").get_text()
				subseq = hsp.find("hsp_hseq").get_text()
				seqs[accession + ":" + start + "-" + end] = subseq
		seqs = tools.agglomerate_seqs(seqs)
		seqs[genename] = seq
		if not quiet: 
			print("Number of sequences returned from BLAST for " + genename + " with mode " + str(mode) + ": " + str(len(hits)) + (" "*20))
	else:
		seqs = {genename:seq}
	#write out results
	if write:
		tools.write_fasta(seqs, os.path.join(path, mode + "_msa.fasta"))
	return(seqs)

#Input: Any sequence, name of sequence, mode ("amino acid" or "nucleotide")
#Output: sequences returned from BLAST
#Note: full sequences are downloaded (individually), rather than simply hit regions
def run_blast_gb(seq, genename, mode=None, write=False, forbidden=None, path=ram_disk, space="-", size=2000, email=Entrez.email, apikey=Entrez.api_key, quiet=False):
	Entrez.email = email
	Entrez.api_key = apikey

	if mode in ["p", "protein", "amino acid", "aa"] and tools.infer_alphabet(seq) == "nt":
		try:
			seq = translate_seq(seq)
		except KeyError as e:
			pass
	inputs = run_blast(seq, genename, mode=mode, write=write, forbidden=forbidden, path=path, space=space, size=size)
	accids = list(inputs.keys())[:-1]
	seqs = get_accid_seqs(accids, mode=mode, email=email, apikey=apikey)
	seqs[genename] = seq
	if not quiet: 
		print("Number of sequences returned from BLAST with mode " + str(mode) + " for " + genename + ": " + str(len(hits)) + (" "*20))
	
	#write out results
	if write:
		tools.write_fasta(seqs, path + mode + "_msa.fasta")
	return(seqs)

#Input: List of RefSeq accessions, mode to search in ("protein" or "nucleotide"), Entrez email and key
#Output: sequences for each accession ID
def get_accid_seqs(accids, mode, email=Entrez.email, apikey=Entrez.api_key):
	Entrez.email = email
	Entrez.api_key = apikey
	seqs = {}
	while len(seqs) < 0.8*len(accids):
		try:
			for accid in accids:
				if accid not in seqs:
					if mode in ["p", "protein", "amino acid", "aa"]:
						handle = Entrez.efetch(db='protein',id=accid, rettype='fasta', retmode='text')
					else:
						handle = Entrez.efetch(db='nucleotide',id=accid, rettype='fasta', retmode='text')
					hseq = handle.read().split("\n")
					hseq = "".join([hseq[i].strip() for i in range(len(hseq)) if not hseq[i].startswith(">")])
					seqs[accid] = hseq
					time.sleep(1)
		except urllib.error.HTTPError as e:
			print(e.code)
			print(e.reason)
			if e.code == 429:
				time.sleep(5)
			else:
				break
	return(seqs)

#Input: multiple sequences (output from run_blast), name of focus sequence
#Output: multiple sequence alignment
def align_blast(seqs, genename, space="-", minlength=0.25, minsubset=0.5):
	'''
	seqs = {k:v.replace(space, "") for k,v in seqs}
	seqs = sorted(seqs.items(), key=lambda kv: len(kv[1]))
	seqs = {k:v for i, (k,v) in enumerate(seqs) if len(v) >= minlength*len(seqs[genename].replace(space, "")) or i >= minsubset*len(seqs) or k == genename}
	'''
	seqs = tools.align_sequences(seqs)
	seqs = {k:v for k,v in seqs.items() if "query" not in k}
	return(seqs)

#Input: sequences from BLAST, name of focus sequence
#Output: multiple sequence alignment, percent matching, entropy, variance, and also sum_of_pairs conservation scores
def compute_conservation(seqs, genename, space="-", mode="nt"):
	overall_p = {}
	for name, seq in seqs.items():
		counter = Counter(seq)
		for key in counter:
			if key != space:
				if key in overall_p:
					overall_p[key] += counter[key]
				else:
					overall_p[key] = counter[key]
	total = sum(overall_p.values())
	overall_p = {alpha:overall_p[alpha]/total for alpha in overall_p}
	alphabetsize = len(overall_p.keys())

	seqs = align_blast(seqs, genename)
	inverted_seqs = {i:"".join([seqs[name][i] for name in seqs]) for i in range(len(seqs[genename]))}
	
	identity = []
	entropy = []
	variance = []
	sum_of_pairs = []
	for i in range(len(seqs[genename])):
		if seqs[genename][i] != space:
			#compute percent identity
			pcnt = inverted_seqs[i].count(seqs[genename][i])/len(inverted_seqs[i].replace(space, ""))
			identity.append(pcnt)

			#compute entropy
			counts = Counter(inverted_seqs[i])
			del counts[space]
			alphabet = set(counter.keys()).difference(set([space]))
			p = {alpha:0 for alpha in alphabet}
			H = 0
			for i, alpha in enumerate(alphabet):
				p[alpha] = counts[alpha]/sum(counts.values())
				if p[alpha] != 0:
					H += p[alpha] * math.log(p[alpha]*alphabetsize)
			entropy.append(H)
		
			#compute variance
			sigma = 0
			for alpha in alphabet:
				sigma += (p[alpha]-overall_p[alpha])**2
			sigma = math.sqrt(sigma)
			variance.append(sigma)

			if mode in ["p", "protein", "amino acid", "aa"]:
				S = 0
				for alpha in alphabet:
					for beta in alphabet:
						if alpha in blosum and beta in blosum[alpha]:
							S += p[alpha]*p[beta]*blosum[alpha][beta]
						else:
							warnings.warn("\033[93m" + "Amino acids not found in BLOSUM dict: " + alpha + " " + beta + "\033[0m")
				sum_of_pairs.append(S)

	if mode in ["p", "protein", "amino acid", "aa"]:
		return(seqs, identity, entropy, variance, sum_of_pairs)
	else:
		return(seqs, identity, entropy, variance)	

#Input: position of variant in coding sequence, nucleotides changes (reference and mutant), and name of gene
#Output: prediction scores for variant, also MAFs, and any clinical significance from dbSNP
def get_vep(pos, nts, genename, epsilon=0.001):
	mutscores = {"Polyphen":[], "SIFT":[], "GnomAD":[], "dbSNP":None, "rsid":None}
	mutlist = [genename + ":c." + str(pos) + nts[0] + ">" + nts[1]]
	server = "https://rest.ensembl.org"
	ext = "/vep/human/hgvs"
	headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
	r = tools.retry_request(requests.post, positional_arguments=[server+ext], keyword_arguments={'headers':headers, 'data':'{ "hgvs_notations" : ' + str(mutlist).replace("\'", "\"") + ' }'})
	 
	if r != None:
		decoded = r.json()
		if "transcript_consequences" in decoded[0]:
			for j in range(len(decoded[0]["transcript_consequences"])):
				if "polyphen_score" in decoded[0]["transcript_consequences"][j]:
					mutscores["Polyphen"].append(1-decoded[0]["transcript_consequences"][j]["polyphen_score"])
				if "sift_score" in decoded[0]["transcript_consequences"][j]:
					mutscores["SIFT"].append(decoded[0]["transcript_consequences"][j]["sift_score"])
		if "colocated_variants" in decoded[0]:
			for j in range(len(decoded[0]["colocated_variants"])):
				if "frequencies" in decoded[0]["colocated_variants"][j]:
					try:
						mutscores["GnomAD"].append(decoded[0]["colocated_variants"][j]["frequencies"][nts[1]]["gnomad"])
					except KeyError as e:
						pass
				if "id" in decoded[0]["colocated_variants"][j] and "rs" in decoded[0]["colocated_variants"][j]["id"]:
					try:
						mutscores["rsid"] = decoded[0]["colocated_variants"][j]["id"]
					except:
						pass
		try:
			rsid = max(re.findall("\d+", mutscores["rsid"]))
			r = tools.retry_request(requests.get, positional_arguments=["https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/" + rsid])
			mutscores["dbSNP"] = score_variants.parse_dict(r.json(), "clinical_significances")
			mutscores["dbSNP"] = [score_variants.sigdict(val) for val in mutscores["dbSNP"]]
		except:
			mutscores["dbSNP"] = [0.5]

		for t in ["Polyphen","SIFT","GnomAD", "dbSNP"]:
			mutscores[t] = func_or_default(mutscores[t], func=statistics.mean, default=0.5, verbose=False)
		time.sleep(1)
		
		return(mutscores)

	return({"Polyphen":"", "SIFT":"", "GnomAD":"", "dbSNP":"", "rsid":""})

#Intput: input to a function, function to use, default alternate value
#Output: output of funciton in input, otherwise if errors are encountered, the default value
def func_or_default(l, func=statistics.mean, default=0.5, verbose=False):
	try:
		return(func(l))
	except Exception as e:
		if verbose:
			print(e)
		return(default)

#Input: Name of gene of interest, ORF sequence of gene, any Uniprot accession ID, taxid of organism of interest
#Output: Data from Uniprot entry that most closely matches the gene name and taxid, aligned to the sequence
#looks up a gene on Uniprot, gets the best match (for human genes), and reads through secondary structures, and binding sites, and post-translational modifications
def parse_uniprot(genename, ntseq, accession=None, taxid=9606):
	preferences = ["Domain", "Region", "Topological domain", "Transmembrane"]
	tmpdomain = {pref:{} for pref in preferences}
	data = {"Helix":[], "Beta strand":[], "Turn":[], "Binding site":[], "Metal binding":[], "Active site":[], "Glycosylation (N-linked)":[], "Glycosylation (O-linked)":[], "Disulfide bond":[], "Lipidation":[], "Cross-link":[], "Phosphorylation site":[], "Cleavage site":[], "Other modification":[], "Domain":{}}
	if accession == None or len(accession) < 1:
		response = tools.retry_request(requests.get, ["https://www.uniprot.org/uniprot/?query=" + genename + "+AND+organism:" + str(taxid) + "&columns=id,entry%20name,protein%20names,genes,organism,length&format=tab"])

		query = response.text
		accession = ""
		for line in query.split("\n"):
			if len(line.split("\t")) > 3 and genename in line.split("\t")[3].split():
				accession = line.split("\t")[0].strip()
				break
	if accession != "" and accession != None:
		response = tools.retry_request(requests.get, ["https://www.uniprot.org/uniprot/" + accession + ".fasta"])
		uniprotseq = "".join(response.text.split("\n")[1:])
		aaseq = translate_seq(ntseq)
		aligned_seqs = tools.align_sequences({"Uniprot":uniprotseq, "ORF":aaseq})
		response = tools.retry_request(requests.get, ["https://www.uniprot.org/uniprot/" + accession + ".gff"])
		table = response.text
		for line in table.split("\n"):
			pieces = line.split("\t")
			if len(pieces) >= 10:
				feat = pieces[2].strip()
				pos1, error1 = convert_position(aligned_seqs["Uniprot"], aligned_seqs["ORF"], pieces[3])
				pos2, error2 = convert_position(aligned_seqs["Uniprot"], aligned_seqs["ORF"], pieces[4])
				if not error1 and not error2:
					feat_lower = feat.lower()
					additional_lower = pieces[8].lower()
					r = [pos1, pos2]
					if "glycosylation" in feat_lower:
						if "o-linked" in additional_lower:
							data["Glycosylation (O-linked)"].append(r)
						elif "n-linked" in additional_lower:
							data["Glycosylation (N-linked)"].append(r)
					elif "modified residue" in feat_lower:
						if "phospho" in additional_lower:
							data["Phosphorylation site"].append(r)
						else:
							data["Other modification"].append(r)
					elif "disulfide bond" in feat_lower:
						data["Disulfide bond"].append(r)
					elif "active site" in feat_lower:
						data["Active site"].append(r)
					elif "metal binding" in feat_lower:
						data["Metal binding"].append(r)
					elif "site" in feat_lower and "cleavage" in additional_lower:
						data["Cleavage site"].append(r)
					elif "binding site" in feat_lower or ("site" in feat_lower and "binding" in additional_lower):
						data["Binding site"].append(r)
					elif "lipidation" in feat_lower:
						data["Lipidation"].append(r)
					elif "cross-link" in feat_lower:
						data["Cross-link"].append(r)
					if "domain" in feat_lower or "region" in feat_lower:
						try:
							name = pieces[8].split("=")[1].split(";")[0]
							for (feat2, oldname, r2) in sorted([(feat2, k, v) for feat2 in tmpdomain for k,v in tmpdomain[feat2].items()], key=lambda kv: kv[2][0]): #check if duplicate names
								if oldname == name:
									name = name + "_2"
								elif name in oldname and re.match(".+\_\d+", oldname):
									num = max(re.findall(".+\_(\d+)", oldname), key=len)
									name = name + "_" + str(int(num)+1)
							if feat in tmpdomain:
								tmpdomain[feat][name] = r
							else:
								tmpdomain[feat] = {name:r}
						except IndexError as e:
							print(e)
					elif feat in data.keys():
						data[feat].append(r)
		for feat in sorted(tmpdomain.keys(), key=lambda kv: preferences.index(kv) if kv in preferences else np.inf):
			for name1, r1 in sorted(tmpdomain[feat].items(), key=lambda kv: kv[1][0]):
				overlap = False
				for name2, r2 in sorted(data["Domain"].items(), key=lambda kv: kv[1][0]):
					if (r1[0] < r2[0] and r1[1] > r2[0]) or (r1[0] >= r2[0] and r1[1] <= r2[1]) or (r1[0] < r2[1] and r1[1] > r2[1]):
						overlap = True
				if not overlap:
					data["Domain"][name1] = r1					
	else:
		print("No matches found.")
	data["Accession"] = accession
	return(data)

#Input: a query statement for Protein Data Bank
#Output: the match score associated with each PDB/chain found
def query_pdb(query, quiet=False):
	url = 'https://search.rcsb.org/rcsbsearch/v1/query'
	req = tools.retry_request(requests.post, [url], {'data':json.dumps(query)})
	data = {}
	if req.status_code == 200:
		d = json.loads(req.text)["result_set"]
		for result in d:
			if "identifier" in result and len(result["identifier"].split(".")) > 1:
				pdb = result["identifier"].split(".")[0]
				chain = result["identifier"].split(".")[1]
				if pdb not in data.keys():
					data[pdb] = {chain:float(result["score"])}
				else:
					data[pdb][chain] = float(result["score"])
	return(data)

#Input: Uniprot accession ID, gene name, nucleotide ORF sequence
#Output: List of Protein Data Bank entries matching the query, in order of best to worst match
#searches PDBs based on a Uniprot accession, and return the best match based on sequence length, additional chains (proteins), and number of ligands
def parse_pdb(acc, genename, ORF="", quiet=False, seqsearch_min=10):
	aaseq = translate_seq(ORF)
	options = ["Uniprot", "Genename", "Sequence"]
	pdbs = {}
	for option in options:
		try:
			if option == "Uniprot":
				jsonval = {"query": {"type": "terminal","service": "text","parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession","operator": "exact_match","value": acc}},"request_options": {"return_all_hits": True, "sort": [{"sort_by": "rcsb_polymer_entity_feature_summary.maximum_length","direction": "desc"}]},"return_type": "polymer_instance"}
				tmppdbs = query_pdb(jsonval, quiet=quiet)
			elif option == "Genename":
				jsonval = {"query": {"type": "terminal","service": "text","parameters": {"attribute": "rcsb_entity_source_organism.rcsb_gene_name.value","operator": "exact_match","value": genename}},"request_options": {"return_all_hits": True, "sort": [{"sort_by": "rcsb_polymer_entity_feature_summary.maximum_length","direction": "desc"}]},"return_type": "polymer_instance"}
				tmppdbs = query_pdb(jsonval, quiet=quiet)
			elif option == "Sequence" and ORF != None and len(ORF) > 0:
				tmppdbs = {}
				evalue = 0.01
				while len(tmppdbs) < seqsearch_min and evalue <= 0.05: #increment evalue until match found
					pct_iden = 0.9
					while len(tmppdbs) < seqsearch_min and pct_iden >= 0.1: #decrement identity cutoff until match found
						aaseq = translate_seq(ORF).replace("*", "")
						jsonval = {"query": {"type": "terminal", "service": "sequence","parameters": {"evalue_cutoff": evalue, "identity_cutoff": pct_iden,"target": "pdb_protein_sequence","value": aaseq}},"request_options": {"return_all_hits": True, "scoring_strategy": "sequence"},"return_type": "polymer_instance"}
						tmppdbs = query_pdb(jsonval, quiet=quiet)
						pct_iden -= 0.05
					evalue += 0.01
				print("BLASTing against PDB: Min pct identity: " + str(round(pct_iden, 2)) + "\tMax e-value: " + str(round(evalue, 2)))
			for pdb in tmppdbs:
				if pdb not in pdbs:
					pdbs[pdb] = tmppdbs[pdb]
		except Exception as e:
			traceback.print_exc()
	pdbs = [(k, sorted(pdbs[k].keys(), key=lambda kv: pdbs[k][kv])) for k in pdbs]
	if not quiet:
		if len(pdbs) > 0:
			print("PDB selected for " + genename + ": " + pdbs[0][0] + ", with chain(s) " + ", ".join(pdbs[0][1]))
		else:
			print("No PDBs found for " + acc + " or " + genename)
	return(pdbs)

#Input: HTML from Dali query, name of gene, location to write file to
#Output: data from Dali query
def parse_dali(html, genename, path=ram_disk, write=True):
	if not html.endswith(".html") or html.endswith("/"):
		html = html + "s001A-25.html"
	matches = tools.retry_request(requests.get, [html])
	if matches == None or matches.status_code != 200:
		return(None)
	pdbs = re.findall("Sbjct\=\w{5}\sZ\-score\=\d+\.?\d*", matches.text)
	pdbs = {max(re.findall("Sbjct\=(\w{5})", i), key=len):float(max(re.findall("Z\-score\=(\d+\.?\d*)", i), key=len)) for i in pdbs}
	data = {}
	for k in pdbs:
		try:
			data[k] = max(re.findall(k[:-1].lower() + "\-" + k[-1].upper() + "\s+(\d+\.?\d*\s+\d+\.?\d*\s+\d+\s+\d+\s+\d+)", matches.text), key=len)
		except ValueError as e:
			pass
	newdata = {}
	for k,v in data.items():
		vstr = v.split()
		if len(vstr) >= 5:
			newdata[k[:-1]] = {"chain":k[-1], "Z-score":float(vstr[0]), "rmsd":float(vstr[1]), "lali":int(vstr[2]), "nres":int(vstr[3]), "%id":int(vstr[4])} 
	if write:
		df = pd.DataFrame.from_dict(newdata, dtype=str)
		df.T.to_csv(path + genename + "_dali.tsv", sep="\t")
	return(newdata)
	
#Input: Name of PDB entry (ABCDE, where E is the name of the chain and ABCD is the PDB accession), name of the gene)
#Output: Writes the output from dali to genename + "_dali.tsv"
#This gives structurally similar protein structures
def query_dali(pdb, genename, match=25, email=Entrez.email, path=ram_disk, space="-", wait_time=60, give_up=12):
	pdbstr = pdb[0].lower() + pdb[1][0].upper()
	headers = {"Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9", "Accept-Encoding": "gzip, deflate", "Accept-Language": "en-US,en;q=0.9", "Cache-Control": "max-age=0", "Connection": "keep-alive", "Content-Length": "703", "Content-Type": "multipart/form-data; boundary=----WebKitFormBoundaryOVadri1CnKBpmrIC", "Host": "ekhidna2.biocenter.helsinki.fi", "Origin": "http://ekhidna2.biocenter.helsinki.fi", "Referer": "http://ekhidna2.biocenter.helsinki.fi/dali/", "Upgrade-Insecure-Requests": "1", "User-Agent": "Mozilla/5.0 (Windows NT 6.2; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.163 Safari/537.36"}
	r = tools.retry_request(requests.post, ["http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/dump.cgi"], {'data':{'cd1':pdbstr, 'file1':'(binary)', 'title':genename, 'address':email, 'method':'search', 'submit':'Submit'}})

	#wait until success or time out (defaults to 12 hours)
	def wait_for_dali(pdbstr, genename, match=25, path=ram_disk, space="-", wait_time=60, give_up=12, quiet=True):
		retry_counter = 1
		success = False
		while retry_counter < give_up*60*(60/wait_time) and not success:
			results = tools.retry_request(requests.get, ["http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//" + pdbstr + "/"])
			if results != None and "FAILED" in results.text.upper(): #failure
				retry_counter += give_up*60*(60/wait_time)
				warnings.warn("\033[93m" + "Dali failed on " + genename + "\033[0m")
				return(False)
			if match in [25, 50, 90]:
				html = "http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//" + pdbstr + "/" + pdbstr + "-" + str(match) + ".html"
			else: 
				html = "http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//" + pdbstr + "/" + pdbstr + ".html"
			newdata = parse_dali(html, genename, write=True)
			if newdata != None: #success!
				success = True
			else:
				retry_counter += 1
				time.sleep(wait_time)
				elapsed = retry_counter*wait_time
				if not quiet:
					print("Waiting on Dali for " + str(elapsed/60) + " minutes, " + str(elapsed % 60) + " seconds" + " " * 20, end="\r")
		return(success)

	#start the process and return (we can do other stuff while it waits)
	p = Process(target=wait_for_dali, args=(pdbstr, genename, match, path, space, wait_time, give_up, True))
	p.start()
	return(p)
	

#Input: a list of objects
#Output: a set
#creates a list without duplicates while maintaining order or original list
def sorted_set(l):
	new_l = []
	for li in l:
		if li not in new_l:
			new_l.append(li)
	return(new_l)

#Input: Any dictionary of chain:index:data for PDB-specific data, the nucleotide ORF for the protein
#Output: data aligned to the protein sequence
#aligns data (keys are chains, next level of keys are positions, one of the included subkeys is aa), to ORF and generates new data
#generates data for each chain based on best alignment with sequence (we can't assume anything about chain ordering)
def get_aligned_data(data, ORF):
	wtseq = translate_seq(ORF)
	newdata = {}
	for chain in data:
		pdbseq = ""
		for i in range(min(data[chain].keys()), max(data[chain].keys())):
			if i in data[chain].keys():
				pdbseq += data[chain][i]["aa"]
			else: #between indices, mark with ambiguity
				pdbseq += "X"
		aligned_seqs = tools.align_sequences({"wtseq":wtseq, "pdbseq":pdbseq})
		for i in range(len(wtseq)):
			j, error = convert_position(aligned_seqs["wtseq"], aligned_seqs["pdbseq"], i) #position i in wtseq corresponds to j in pdbseq
			if j + min(data[chain].keys()) in data[chain] and wtseq[i] == data[chain][j + min(data[chain].keys())]["aa"] and not error: #position j in pdbseq corresponds to j+min(data[chain].keys() in the original data numbering
				newdata[i] = data[chain][j + min(data[chain].keys())]
	for i in range(len(ORF)//3): #fill out newdata with blanks
		if i not in newdata:
			newdata[i] = {key:"" for key in newdata[min(newdata.keys())]}
	return(newdata)	

#Input: name of best matching PDB entry and associated chains, nucleotide ORF sequence
#Output: conservation scores for each position in PDB entry for all chains given
def get_consurf(pdb, chains, ORF):
	data = {}
	prevpos = -1
	for chain in chains:
		data[chain] = {}
		try:
			response = requests.get("https://consurfdb.tau.ac.il/DB/" + pdb.upper() + chain.upper() + "/consurf_summary.txt")
		except requests.exceptions.SSLError as e:
			response = requests.get("https://consurfdb.tau.ac.il/DB/" + pdb.upper() + chain.upper() + "/consurf_summary.txt", verify=False)
		for i, line in enumerate(response.text.split("\n")):
			pieces = line.split("\t")
			if i >= 15 and len(pieces) >= 9:
				try:
					aa = pieces[1].strip()
					pos = int(max(re.findall(r'\d+', pieces[2].strip()), key=len)) #read amino acid position
					data[chain][pos] = {"aa": aa, "Conservation score":float(pieces[3].strip())}
					prevpos = pos
				except ValueError as e: #aa position was not parsed from line
					if prevpos >= 0: #previous line was associated with an aa position
						pos = prevpos + 1
						data[chain][pos] = {"aa": aa, "Conservation score":float(pieces[3].strip())}
						prevpos = pos
	newdata = get_aligned_data(data, ORF)
	scores = [newdata[i]["Conservation score"] for i in range(len(ORF)//3)]
	return(scores)

#Input: nucleotide ORF sequence, gene name, location to write to
#Output: conservation score and buried/exposed score at all positions in sequence
#Currently doesn't work, as Consurf has changed their website
def get_consurf_seq(ORF, genename, path=ram_disk, space="-", wait_time=60, give_up=12, quiet=True):
	aaseq = max(re.findall("[ACDEFGHIKLMNPQRSTVWY]+", translate_seq(ORF)), key=len) #maximal length protein sequence with only aa symbols
	response = tools.retry_request(requests.post, ["http://consurf.tau.ac.il/cgi-bin/consurf.cgi"], {"data":{"DNA_AA": "AA", "pdb_ID": "", "pdb_FILE": "", "PDB_chain": "", "PDB_yes_no": "yes", "modeller_key": "", "MSA_yes_no": "yes", "msa_FILE_ajax": "(binary)", "MAX_FILE_SIZE": "300000", "msa_SEQNAME_ajax": "", "msa_FILE": "(binary)", "msa_SEQNAME": "", "FASTA_text": genename + "\n" + aaseq, "Homolog_search_algorithm": "HMMER", "ITERATIONS": "1", "E_VALUE": "0.01", "proteins_DB": "UNIREF90", "user_select_seq": "no", "MAX_NUM_HOMOL": "250", "best_uniform_sequences": "uniform", "MAX_REDUNDANCY": 95, "MIN_IDENTITY": 35, "MSAprogram": "MAFFT", "tree_FILE": "(binary)", "ALGORITHM": "Bayes", "SUB_MATRIX": "BEST", "JOB_TITLE": genename, "user_email": "", "submitForm": "Submit"}})
	if response.status_code == 200: #response went through
		jobnum = max(re.findall("ConSurf is now processing your job number (\d+)", response.text), key=len)
	else: #failed for some reason
		warnings.warn("\033[93m" + "Consurf failed on " + genename + "\033[0m")
		return(False)
	#wait until we get a response
	print("Running Consurf on " + genename + ". Job num is " + jobnum)

	#wait until success or time out (defaults to 12 hours)
	def wait_for_consurf(jobnum, genename, path=ram_disk, space="-", wait_time=60, give_up=12, quiet=True):
		retry_counter = 1
		success = False
		while retry_counter < give_up*60*(60/wait_time) and not success:
			results = tools.retry_request(requests.get, ["http://consurf.tau.ac.il/results/" + jobnum + "/output.php"])
			if results != None and "FAILED" in results.text: #failure
				retry_counter += give_up*60*(60/wait_time)
				warnings.warn("\033[93m" + "Consurf failed on " + genename + "\033[0m")
				return(False)
			conservation_page = requests.get("http://consurf.tau.ac.il/results/" + jobnum + "/consurf.grades")
			aln_page = requests.get("http://consurf.tau.ac.il/results/" + jobnum + "/query_msa.aln")
			if conservation_page.status_code == 200 and aln_page.status_code == 200: #success!
				#get the MSA
				alnment = {}
				seqs = aln_page.text.split(">")
				for entry in seqs:
					try:
						seq = max(re.findall("[ACDEFGHIKLMNPQRSTVWY\\" + space + "\s]+", entry), key=len)
						name = entry.split(seq)[0]
						seq = re.sub("\s+", "", seq)
						alnment[name] = seq
					except ValueError as e:
						pass
		
				#read the conservation scores
				data = parse_consurf(conservation_page.text)
				
				success = True
			else:
				retry_counter += 1
				time.sleep(wait_time)
				elapsed = retry_counter*wait_time
				if not quiet:
					print("Waiting on Consurf for " + str(elapsed/60) + " minutes, " + str(elapsed % 60) + " seconds" + " " * 20, end="\r")
		if success: #write results to temp file
			df = pd.DataFrame.from_dict(data, dtype=str)
			df.T.to_csv(path + genename + "_consurf.tsv", sep="\t")
			tools.write_fasta(alnment, path + genename + "_consurf.fasta")
			return(True)
		else:
			return(False)

	#start the process and return (we can do other stuff while it waits)
	p = Process(target=wait_for_consurf, args=(jobnum, genename, path, space, wait_time, give_up, True))
	p.start()
	return(p)

#Input: text string from Consurf output (either from HTML or locally copied file)
#Output: conservation score and buried/exposed score from Consurf
def parse_consurf(text):
	data = {}
	for i, line in enumerate(text.split("\n")):
		pieces = re.split(" {2,}|\t+", line.strip())
		if i >= 16 and len(pieces) >= 11:
			try:
				score = float(pieces[3])
				b_e = ""
				if pieces[10].strip() == "b":
					b_e = 0
				elif pieces[10].strip() == "e":
					b_e = 1
				data[int(pieces[0])] = {"char":pieces[2], "score":score, "buried/exposed":b_e}
			except ValueError as e:
				pass
	return(data)

#queries Consurf based on the AA sequence, returns a process which waits for the data and eventually writes to file when Consurf is finished
def get_rnafold(seq, genename, path=ram_disk, space="-", wait_time=60, give_up=12, quiet=True):
	response = tools.retry_request(requests.post, ["http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAfold.cgi"], {"data":{"PAGE":2, "SCREEN":seq, "CONSTRAINT":"", "FILE": "(binary)", "method": "mfe", "noLP": "on", "dangling": "d2", "param": "rna2004", "SHAPEDATA": "", "SHAPEFILE": "(binary)", "shapemethod": "deigan", "shape_slope": "1.9", "shape_intercept": "-0.7", "shape_beta": "0.8", "deigan_conversion": "linearlog", "shape_conv_cutoff": "0.25", "shape_conv_linear_s": "0.68", "shape_conv_linear_i": "0.2", "shape_conv_linearlog_s": "1.6", "shape_conv_linearlog_i": "-2.29", "Temp": "37", "svg": "on", "reliability": "on", "mountain": "on", "EMAIL": "", "proceed": ""}})
	if response.status_code == 200: #response went through
		jobnum = max(re.findall("ConSurf is now processing your job number (\d+)", response.text), key=len)
	else: #failed for some reason
		warnings.warn("\033[93m" + "Consurf failed on " + genename + "\033[0m")
		return(False)
	#wait until we get a response
	print("Running Consurf on " + genename + ". Job num is " + jobnum)

	#wait until success or time out (defaults to 12 hours)
	def wait_for_consurf(jobnum, genename, path=ram_disk, space="-", wait_time=60, give_up=12, quiet=True):
		retry_counter = 1
		success = False
		while retry_counter < give_up*60*(60/wait_time) and not success:
			results = tools.retry_request(requests.get, ["http://consurf.tau.ac.il/results/" + jobnum + "/output.php"])
			if results != None and "FAILED" in results.text: #failure
				retry_counter += give_up*60*(60/wait_time)
				warnings.warn("\033[93m" + "Consurf failed on " + genename + "\033[0m")
				return(False)
			conservation_page = requests.get("http://consurf.tau.ac.il/results/" + jobnum + "/consurf.grades")
			aln_page = requests.get("http://consurf.tau.ac.il/results/" + jobnum + "/query_msa.aln")
			if conservation_page.status_code == 200 and aln_page.status_code == 200: #success!
				#get the MSA
				alnment = {}
				seqs = aln_page.text.split(">")
				for entry in seqs:
					try:
						seq = max(re.findall("[ACDEFGHIKLMNPQRSTVWY\\" + space + "\s]+", entry), key=len)
						name = entry.split(seq)[0]
						seq = re.sub("\s+", "", seq)
						alnment[name] = seq
					except ValueError as e:
						pass
		
				#read the conservation scores
				data = parse_consurf(conservation_page.text)
				
				success = True
			else:
				retry_counter += 1
				time.sleep(wait_time)
				elapsed = retry_counter*wait_time
				if not quiet:
					print("Waiting on Consurf for " + str(elapsed/60) + " minutes, " + str(elapsed % 60) + " seconds" + " " * 20, end="\r")
		if success: #write results to temp file
			df = pd.DataFrame.from_dict(data, dtype=str)
			df.T.to_csv(path + genename + "_consurf.tsv", sep="\t")
			tools.write_fasta(alnment, path + genename + "_consurf.fasta")
			return(True)
		else:
			return(False)

	#start the process and return (we can do other stuff while it waits)
	p = Process(target=wait_for_consurf, args=(jobnum, genename, path, space, wait_time, give_up, True))
	p.start()
	return(p)

#aligns the NT MSA file based on codons (codons align with ORF sequence), then computes standard conservation measures based on alignment
def codon_conservation(filename, genename, alphabet=["T","C", "A", "G", "t", "c", "a", "g"], path=ram_disk, space="-"):
	output = path + os.path.splitext(filename)[0] + "_aligned.fasta"
	align_codons.pre_align_codons(path+filename, genename, output = output, space=space)
	seqs = tools.read_fasta(output, aformat="WHOLE")
	for name, seq in list(seqs.items()):
		if not re.match("^["+"".join(alphabet)+"]+$", seq):
			del seqs[name]
	tools.write_fasta(seqs, output)
	align_codons.align_codons(output, output=output)
	seqs = tools.read_fasta(output, aformat="WHOLE")
	
	inverted_seqs = {i:{name:seqs[name][3*i:3*(i+1)] for name in seqs} for i in range(len(seqs[genename])//3)}
	
	#compute the overall distribution of codons in the alignment
	overall_p = {}
	for name, seq in seqs.items():
		for i in range(len(seq)//3):
			codon = seq[3*i:3*(i+1)]
			if codon != space * 3:
				if codon in overall_p:
					overall_p[codon] += 1
				else:
					overall_p[codon] = 1
	total = sum(overall_p.values())
	overall_p = {alpha:overall_p[alpha]/total for alpha in overall_p}
	alphabetsize = len(overall_p.keys())
	
	identity = []
	entropy = []
	variance = []
	numseqs = []
	for i in range(len(seqs[genename])//3):
		if space not in seqs[genename][3*i:3*(i+1)]:
			#compute percent identity
			pcnt = 0
			total = 0
			for name in seqs:
				if (inverted_seqs[i][name] == inverted_seqs[i][genename]):
					pcnt += 1
				if (space not in inverted_seqs[i][name]):
					total += 1
			pcnt = pcnt / total
			identity.append(pcnt)

			#compute entropy
			counts = {}
			for codon in set(inverted_seqs[i].values()):
				count = 0
				for name in seqs:
					if inverted_seqs[i][name] == codon:
						count += 1
				counts[codon] = count
			if space * 3 in counts.keys():
				del counts[space * 3]
			alphabet = set(counts.keys())
			p = {alpha:0 for alpha in alphabet}
			H = 0
			for i, alpha in enumerate(alphabet):
				p[alpha] = counts[alpha]/sum(counts.values())
				if p[alpha] != 0:
					H += p[alpha] * math.log(p[alpha]*alphabetsize)
			entropy.append(H)
		
			#compute variance
			sigma = 0
			for alpha in alphabet:
				sigma += (p[alpha]-overall_p[alpha])**2
			sigma = math.sqrt(sigma)
			variance.append(sigma)
			'''
			nonspace = []
			for gi in seqs:
				if space not in seqs[genename][3*i:3*(i+1)]:
					nonspace.append(gi)
			numseqs.append(len(nonspace))
			'''
	return(identity, entropy, variance)

#unfinished
'''
def mutation_taster(genename, pos, mutnt):
	response = tools.retry_request(requests.post, ["http://www.mutationtaster.org/cgi-bin/MutationTaster/MutationTaster69.cgi"], {"data":{"gene":genename, "ajax":1}})
	transcripts = re.findall("\"(ENST\d+)\"", response.text)

	response = None
	i = 0
	while not response and i < len(transcripts):
		response = tools.retry_request(requests.post, ["http://www.mutationtaster.org/cgi-bin/MutationTaster/MutationTaster69.cgi"], {"data":{"gene":genename, "transcript_stable_id_text":transcripts[i], "transcript_stable_id_radio":transcripts[i], "sequence_type":"CDS", "position_be":pos, "new_base":mutnt}})
		i += 1

	if "disease causing" in response.text:
		invert = 1
	elif "polymorphism" in response.text:
		invert = 0
	prob = float(re.findall("prob\: ([01]\.\d+)", response.text)[0])
	if invert:
		prob = 1 - prob
	return({"Pr Neutral":prob, "Disease":invert})
'''

#unfinished
'''
def get_PrDSM(genename, ntvars, PrDSM_loc="/media/home/software/PrDSM/PrDSM_Grch37_20190304.vcf.gz.tbi"):
	accid_to_chrom = {'NC_000001': '1', 'NC_000002': '2', 'NC_000003': '3', 'NC_000004': '4', 'NC_000005': '5', 'NC_000006': '6', 'NC_000007': '7', 'NC_000008': '8', 'NC_000009': '9', 'NC_000010': '10', 'NC_000011': '11', 'NC_000012': '12', 'NC_000013': '13', 'NC_000014': '14', 'NC_000015': '15', 'NC_000016': '16', 'NC_000017': '17', 'NC_000018': '18', 'NC_000019': '19', 'NC_000020': '20', 'NC_000021': '21', 'NC_000022': '22', 'NC_000023': 'X', 'NC_000024': 'Y'}
	seqs, cdses, mRNA, protein = bio_tools.pipeline(genename, outdir="./gene_data/", organism="homo sapiens", write=False)
	data = {}
	for (pos, nts) in ntvars.items():
		nts = tuple(nts)
		genvar = bio_tools.cds_to_genomic(mRNA, (pos, nts), assembly="GRCh37")
		version = int(mRNA.split(".")[-1])
		mRNA_wo_version = ".".join(mRNA.split(".")[:-1])
		while genvar == None:
			try:
				hp = hgvs.parser.Parser()
				var_c = hp.parse_hgvs_variant(mRNA_wo_version + "." + str(version) + ":c."+str(pos) + nts[0] + ">" + nts[1])
				hdp = hgvs.dataproviders.uta.connect()
				am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', alt_aln_method='splign', replace_reference=True)
				var_g = am.c_to_g(var_c, mRNA)
				genvar = str(var_g)
			except Exception as e:
				genvar = None
				print(e)
			finally:
				version -= 1
		genomic_acc = genvar.split(":g.")[0]
		location = genvar.split(":g.")[1]
		chromosome = accid_to_chrom[genomic_acc.split(".")[0]]
		gen_pos = int(max(re.findall("\d+", location), key=len))
		gen_nts = tuple(re.findall("[ACGT]", location))
		output = subprocess.getoutput("tabix "+PrDSM_loc+" "+str(chromosome)+":"+str(gen_pos-2)+"-"+str(gen_pos+2))
		lines = output.split("\n")
		lines = [line for line in lines if len(line.split("\t")) > 3 and line.split("\t")[3] == gen_nts[1] and line.split("\t")[1] == str(gen_pos)]
		if len(lines) < 1:
			data[(pos, nts)] = {"TraP":"", "SilVA":"", "FATHMM-MKL":"", "PrDSM":""}
		else:
			pieces = lines[0].split("\t")
			data[(pos, nts)] = {"TraP":pieces[4], "SilVA":pieces[5], "FATHMM-MKL":pieces[6], "PrDSM":pieces[7]}
	return(data)
'''

#generates a matrix of dimension len(WT sequence) x 3, with p(indel), p(miss), p(match) in each row
def hmm_seqs(seqs, genename, space="-", alphabet="-ACDEFGHIKLMNPQRSTVWY"):
	if alphabet == None or len(alphabet) < 1:
		alphabet = list(set(["".join(seqs.values())[i] for i in range(len("".join(seqs.values())))]))
	else:
		alphabet = [alphabet[i] for i in range(len(alphabet))]

	inverted_seqs = {i:"".join([seqs[name][i] for name in seqs]) for i in range(len(seqs[genename]))}
	distribution = []
	for i in sorted(inverted_seqs.keys()):
		if len(seqs[genename]) > i and seqs[genename][i] != space:
			indels = matches = inverted_seqs[i].count(space)/len(inverted_seqs[i])
			matches = inverted_seqs[i].count(seqs[genename][i])/len(inverted_seqs[i])
			misses = 1.0 - indels - matches
			distribution.append([indels, misses, matches])

	hmm = hmmlearn.hmm.GaussianHMM(n_components=len(alphabet)-1)
	hmm.fit(distribution)
	return(hmm)

#encodes any sequence for scoring with hmm, given the WT sequence
def encode_seq(seq, wtseq, space="-"):
	if len(seq) != len(wtseq):
		aligned_seqs = tools.align_sequences({"seq":seq.replace(space, ""), "wtseq":wtseq.replace(space, "")})
		seq = aligned_seqs["seq"]
		wtseq = aligned_seqs["wtseq"]
	seqlist = []
	for i in range(len(seq)):
		if len(wtseq) > i and wtseq[i] != space:
			for j in range(3):
				if seq[i] == space:
					seqlist.append([1,0,0])	
				elif seq[i] == wtseq[i]:
					seqlist.append([0,0,1])
				else:
					seqlist.append([0,1,0])
	return(seqlist)


#runs the Jacobs and Shakhnovich rare codon enrichment calculations
def rc_enrichment(seqs, genename, space="-", alphabet=["T","C", "A", "G", "t", "c", "a", "g"], path=ram_disk, quiet=True):
	for name, seq in list(seqs.items()):
		if not re.match("^["+"".join(alphabet)+"]+$", seq):
			del seqs[name]
	tools.write_fasta(seqs, path + "tmp.fasta")
	#convert the nucleotide sequences to amino acid sequences, align them, and convert them back to nucleotide sequences (pseudo-MSA)
	seqs = align_codons.align_codons(seqs)
	tools.write_fasta(seqs, path+"tmp.fasta")

	#rare codon conservation calculations
	try:
		os.remove(path+"input_codon_usage.p.gz")
		os.remove(path+"codon_usage.p.gz")
	except FileNotFoundError as e:
		pass	
	calc_codon_usage.calc_codon_usage(path+"tmp.fasta", output=path, wt_gi=genename)
	calc_codon_usage.calc_codon_usage(path+"tmp.fasta", output=path, wt_gi=genename)
	calc_rare_enrichment.calc_rare_enrichment(genename, path+"tmp.fasta", path+"codon_usage.p.gz", output=path, null_model="genome", wt_gi=genename)
	
	data = {i:"" for i in range(len(seqs[genename].replace(space, ""))//3)}
	with open(path+genename+"_rc_profile.dat", "r") as inf:
		for i, line in enumerate(inf):
			if i >= 3:
				pieces = line.strip().split()
				data[int(pieces[0])] = float(pieces[7])
	return([data[i] for i in range(len(data))])

#given a PDB ID and associated chains, gets and parses PDBPisa data for PDB and chains, then aligns with sequence (to get correct positions)
def get_PDBPisa(pdb, chains, ORF):
	response = tools.retry_request(requests.get, ["https://www.ebi.ac.uk/pdbe/api/pisa/interface/" + pdb.lower() + "/0"])
	d = json.loads(response.text)
	data = {}
	for interface in d[pdb.lower()]["interface"]["interface_array"]:
		interface = interface["interface"]
		chain = interface["molecule"]["chain_id"]
		if chain in chains:
			if chain not in data:
				data[chain] = {}
			for residue in interface["molecule"]["residues"]:
				try:
					position = int(residue["seq_num"])
					if position not in data[chain]:
						data[chain][position] = {"aa" : triple_to_aa[residue["name"]], "bonds" : residue["bonds"], "accessible_surface_area" : residue["accessible_surface_area"], "buried_surface_area" : residue["buried_surface_area"], "solvation_energy" : residue["solvation_energy"]}
					else: #update for most extreme values
						bonds = data[chain][position]["bonds"] + residue["bonds"] #combine both sets of bonds
						data[chain][position]["bonds"] = "".join(sorted(set([bonds[i] for i in range(len(bonds))])))
						data[chain][position]["accessible_surface_area"] = min(data[chain][position]["accessible_surface_area"], residue["accessible_surface_area"])
						data[chain][position]["buried_surface_area"] = max(data[chain][position]["buried_surface_area"], residue["buried_surface_area"])
						data[chain][position]["solvation_energy"] = max(data[chain][position]["solvation_energy"], residue["solvation_energy"], key=abs)
				except KeyError as e: #non-standard amino acids don't have a 3-letter and 1-letter code
					pass
	
	newdata = get_aligned_data(data, ORF)
	newdata = {key:[newdata[i][key] for i in range(len(ORF)//3)] for key in newdata[min(newdata.keys())]}
	for bond, b in [("Hydrogen bond", "H"), ("Salt bridge", "S"), ("Covalent bond", "C"), ("Disulfide bond", "D")]:
		newdata[bond] = [1 if b in newdata["bonds"][i] else 0 for i in range(len(newdata["bonds"]))]
	return(newdata["accessible_surface_area"], newdata["buried_surface_area"], newdata["solvation_energy"], newdata["bonds"])

def lookup_geneid(genename, taxid=9606):
	mg = mygene.MyGeneInfo()
	hits = tools.retry_func(mg.query, [genename], {'size':10})["hits"]
	for hit in sorted(hits, key=lambda hit: hit["_score"], reverse=True):
		if hit["taxid"] == taxid:
			return(hit["_id"])
			

#parses miRDB for genename
def get_miRDB(genename, genomic_aligned, ORF_aligned, taxid=9606, space="-", quiet=False):
	ORF = ORF_aligned.replace(space, "")
	hitstart = "<fontcolor=#0000FF><b>"
	hitend = "</font></b>"

	#get the gene ID (works much better than gene name)
	geneid = 0
	compiled_data = {}
	geneid = lookup_geneid(genename, taxid=taxid)

	if not quiet:
		print("Gene ID found: " + str(geneid).strip())
	response = tools.retry_request(requests.get, ["http://mirdb.org/cgi-bin/search.cgi?searchType=gene&geneChoice=geneID&searchBox=" + str(geneid) + "&full=1"])
	matches = re.findall("/cgi\-bin/target\_detail\.cgi\?targetID\=\d+", response.text) #search for links to matching miRNAs
	scores = re.findall("\<td align\=\"left\" width\=\"65\"\>\<p align\=\"center\"\>\<font size\=\"2\"\>(\d+)\</font\>\</td\>", response.text) #get the scores
	if len(matches) != len(scores):
		warnings.warn("\033[93m" + "Number of miRNAs does not match number of scores. Possibly incorrectly parsed from HTML" + "\033[0m")
	init_time = time.time()
	for j, match in enumerate(matches):
		accession = max(re.findall("\d+", match), key=len)
		response = tools.retry_request(requests.get, ["http://mirdb.org" + match])
		seqparts = re.findall("\d+\s+[gatc|\s|(?:\<font color\=\#0000FF\>\<b\>)(?:\</font\>\</b\>)]+", response.text) #get the sequence part, including formatted matches
		seqparts = [re.sub("\s+", "", re.sub("^\d+\s+", "", seqparts[i])) for i in range(len(seqparts))] 
		seq2parse = "".join(seqparts)
		
		seq = ""
		increment = 0
		i = 1
		data = {}
		bindsite = False
		#read through the sequence part, finding binding sites between html tags (hitstart and hitend)
		while increment < len(seq2parse):
			if increment+len(hitstart) <= len(seq2parse) and seq2parse[increment:increment+len(hitstart)] == hitstart: #binding site starts
				increment += len(hitstart)
				bindsite = True
			elif increment+len(hitend) <= len(seq2parse) and seq2parse[increment:increment+len(hitend)] == hitend: #binding site ends
				increment += len(hitend)
				bindsite = False
			elif seq2parse[increment] in ["a", "c", "g", "t"]: #actual sequence
				seq += seq2parse[increment]
				if bindsite:
					if i not in data or int(scores[j]) > data[i]:
						data[i] = int(scores[j])
				i += 1
				increment += 1
			else: #something else
				increment += 1

		#align sequences and convert to genomic gene positions
		aligned_seqs = tools.align_sequences({"wtseq": genomic_aligned.replace(space, ""), "mirdbseq": seq.upper()}) #align miRDB sequence with ours to convert positions
		for i in data:
			new_i, error = convert_position(aligned_seqs["mirdbseq"], aligned_seqs["wtseq"], i)
			if new_i not in compiled_data or compiled_data[new_i] < data[i]:
				compiled_data[new_i] = data[i]
		tools.update_time(j, len(matches), init_time, "miRDB")

	#convert positions to ORF positions
	data = []
	for i in range(1, len(ORF)+1):
		new_i, error = convert_position(ORF_aligned, genomic_aligned, i)
		if new_i in compiled_data and not error:	
			data.append(compiled_data[new_i])
		else:
			data.append(0)
	return(geneid, [data[i] for i in range(len(ORF))])

#parses miRDB for sequence (partially functional)
def custom_miRDB(genomic_aligned, ORF_aligned="", space="-", quiet=False):
	seq = genomic_aligned.replace(space, "")
	response = requests.post("http://mirdb.org/cgi-bin/custom.cgi", data={"customSub":seq, "subChoice":"mRNATarget", "searchSpecies":"hsa"})
	tmplink = max(re.findall("\<input type=\"hidden\" name=\"fileName\" value=\"(\S+)\".*/\>", response.text), key=len)
	response = requests.post("http://mirdb.org/cgi-bin/custom.cgi", data={"fileName":tmplink,".submit":"Retrieve Prediction Results"})
	matches = set(re.findall("\<input type=\"hidden\" name=\"(" + tmplink + "\_\d+)\" value=\"\S+\".*/\>", response.text))
	data = {i:0 for i in range(len(ORF_aligned.replace("-", "")))}
	init_time = time.time()
	for alpha, match in enumerate(matches):
		
		score = max(re.findall("\<input type=\"hidden\" name=\"" + match + "\" value=\"(\d+)\".*\>", response.text), key=len)
		mirna = max(re.findall("\<input type=\"hidden\" name=\"" + match + "\" value=\"(hsa\-miR\-\S+)\".*\>", response.text), key=len)
		subseq = max(re.findall("\<input type=\"hidden\" name=\"" + match + "\" value=\"([atgc]+)\".*\>", response.text), key=len)
		mirna_seq = max(re.findall("\<input type=\"hidden\" name=\"" + match + "\" value=\"([agcu]+)\".*\>", response.text), key=len)
		
		files = [(match, int(score)), (match, mirna), (match, "submission"), (match, "submission"), (match, "submission_1"), (match, subseq), (match, mirna_seq), (match, ""), (".submit","Details")]
		response2 = requests.post("http://mirdb.org/cgi-bin/custom_predict/customDetail.cgi", files=files, headers={"Content-Type": "multipart/form-data"})
		
		mirna_seq = mirna_seq.replace("t", "u").replace("a", "t").replace("u", "a").replace("c", "x").replace("g", "c").replace("x", "g").upper()[1:8] #complement
		subseq = subseq.upper()
		score = int(score)
		aligned_seqs = tools.align_sequences({"orig":seq, "subseq":subseq})
		for match in re.finditer(mirna_seq + "|" + mirna_seq[::-1], subseq, re.I):
			for i in range(match.start(), match.end()):
				j, error = convert_position(aligned_seqs["subseq"], aligned_seqs["orig"], i)
				if not error:
					k, error = convert_position(genomic_aligned, ORF_aligned, j)
					if k in data and score > data[k] and not error:
						data[k] = score
		tools.update_time(alpha, len(matches), init_time)
		
	return([data[i] for i in range(len(ORF_aligned.replace("-", "")))])

#searches ESEfinder for a genomic sequence, then converts positions for ORF
def run_ESEfinder(subseq="", orig_seq="", genename="1", space="-", quiet=False):
	seq = subseq.replace(space, "")
	data = {i:0 for i in range(1, 1+len(seq))}
	response = tools.retry_request(requests.post, ["http://krainer01.cshl.edu/cgi-bin/tools/ESE3/esefinder.cgi"], {"data":{"process":"search", "db":"SRProteins", "check_sf2":"ON", "threshold_sf2":"1.956", "check_sf2_igm_brca1":"ON", "threshold_sf2_igm_brca1":"1.867", "check_sc35":"ON", "threshold_sc35": "2.383", "check_srp40":"ON", "threshold_srp40":"2.67", "check_srp55":"ON", "threshold_srp55":"2.676", "seq":">" + genename + "\n" + seq, "submit":"Send"}})

	#parse the data
	soup = BeautifulSoup(response.text, "lxml")
	tables = soup.find_all('table')[3:8] #five tables
	for i, table in enumerate(tables):
		for tr in table.find_all('tr'):
			pos = 0
			val = 0
			for j, td in enumerate(tr.find_all('td')):
				text = re.findall("\>(.*)\<", (td.renderContents()).decode("utf-8"))
				if text != None and len(text) > 0:
					text = max(text, key=len)
					if j == 0:
						pos = int(text.split()[0])
					if j == 2:
						val = float(text.strip())
			if pos != 0 and val != 0 and abs(val) > abs(data[pos]):
				data[pos] = val
	#convert positions to ORF
	newdata = {}
	init_time = time.time()
	for i in data.keys():
		j, error = convert_position(subseq, orig_seq, i)
		if not error:
			newdata[j] = data[i]
	return(newdata)

#runs ESEfinder on overlapping segments of length 4999 (the webserver limit is 5000)
#4999 - maxlen is the length of the overlapping segment (doesn't need to be too long)
def run_ESEfinder_wrapper(genomic_aligned="", ORF_aligned="", genename="1", space="-", maxlen=4900, quiet=False):
	genomic = genomic_aligned.replace(space, "")
	data = {i:0 for i in range(1, 1+len(genomic))}
	init_time = time.time()
	for i in range(math.ceil(len(genomic)/float(maxlen))+1):
		subseq = genomic[max(maxlen*i-2500, 0):min(maxlen*i+2499, len(genomic))]
		pos = genomic.find(subseq)
		aligned_start = convert_position2(genomic_aligned, pos)
		aligned_end = convert_position2(genomic_aligned, pos+4999)
		if len(ORF_aligned[aligned_start:aligned_end].replace(space, "")) > 1:
			subseq_aligned = (pos * "-") + subseq + ("-" * (len(genomic) - len(subseq) - pos))
			aligned_seqs = {"seq":genomic, "subseq":subseq_aligned} #faster and uses less resources to align sequencess
			tmpdata = tools.retry_func(run_ESEfinder, positional_arguments=[aligned_seqs["subseq"], aligned_seqs["seq"]], keyword_arguments={"genename":genename, "space":space, "quiet":quiet})
			for j in tmpdata:
				data[j] = max([tmpdata[j], data[j]], key=abs)
			tools.update_time(i, math.ceil(len(genomic)/maxlen)+1, init_time, func_name="ESEfinder")
	newdata = {}
	for i in range(1, 1+len(ORF_aligned.replace(space,""))):
		j, error = convert_position(ORF_aligned, genomic_aligned, i)
		if j in data and not error:
			newdata[i] = data[j]
		else:
			newdata[i] = 0
	return([newdata[i] for i in range(1, 1+len(ORF_aligned.replace(space, "")))])

#calls NetSurfP2 on a sequence, returns data
def run_netsurfp(ORF, genename, path=ram_disk):
	if tools.infer_alphabet(ORF) == "nt":
		aaseq = translate_seq(ORF)
	else:
		aaseq = ORF
	tools.write_fasta({genename:aaseq}, path+"AA_ORF.fasta")
	cmd = "netsurfp2 --csv " + path + genename + "_netsurf.csv --hhdb " + uniclust + " hhblits " + hhsuite_model + " " + path + "AA_ORF.fasta " + path
	output = subprocess.getoutput(cmd)
	df = pd.read_csv(path + genename + "_netsurf.csv", sep=",", header=0, index_col=2)
	return(list(df["rsa"]), list(df["asa"]), list(df["q3"]), list(df["q8"]), list(df["disorder"]))

def parse_netphos(ORF):
	if tools.infer_alphabet(ORF) == "nt":
		aaseq = translate_seq(ORF)
	else:
		aaseq = ORF
	response = requests.post("http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi", data={"configfile": "/usr/opt/www/pub/CBS/services/NetPhos-3.1/NetPhos.cf", "SEQPASTE": aaseq, "SEQSUB": "(binary)", "res_type": "_", "threshold": "0", "outformat": "classical", "graphics": "ps"})
	jobid = max(re.findall("name=\"jobid\" value=\"(\S+)\"", response.text), key=len)
	timer = 20
	while '<b>Explain</b></a> the output.  Go <a\nhref="javascript:history.back()"><b>back</b></a>.</font>' not in str(response.text) and timer >= 0:
		time.sleep(20)
		timer -= 1
		response = tools.retry_request(requests.get, ["http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi?jobid=" + jobid])

	if timer < 0:
		warnings.warn("\\033[93m" + "NetPhos timed out." + "\\033[0m")
	data = {}
	try:
		lines = str(response.text).split("# -------------------------------------------------------------------")[1].split("\n")
	except:
		lines = []
	for line in lines:
		pieces = line.split()
		if len(pieces) >= 6:
			try:
				data[int(pieces[2])] = float(pieces[5])
			except:
				pass
	return([0 if j+1 not in data else data[j+1] for j in range(len(aaseq))])

#queries and parses NetOGlyc
def parse_netOglyc(ORF):
	if tools.infer_alphabet(ORF) == "nt":
		aaseq = translate_seq(ORF)
	else:
		aaseq = ORF
	response = tools.retry_request(requests.post, ["http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi"], {"data":{"configfile":"/var/www/html/services/NetOGlyc-4.0/webface.cf", "SEQPASTE":aaseq, "SEQSUB":"(binary)"}})
	jobid = max(re.findall("Job status of ([A-Z\d]+)\<", response.text), key=len)
	timer = 20
	while '#seqname	source	feature	start	end	score	strand	frame	comment' not in str(response.text) and timer >= 0:
		time.sleep(20)
		timer -= 1
		response = tools.retry_request(requests.get, ["https://services.healthtech.dtu.dk//cgi-bin/webface2.cgi?jobid=" + jobid])

	if timer < 0:
		warnings.warn("\\033[93m" + "NetOGlyc timed out." + "\\033[0m")
	data = {}
	lines = str(response.text).split("\n")
	for i in range(14, len(lines)):
		line = lines[i]
		pieces = line.split("\t")
		if len(pieces) >= 6:
			try:
				for j in range(int(pieces[3]), int(pieces[4])+1):
					data[j] = float(pieces[5])
			except:
				pass
	return([0 if j+1 not in data else data[j+1] for j in range(len(aaseq))])

def parse_netNglyc(ORF):
	if tools.infer_alphabet(ORF) == "nt":
		aaseq = translate_seq(ORF)
	else:
		aaseq = ORF
	response = tools.retry_request(requests.post, ["http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi"], {"data":{"configfile": "/usr/opt/www/pub/CBS/services/NetNGlyc-1.0/NetNGlyc.cf", "SEQPASTE": aaseq, "SEQSUB": "(binary)", "id": "", "graphics": "on"}})
	jobid = max(re.findall("name=\"jobid\" value=\"(\S+)\"", response.text), key=len)
	timer = 20
	while '<b>Explain</b></a> the output.  Go <a\nhref="javascript:history.back()"><b>back</b></a>.</font>' not in str(response.text) and timer >= 0:
		response = tools.retry_request(requests.get, ["http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi?jobid=" + jobid])
		time.sleep(20)
		timer -= 1

	if timer < 0:
		warnings.warn("\\033[93m" + "NetNGlyc timed out." + "\\033[0m")
	data = {}
	try:
		lines = str(response.text).split("----------------------------------------------------------------------")[2].split("\n")
	except:
		lines = []
	for line in lines:
		pieces = line.split()
		if len(pieces) >= 4:
			data[int(pieces[1])] = float(pieces[3])
	return([0 if j+1 not in data else data[j+1] for j in range(len(aaseq))])

def run_netphos(ORF, path=ram_disk):
	if tools.infer_alphabet(ORF) == "nt":
		aaseq = translate_seq(ORF)
	else:
		aaseq = ORF
	tools.write_fasta({"1":aaseq}, path + "AA_ORF.fasta")
	cmd = "netphos -b " + path + "AA_ORF.fasta"
	output = subprocess.getoutput(cmd)
	data = {}
	for i, line in enumerate(output.split("\n")):
		pieces = line.split()
		if i > 5 and len(pieces) > 5:
			data[int(pieces[2])] = float(pieces[5])
	for i in range(len(aaseq)):
		if i+1 not in data:
			data[i+1] = 0
			
	return([data[i] for i in range(min(data.keys()), max(data.keys()))])

def run_netNglyc(ORF, path=ram_disk):
	if tools.infer_alphabet(ORF) == "nt":
		aaseq = translate_seq(ORF)
	else:
		aaseq = ORF
	tools.write_fasta({"1":aaseq}, path + "AA_ORF.fasta")
	cmd = "netNglyc " + path + "AA_ORF.fasta"
	output = subprocess.getoutput(cmd)
	data = {}
	section = 0
	for i, line in enumerate(output.split("\n")):
		line = line.strip()
		if len(line) > 50 and len(line.replace("-", "")) == 0:
			section += 1
		if section == 2:
			pieces = line.split()
			if len(pieces) > 3:
				data[int(pieces[1])] = float(pieces[3])
	for i in range(len(aaseq)):
		if i+1 not in data:
			data[i+1] = 0
			
	return([data[i] for i in range(min(data.keys()), max(data.keys()))])

def run_netOglyc(ORF, path=ram_disk):
	if tools.infer_alphabet(ORF) == "nt":
		aaseq = translate_seq(ORF)
	else:
		aaseq = ORF
	tools.write_fasta({"1":aaseq}, path + "AA_ORF.fasta")
	cmd = "netOglyc " + path + "AA_ORF.fasta"
	output = subprocess.getoutput(cmd)
	data = {}
	section = 0
	for i, line in enumerate(output.split("\n")):
		line = line.strip()
		if len(line) > 50 and len(line.replace("-", "")) == 0:
			section += 1
		if section == 1:
			pieces = line.split()
			if len(pieces) > 6:
				data[int(pieces[2])] = float(pieces[3])
	for i in range(len(aaseq)):
		if i+1 not in data:
			data[i+1] = 0
	return([data[i] for i in range(min(data.keys()), max(data.keys()))])

def fas_ess(seq):
	start = "<font color=\"red\">"
	end = "</font>"
	response = tools.retry_request(requests.post, ["http://hollywood.mit.edu/cgi-bin/fas-ess.pl"], {"data":{"sequence":seq, "set":"FAS-hex2"}})
	lines = str(response.text).split("\n")
	hits = []
	for i, line in enumerate(lines):
		if i >= 13 and i <= len(lines) - 6:
			pos = 0
			j = 0
			while j < len(line):
				if len(line[j:]) >= 5 and line[j:j+5] == "<font":
					hits.append([pos+1])
					j += len(start) - 2
				elif len(line[j:]) >= 6 and line[j:j+6] == "</font":
					hits[-1].append(pos)
					j += len(end) - 2
				elif line[j] in ["A", "G", "C", "T", " "]:
					pos += 1
				j += 1
	return(sorted(hits, key=lambda k: k[0]))

def get_exonscan(seq):
	response = tools.retry_request(requests.post, ["http://hollywood.mit.edu/cgi-bin/exonscan-results.pl"], {"data":{"sequence":seq}})
	lines = str(response.text).split("\n")
	data = {}
	for i, line in enumerate(lines):
		nums = re.findall("\-?\d+", line)
		if len(nums) >= 7:
			data[(nums[0], nums[1])] = {"3’ss":nums[2], "5’ss":nums[3], "ESE":nums[4], "ESS":nums[5], "GGG":nums[6], "total":nums[7]}
	return(data)
'''
def get_HSF(genomic_aligned, ORF_aligned, genename="", space="-", quiet=False):
	response = requests.post("http://www.umd.be/HSF/4DACTION/input_SSF", data = [("choix_analyse","ssf_sequence"), ("autoselect","yes"), ("snp_select","no"), ("showonly","yes"), ("geneStatus","all"), ("transcriptStatus","all"), ("nuclposition","100"), ("choix_bdd","gene_name"), ("champlibre",genename), ("paramfulltables","onlyvariants"), ("fenetreintron","yes"), ("fenetretaille","24"), ("paramimages","no"), ("matrice_3","no"), ("Matrice","PSS"), ("Matrice","maxent"), ("seuil_maxent5","3"), ("seuil_maxent3","3"), ("seuil_nnsplice5","0.4"), ("seuil_nnsplice3","0.4"), ("Matrice","BPS"), ("Matrice","ESE finder"), ("seuil_sf2","72.98"), ("seuil_sf2_esef","1.867"), ("seuil_sc35","75.05"), ("seuil_sc35_esef","2.383"), ("seuil_srp40","78.08"), ("seufil_srp40_esef","2.67"), ("seufil_srp55","73.86"), ("seufil_srp55_esef","2.676"), ("Matrice","RESCUE ESE"), ("Matrice","ESE New"), ("seuil_9g8","59.245"), ("seuil_tra2","75.964"), ("Matrice","Sironi"), ("seufil_sironi1","60"), ("seuil_sironi2","60"), ("seuil_sironi3","60"), ("Matrice","Decamers"), ("Matrice","ESS hnRNP"), ("seuil_hnrnpa1","65.476"), ("Matrice","PESE"), ("Matrice","ESR"), ("Matrice","EIE")])
	ensembles = re.findall("name=\\\'choix_transcrit\\\' value=\\\'(ENST\d+)\\\'", response.text)
	if len(ensembles) > 0:
		response = requests.post("http://www.umd.be/HSF/4DACTION/input_AnalyseChoixGene", data=[("choix_transcrit", ensembles[0]), ("ESE Finde", "ESE Finder"), ("RESCUE ESE", "RESCUE ESE"), ("PESE", "PESE"), ("Sironi", "Sironi"), ("Decamers", "Decamers"), ("ESS hnRNP", "ESS hnRNP"), ("ESE New", "ESE New"), ("BPS", "BPS"), ("PSS", "PSS"), ("ESR", "ESR"), ("EIE", "EIE"), ("maxent", "maxent"), ("nuclposition", "100"), ("nuclposition5", "0"), ("nuclposition3", "0"), ("intron_startD", "0"), ("intron_endD", "0"), ("seuil_sf2", "72.98"), ("seuil_sf2ig", "70.51"), ("seuil_srp40", "78.08"), ("seuil_sc35", "75.05"), ("seuil_srp55", "73.86"), ("seuil_tra2", "75.964"), ("seuil_9g8", "59.245"), ("seuil_hnrnpa1", "65.476"), ("seuil_maxent5", "3"), ("seuil_maxent3", "3"), ("snp_select", "no"), ("choix_analyse", "ssf_sequence"), ("matrice_3", "no"), ("showonly", "no"), ("exon_number", "0")])
		exons = re.findall("option value=\\\'(\S+)\\\'\>Exon", response.text)
		for exon in exons:
			response = requests.post("http://www.umd.be/HSF/4DACTION/input_AnalyseChoixEnsembl", data=[("nb_exon", exon), ("ESE Finder", "ESE Finder"), ("RESCUE ESE", "RESCUE ESE"), ("PESE", "PESE"), ("Sironi", "Sironi"), ("Decamers", "Decamers"), ("ESS hnRNP", "ESS hnRNP"), ("ESE New", "ESE New"), ("BPS", "BPS"), ("PSS", "PSS"), ("ESR", "ESR"), ("EIE", "EIE"), ("maxent", "maxent"), ("nuclposition", "100"), ("nuclposition5", "0"), ("nuclposition3", "0"), ("intron_startD", "0"), ("intron_endD", "0"), ("seuil_sf2", "72.98"), ("seuil_sf2ig", "70.51"), ("seuil_srp40", "78.08"), ("seuil_sc35", "75.05"), ("seuil_srp55", "73.86"), ("seuil_tra2", "75.964"), ("seuil_9g8", "59.245"), ("seuil_hnrnpa1", "65.476"), ("seuil_maxent5", "3"), ("seuil_maxent3", "3"), ("snp_select", "no"), ("choix_analyse", "ssf_sequence"), ("matrice_3", "no"), ("showonly", "no"), ("exon_number", "0"), ("transcript_id", ensembles[0])])
			soup = BeautifulSoup(response.text, "html.parser")
'''

#computes EVmutation on each domain for either the amino acid or nucleotide sequence
#useful when full sequence is too long to run EVmutation, but may not produce complete data
#domains are provided in AA sequence coordinates
def domain_EV(ORF,  genename, domains = {}, input_seqs={}, mode="nt", space="-", path=ram_disk):
	if domains == None or len(domains) == 0:
		if mode in ["p", "protein", "amino acid", "aa"]:
			domains = {"entire":[1, len(ORF)//3]}
		else:
			domains = {"entire":[1, len(ORF)]}
	if mode in ["p", "protein", "amino acid", "aa"]:
		data = {"EVmutation " + aa: {} for aa in rev_translate.keys() if aa != "*"}
	else:
		data = {"EVmutation(nt) A":{}, "EVmutation(nt) C":{}, "EVmutation(nt) G":{}, "EVmutation(nt) T":{}}
	alphabet = [key.split()[-1].strip() for key in data.keys()]
	for name, domain in sorted(domains.items(), key=lambda kv: kv[1][0]):
		try:
			if mode in ["p", "protein", "amino acid", "aa"]:
				subseq = translate_seq(ORF)[domain[0]:domain[1]]
			else:
				subseq = ORF[3*domain[0]:3*domain[1]]
			available_memory = (psutil.virtual_memory().available/1073741824)
			if (mode in ["p", "protein", "amino acid", "aa"] and available_memory <= (1/486.0)**2*((len(subseq)*3)**2)) or available_memory <= (1/806.0)**2*(len(subseq)**2):
				raise RuntimeError("Empirical memory limit reached for PLMC for domain " + name)
			print("Computing EVmutation for domain " + str(name) + (" "*20), end="\r")
			if genename not in input_seqs:
				seqs = run_blast(subseq, genename + "_" + name, mode=mode, write=False, space=space, quiet=True)
				alignment = align_blast(seqs, genename + "_" + name)
				tools.write_fasta(alignment, path + "tmp_msa.fasta")
			else:
				if mode in ["p", "protein", "amino acid", "aa"]:
					start = convert_position2(input_seqs[genename], domain[0])
					end = convert_position2(input_seqs[genename], domain[1])
				else:
					start = convert_position2(input_seqs[genename], 3*domain[0])
					end = convert_position2(input_seqs[genename], 3*domain[1])
				seqs = {name:seq[start:end] for name, seq in input_seqs.items()}
				tools.write_fasta(seqs, path + "tmp_msa.fasta")
			create_analysis.create_analysis(fasta=path + "tmp_msa.fasta", focus=genename + "_" + name, skip_align_muts=True, model_params="./tmp.model_params", mode="matrix", alphabet=None, space=space, quiet=True)
			ev_data = pd.read_csv(path + genename + "_" + name + "_outfile.csv", sep=",", index_col=0)
			for i, col in enumerate(ev_data):
				if col.strip() in alphabet and i > 0 and mode in ["p", "protein", "amino acid", "aa"]:
					for j in range(domain[0],domain[1]):
						data["EVmutation " + col.strip()][j] = list(ev_data[col])[j-domain[0]]
				elif col.strip() in alphabet and i > 0:
					for j in range(3*domain[0],3*domain[1]):
							data["EVmutation(nt) " + col.strip()][j] = list(ev_data[col])[j-3*domain[0]]
			
		except RuntimeError as e:
			print("\033[93m" + str(e) + "\033[0m")

	if mode in ["p", "protein", "amino acid", "aa"]:
		for i in range(len(ORF)//3):
			for aa in alphabet:
				if i not in data["EVmutation " + aa]:
					data["EVmutation " + aa][i] = ""
	else:
		for i in range(len(ORF)):
			for nt in alphabet:
				if i not in data["EVmutation(nt) " + nt]:
					data["EVmutation(nt) " + nt][i] = ""
	for col in data:
		data[col] = [data[col][i] for i in range(len(data[col]))]
	return(data)

#returns domains from pfam
def get_pfam(acc, ORF):
	aaseq = translate_seq(ORF)
	response = tools.retry_request(requests.get, ["https://pfam.xfam.org/protein/" + acc + "?output=xml"])
	soup = BeautifulSoup(response.text, "html.parser")
	seq = str(soup.find("sequence").text)
	aligned_seqs = tools.align_sequences({"Uniprot":seq, "ORF":aaseq})

	domains = {}
	matches = soup.find_all("match")
	for match in matches:
		name = max(re.findall(" id=\"(\S+)\"", str(match)), key=len)
		domaintype = max(re.findall(" type=\"(\S+)\"", str(match)), key=len)
		start = max(re.findall(" start=\"(\d+)\"", str(match)), key=len)
		end = max(re.findall(" end=\"(\d+)\"", str(match)), key=len)
		pos1, error1 = convert_position(aligned_seqs["Uniprot"], aligned_seqs["ORF"], int(start))
		pos2, error2 = convert_position(aligned_seqs["Uniprot"], aligned_seqs["ORF"], int(end))
		if not error1 and not error2 and "pfam" in domaintype.lower():
			if name in domains:
				j = 2
				while name + "(" + str(j) + ")" in domains:
					j += 1
				name = name + "(" + str(j) + ")"
			domains[name] = [pos1, pos2]
	return(domains)	

#given two aligned sequences and a position in the first sequence (before alignment), computes the corresponding position in the second sequence (before alignment)
#position1 is assumed to be in 1-indexed position (for position1=0, prints a warning)
def convert_position(seq1, seq2, position1, space="-"):
	error = None
	if position1 == 0:
		warnings.warn("\033[93m" + "Position given is not 1-indexed in function " + str(inspect.stack()[1].function) + "\033[0m")
		error = "Position given is not 1-indexed\n" + str(inspect.stack()[1].function)

	i1 = 0
	i2 = 0
	increment = 0
	while i1 < int(position1) and increment < len(seq1) and increment < len(seq2):
		if seq1[increment] != space:
			i1 += 1
		if seq2[increment] != space:
			i2 += 1
		increment += 1
	if not seq1[increment-1] == seq1.replace(space, "")[i1-1]:
		warnings.warn("\033[93m" + "Positions improperly converted (" + str(inspect.stack()[1].function) + ") at position " + str(position1) + ": " + seq1[increment-1] + " " + seq1.replace(space, "")[i1-1] + "\033[0m")
	elif seq2[increment-1] != space and not seq2[increment-1] == seq2.replace(space, "")[i2-1]:
		warnings.warn("\033[93m" + "Positions improperly converted (" + str(inspect.stack()[1].function) + ") at position " + str(position1) + ": " + seq2[increment-1] + " " + seq2.replace(space, "")[i2-1] + "\033[0m")
	if seq2[increment-1] == space:
		error = "Sequence 1 position aligns with a gap in sequence 2\n" + str(inspect.stack()[1].function)
	return(i2, error)

def convert_position2(seq, position, space="-"):
	charcount = 0
	increment = 0
	while charcount <= position and increment < len(seq):
		if seq[increment] != space:
			charcount += 1
		increment += 1
	return(increment-1)

def hex_score(pos, seq, scores, index=0):
	val = 0
	valid = 0
	for i in range(6):
		subseq = seq[pos-index-i:pos-index-i+6]
		try:
			val += scores[subseq]
			valid += 1
		except KeyError as e:
			pass
	if valid > 0:
		return(val/valid)
	else:
		return(None)

#reads a fasta file in seqdir corresponding to genename
#assumes 3 sequences: ORF, transcript, and genomic
#ORF should only include CDS's
#transcript should be ORF + UTRs
#genomic should splice into transcript/ORF
def parse_seqs(seqdir, genename, space="-"):
	seqs = {"ORF": "", "genomic": "", "transcript": "", "ORF_aligned": "", "genomic_aligned": "", "transcript_aligned":""}
	f = os.path.join(seqdir, genename + ".fasta")
	newseqs = tools.read_fasta(f, aformat="WHOLE")
	if set(seqs.keys()) == set(newseqs.keys()) and all([True if "aligned" in k or space not in v else False for k,v in newseqs.items()]): #newseqs is already properly annotated
		return(newseqs)
	else: #assume sequences are aligned
		sorted_seqs = sorted(newseqs.items(), key=lambda kv: len(kv[1].replace(space, "")))
		seqs["ORF_aligned"] = sorted_seqs[0][1]
		seqs["ORF"] = seqs["ORF_aligned"].replace(space, "")
		seqs["transcript_aligned"] = sorted_seqs[1][1]
		seqs["transcript"] = seqs["transcript_aligned"].replace(space, "")
		seqs["genomic_aligned"] = sorted_seqs[2][1]
		seqs["genomic"] = seqs["genomic_aligned"].replace(space, "")
	return(seqs)

def impute_ORFs(seqs, focus, identifiers, space="-", searchrange=100, matches=[["ATG"], ["TAA", "TAG", "TGA"]]):
	if isinstance(identifiers, str): #if you feed an ORF as a start/end sequence identifier, split it into two
		identifiers = [identifiers[:len(identifiers)//2], identifiers[-len(identifiers)//2:]]
	r = {k:[None, None] for k in seqs.keys()}
	r[focus] = [seqs[focus].replace(space, "").index(identifiers[0]), seqs[focus].replace(space, "").index(identifiers[1]) + len(identifiers[1])] #find location of WT ORF in sequence
	print(r[focus])	
	r[focus] = [convert_position2(seqs[focus], r[focus][i], space) for i in range(len(r[focus]))] #convert location to MSA coords
	print(r[focus])
	for k in seqs.keys():
		for j in [0,1]:
			possible = ["", ""]
			for i in range(searchrange):
				if r[focus][j] - i >= 0 and seqs[k][r[focus][j] - i] != space: #search left of region
					if len(possible[0]) < 3:
						possible[0] = seqs[k][r[focus][j] - i] + possible[0]
					else:
						possible[0] = seqs[k][r[focus][j] - i] + possible[0][:2]
				if r[focus][j] + i < len(seqs[k]) and seqs[k][r[focus][j] + i] != space: #search right of region
					if len(possible[1]) < 3:
						possible[1] = possible[1] + seqs[k][r[focus][j] + i]
					else:
						possible[1] = possible[1][-2:] + seqs[k][r[focus][j] + i]
				if j == 1:
					if possible[0] in matches[j]: #start (stop) codon found
						r[k][j] = r[focus][j] - i
					elif possible[1] in matches[j]: #start (stop) codon found
						r[k][j] = r[focus][j] + i
					if r[k][j] != None:
						break
				elif j == 0:
					if possible[0] in matches[j]: #start (stop) codon found
						r[k][j] = r[focus][j] - i
					elif possible[1] in matches[j]: #start (stop) codon found
						r[k][j] = r[focus][j] + i
					if r[k][j] != None:
						tmploc = len([i for i in range(r[k][j]+1) if seqs[k][i] != "-"]) - 2
						r[k][j] = convert_position2(seqs[k], r[k][j], space)
						break
	cut = {}
	for k in seqs.keys():
		if tools.is_numeric(r[k][0]) and tools.is_numeric(r[k][1]):
			cut[k] = seqs[k][r[k][0]:r[k][1]+1]
		elif not tools.is_numeric(r[k][0]) and tools.is_numeric(r[k][1]):
			cut[k] = seqs[k][:r[k][1]+1]
			modval = len(cut[k].replace("-", "")) % 3
			if modval != 0:
				val = len(cut[k].replace("-", "")) - modval
				converted_pos = convert_position2(cut[k][::-1], val, space)
				cut[k] = cut[k][::-1][:converted_pos][::-1]
		elif tools.is_numeric(r[k][0]) and not tools.is_numeric(r[k][1]):
			cut[k] = seqs[k][r[k][0]:]
			modval = len(cut[k].replace("-", "")) % 3
			if modval != 0:
				val = len(cut[k].replace("-", "")) - modval
				converted_pos = convert_position2(cut[k], val, space)
				cut[k] = cut[k][:converted_pos]
	return(r, cut)

def quality_ORFs(seqs, focus, pctlen = 0.8):
	seqs = {k:v for k,v in seqs.items() if len(v.replace("-", "")) >= pctlen * len(seqs[focus].replace("-", "")) and len(v.replace("-", "")) % 3 == 0}
	align_codons.align_codons(seqs, write=False)
	return(tools.read_fasta(os.path.join(ram_disk, "ORFs.fasta")))

def align_against_focus(seqs, focus):
	if focus not in seqs.keys():
		focus = list(seqs.keys())[-1]
	aligned = {}
	init_time = time.time()
	for i, (k,v) in enumerate(seqs.items()):
		newseqs = tools.align_sequences({focus:seqs[focus], k:v})
		start = convert_position2(newseqs[focus], 1)
		end = convert_position2(newseqs[focus], len(seqs[focus]))
		aligned[k] = newseqs[k][start-1:end+1]
		tools.update_time(i, len(seqs), init_time)
	aligned = tools.align_sequences(aligned)
	return(aligned)

def get_memory(state, maxmem):
	minmem = sys.float_info.max #background memory usage
	while state.value == 0:
		mem = psutil.virtual_memory().used/1073741824
		maxmem.value = max(maxmem.value, mem)
		minmem = min(minmem, mem)
		time.sleep(5)
	maxmem.value -= minmem

#reads all fasta files in directory, then finds the ORF and genomic (pre-spliced) sequences and computes gene-specific data
def precompute(genename, seqs, ids={}, taxid=9606, path=ram_disk, outdir="./gene_data/", space="-", quiet=False, outfile=None):
	if outfile == None or len(outfile) < 1:
		outfile = genename
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	init_time = time.time()
	if not quiet:
		print("\033[95m" + "Computing gene features for " + genename + "\033[0m")
	data = {}
	
	try:
		p_consurf = get_consurf_seq(seqs["ORF"], genename, path=path, space=space, wait_time=60, give_up=12, quiet=quiet)	
	except Exception as e:
		p_consurf = None
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	
	try:
		p_aa = Process(target=run_blast, args=(seqs["ORF"], genename, "aa", "nr", "blastp", True, None, path, space))
		p_nt = Process(target=run_blast, args=(seqs["ORF"], genename, "nt", "nt", "blastn", True, None, path, space))
		p_aa.start()
		p_nt.start()
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	
	try:
		if "Uniprot" in ids and not pd.isnull(ids["Uniprot"]):
			uniprot = parse_uniprot(genename, seqs["ORF"], accession=ids["Uniprot"].strip(), taxid=taxid)
		else:
			uniprot = parse_uniprot(genename, seqs["ORF"], taxid=taxid)
		for ss in uniprot:
			if ss != "Accession" and ss != "Domain":
				data[ss] = []
				for i in range(len(seqs["ORF"])//3):
					flag = 0
					for r in uniprot[ss]:
						if ss != "Disulfide Bond" and r[0] <= i+1 and r[1] >= i+1:
							flag = 1
						elif ss == "Disulfide Bond" and (r[0] == i+1 or r[1] == i+1):
							flag = 1
					data[ss].append(flag)
		ids["Uniprot"] = uniprot["Accession"]
	except Exception as e:
		ids["Uniprot"] = None
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	
	try:
		data["Domain"] = []
		domains = copy.deepcopy(uniprot["Domain"])
		for i in range(len(seqs["ORF"])//3):
			flag = ""
			for domain, r in domains.items():
				if re.match("\(\d+\)?", domain): #domain name is a duplicate
					domain = max(re.findall("(.+)\(\d+\)?", domain), key=len)
				if r[0] <= i+1 and r[1] >= i+1:
					flag = domain
			data["Domain"].append(flag)
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	
	if not quiet:
		print("Uniprot accession found: " + uniprot["Accession"])

	pdbdat = parse_pdb(uniprot["Accession"], genename, ORF=seqs["ORF"])
	
	if "ConsurfDB" in ids and not pd.isnull(ids["ConsurfDB"]): #append given ConsurfDB PDB ID to front of list
		pdb = ids["ConsurfDB"].split(":")[0]
		chains = ids["ConsurfDB"].split(":")[1].split(",")
		pdbdat.insert(0, [pdb, chains])
	if "PDBePISA" in ids and not pd.isnull(ids["PDBePISA"]): #append given PDBePISA PDB ID to front of list
		pdb = ids["PDBePISA"].split(":")[0]
		chains = ids["PDBePISA"].split(":")[1].split(",")
		pdbdat.insert(0, [pdb, chains])
	#the webservers may not have data for some PDBs, so we try until we find a winner
	pdbpisa_flag = False
	consurf_flag = False
	for pdb in pdbdat: 
		if not pdbpisa_flag:
			try:
				data["Accessible surface area (PDBPisa)"], data["Buried surface area (PDBPisa)"], data["Solvation energy (PDBPisa)"], data["Bonds"] = get_PDBPisa(pdb[0], pdb[1], seqs["ORF"])
				if not quiet:
					print("PDBPisa data computed using 3D info from " + pdb[0] + " with chain(s) " + ", ".join(pdb[1]))
				ids["PDBePISA"] = pdb[0] + ":" + ",".join(pdb[1])
				pdbpisa_flag = True
			except (ValueError, KeyError):
				ids["PDBePISA"] = None
				pass
		
		if not consurf_flag:
			try:
				data["Conservation (with PDB)"] = get_consurf(pdb[0], pdb[1], seqs["ORF"])
				if not quiet:
					print("Consurf data computed using 3D info from " + pdb[0] + " with chain(s) " + ", ".join(pdb[1]))
				ids["ConsurfDB"] = pdb[0] + ":" + ",".join(pdb[1])
				consurf_flag = True
			except ValueError:
				ids["ConsurfDB"] = None
				pass
		
	if not pdbpisa_flag:
		print("\033[93m" + "Unable to find a PDBPisa entry" + "\033[0m")
	if not consurf_flag:
		print("\033[93m" + "Unable to find a Consurf entry with 3D info" + "\033[0m")

	try:
		data["Phosphorylation potential"] = parse_netphos(seqs["ORF"])
		data["N-linked Glycosylation potential"] = parse_netNglyc(seqs["ORF"])
		data["O-linked Glycosylation potential"] = parse_netOglyc(seqs["ORF"])
	except Exception as e:
		data["Phosphorylation potential"] = run_netphos(seqs["ORF"], path=path)
		data["N-linked Glycosylation potential"] = run_netNglyc(seqs["ORF"], path=path)
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	
	try:
		data["ESEfinder"] = run_ESEfinder_wrapper(seqs["genomic_aligned"], seqs["ORF_aligned"], genename=genename, space=space, quiet=quiet)
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	
	if not quiet:
		print("Computing miRNA binding with miRDB")
	try:
		pholder, data["miRNA binding (miRDB)"] = get_miRDB(genename, seqs["genomic_aligned"], seqs["ORF_aligned"], taxid=taxid)
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	
	if not quiet:
		print("Computing surface area with NetSurfP2")
	try:
		data["Relative surface area (NetsurfP2)"], data["Relative surface area (NetsurfP2)"], data["SS3 (NetsurfP2)"], data["SS8 (NetsurfP2)"], data["Disorder"] = run_netsurfp(seqs["ORF"], genename)
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	'''
	try:
		if not quiet:
			print("Computing relative mRNA MFE values")
		data["mRNA_energies"] = relative_mRNA_MFE(seqs["ORF_aligned"], seqs["transcript_aligned"], space=space)
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	'''
	if not quiet:
		print("Finding and computing MSAs for AA and NT sequences")
	try:
		p_aa.join()
		aa_seqs = tools.read_fasta(path + "aa_msa.fasta", aformat="WHOLE")
		if not genename in aa_seqs or len(aa_seqs) < 5: #if BLAST ran into problems
			aa_seqs = tools.read_fasta(os.path.join(path,"0000_" + genename, "0000_" + genename + "_MSA.a3m"), aformat="WHOLE") #use the NetSurfP2 output fasta
			aa_seqs[genename] = copy.deepcopy(aa_seqs["0000_" + genename])
			del aa_seqs["0000_" + genename]
		aaalignment, data["AA conservation"], data["AA entropy"], data["AA variance"], data["BLOSUM conservation"] = compute_conservation(aa_seqs, genename, mode="p")
		tools.write_fasta(aaalignment, path + "aa_msa.fasta")
		dist = compute_distribution(aaalignment, genename)
		tools.pickledump(dist, outdir + genename + "_aa_dist.pkl")
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")

	try:
		p_nt.join()
		nt_seqs = tools.read_fasta(path + "nt_msa.fasta", aformat="WHOLE")
		ntalignment, data["NT conservation"], data["NT entropy"], data["NT variance"] = compute_conservation(nt_seqs, genename, mode="nt")
		tools.write_fasta(ntalignment, path + "nt_msa.fasta")
		dist = compute_distribution(ntalignment, genename)
		tools.pickledump(dist, outdir + genename + "_nt_dist.pkl")
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	
	model = hmm_seqs(aaalignment, genename, space=space, alphabet="-ACDEFGHIKLMNPQRSTVWY")
	with open(outdir + genename + "_AA_HMM.pkl", "wb") as outf:
		pickle.dump(model, outf)
	model = hmm_seqs(ntalignment, genename, space=space, alphabet="-TAGC")
	with open(outdir + genename + "_NT_HMM.pkl", "wb") as outf:
		pickle.dump(model, outf)
	
	try:
		data["Codon conservation"], data["Codon entropy"], data["Codon variance"] = codon_conservation("nt_msa.fasta", genename, path=ram_disk, space="-")
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")

	try:
		data["Rare codon enrichment"] = rc_enrichment(nt_seqs, genename, quiet=quiet)
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	
	plmc_memory = {"p":0, "nt":0}
	
	state = Value("i", 0)
	memused = Value("d", 0.0)
	p_mem = Process(target=get_memory, args=(state, memused))
	p_mem.start()
	if not quiet:
		print("Computing EVMutation scores based on MSAs under mode p")
	#run EVmutation on the AA sequence
	try:
		available_memory = (psutil.virtual_memory().available/1073741824)
		if available_memory <= (1/486.0)**2*(len(seqs["ORF"])**2):
			raise RuntimeError("Empirical memory limit reached for PLMC.")
		create_analysis.create_analysis(fasta=path + "aa_msa.fasta", focus=genename, skip_align_muts=True, model_params="./aa.model_params", mode="matrix", alphabet=None, space="-")
		ev_data = pd.read_csv(path + genename + "_outfile.csv", sep=",", index_col=0)
		for i, col in enumerate(ev_data):
			if col.strip() in rev_translate.keys():
				data["EVmutation " + col.strip()] = list(ev_data[col])
	except RuntimeError as e:
		print("\033[93m" + str(e) + "\033[0m")
		tmpdata = domain_EV(seqs["ORF"], genename, domains=domains, mode="p", space=space, path=path)
		data.update(tmpdata)
	state.value = 1
	p_mem.join()
	plmc_memory["p"] = memused.value

	state = Value("i", 0)
	memused = Value("d", 0.0)
	p_mem = Process(target=get_memory, args=(state, memused))
	p_mem.start()

	if not quiet:
		print("Computing EVMutation scores based on MSAs under mode nt")
	#run EVmutation on the NT sequence
	try:
		available_memory = (psutil.virtual_memory().available/1073741824)
		if available_memory <= (1/806.0)**2*(len(seqs["ORF"])**2):
			raise RuntimeError("Empirical memory limit reached for PLMC.")
		create_analysis.create_analysis(fasta=path + "nt_msa.fasta", focus=genename, skip_align_muts=True, model_params="./nt.model_params", mode="matrix", alphabet=None, space="-")
		ev_data = pd.read_csv(path + genename + "_outfile.csv", sep=",")
		for i, col in enumerate(ev_data):
			if col.strip() in ["A", "G", "C", "T"] and i > 0:
				data["EVmutation(nt) " + col.strip()] = list(ev_data[col])
	except RuntimeError as e:
		print("\033[93m" + str(e) + "\033[0m")
		tmpdata = domain_EV(seqs["ORF"], genename, domains=domains, mode="nt", space=space, path=path)
		data.update(tmpdata)
	state.value = 1
	p_mem.join()
	plmc_memory["nt"] = memused.value
	print("EVmutation calculations finished")
	
	#read the Consurf output
	if not quiet:
		print("Waiting on Consurf...")
	try:
		if p_consurf != None:
			p_consurf.join()
		df = pd.read_csv(path + genename + "_consurf.tsv", sep="\t", index_col=0)
		data["Consurf (wo 3D info)"] = list(df["score"])
		data["Buried/Exposed"] = list(df["buried/exposed"])
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	
	#add any necessary padding
	for key in data:
		data[key] += ["" for i in range(len(seqs["ORF"])-len(data[key]))]
	#print to file
	try:
		df = pd.DataFrame(data)
		df.index += 1
		df.to_csv(os.path.join(outdir, outfile + ".tsv"), sep="\t")
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")

	elapsed_time = time.time() - init_time
	print(str(len(seqs["ORF"])) + " nts: " + str(round(elapsed_time/60.0, 2)) + " minutes, " + str(round(plmc_memory["p"],5)) + "(aa)/" + str(round(plmc_memory["nt"],5)) + "(nt)GB memory used...\n")
	subprocess.getoutput("beep -l 1000")

	return(data, ids)

def get_codon_counts(seq):
	codondata = {"counts":{codon:0 for codon in translate.keys()},"RSCU":{codon:0 for codon in translate.keys()}}
	codonpairdata = {"counts":{codon1+codon2:0 for codon1 in translate.keys() for codon2 in translate.keys()}, "RSCPU":{codon1+codon2:0 for codon1 in translate.keys() for codon2 in translate.keys()}, "noln CPS":{}, "CPS":{}}
	for i in range(len(seq)//3):
		codondata["counts"][seq[3*i:3*i+3]] += 1
		if 3*i+6 < len(seq):
			codonpairdata["counts"][seq[3*i:3*i+6]] += 1
	for codon1 in codondata["counts"].keys():
		syn_codons = 0
		numsyn = 0
		for codon2 in rev_translate[translate[codon1]]:
			numsyn += codondata["counts"][codon2]
			syn_codons += 1
		try:
			codondata["RSCU"][codon1] = codondata["counts"][codon1]/(numsyn/syn_codons)
		except ZeroDivisionError as e:
			codondata["RSCU"][codon1] = np.nan
	for cp1 in codonpairdata["counts"].keys():
		syn_cp = 0
		numsyn = 0
		for cp2_a in rev_translate[translate[cp1[:3]]]:
			for cp2_b in rev_translate[translate[cp1[3:]]]:
				cp2 = cp2_a+cp2_b
				numsyn += codonpairdata["counts"][cp2]
				syn_cp += 1
		try:
			codonpairdata["RSCPU"][cp1] = codonpairdata["counts"][cp1]/(numsyn/syn_cp)
		except ZeroDivisionError as e:
			codonpairdata["RSCPU"][cp1] = np.nan
	codondata["W"] = {codon1:codondata["RSCU"][codon1]/max([codondata["RSCU"][codon2] for codon2 in rev_translate[translate[codon1]]]) for codon1 in codondata["RSCU"].keys()}
	codonpairdata["W_CP"] = {cp1:codonpairdata["RSCPU"][cp1]/max([codonpairdata["RSCPU"][cp2_a+cp2_b] for cp2_a in rev_translate[translate[cp1[:3]]] for cp2_b in rev_translate[translate[cp1[3:]]]]) for cp1 in codonpairdata["RSCPU"].keys()}
	for cp1 in codonpairdata["counts"].keys():
		f_x = sum([codondata["counts"][codon] for codon in rev_translate[translate[cp1[:3]]]])
		f_y = sum([codondata["counts"][codon] for codon in rev_translate[translate[cp1[3:]]]])
		f_xy = sum([codonpairdata["counts"][cpa+cpb] for cpa in rev_translate[translate[cp1[:3]]] for cpb in rev_translate[translate[cp1[3:]]]])
		try:
			codonpairdata["noln CPS"][cp1] = codonpairdata["counts"][cp1]/(((codondata["counts"][cp1[:3]] * codondata["counts"][cp1[3:]])/(f_x*f_y))*f_xy)
			codonpairdata["CPS"][cp1] = math.log(codonpairdata["noln CPS"][cp1])
		except (ZeroDivisionError, ValueError) as e:
			codonpairdata["noln CPS"][cp1] = np.nan
			codonpairdata["CPS"][cp1] = np.nan
	return(codondata, codonpairdata)

def get_exon_locs(seq, space="-"):
	exons = {}
	count = 0
	for i in range(len(seq)):
		if seq[i] != space:
			count += 1
		if i > 0 and seq[i-1] == space and seq[i] != space:
			if len(exons) < 1:
				exons[1] = [count]
			else:
				exons[max(exons.keys())+1] = [count]
		if i < len(seq)-1 and seq[i+1] == space and seq[i] != space:
			exons[max(exons.keys())].append(count)
	return(exons)

def compute_features(mutfile, codonfile, codonpairfile, seqdir, space="-"):
	aafeats = ["Helix", "Beta strand", "Turn", "Binding site", "Metal binding", "Active site", "Glycosylation (N-linked)", "Glycosylation (O-linked)", "Disulfide bond", "Phosphorylation site", "Cleavage site", "Other modification", "Domain", "Accessible surface area (PDBPisa)", "Buried surface area (PDBPisa)", "Solvation energy (PDBPisa)", "Bonds", "Conservation (with PDB)", "Phosphorylation potential", "N-linked Glycosylation potential", "O-linked Glycosylation potential", "Relative surface area (NetsurfP2)", "SS3 (NetsurfP2)", "SS8 (NetsurfP2)", "Disorder", "AA conservation", "AA entropy", "AA variance", "BLOSUM conservation","Codon conservation", "Codon entropy", "Codon variance", "Rare codon enrichment", "Consurf (wo 3D info)", "Buried/Exposed"]
	ntfeats = ["NT conservation", "NT entropy", "NT variance", "ESEfinder", "miRNA binding (miRDB)"]

	muts = pd.read_csv(mutfile, sep="\t", header=0, index_col=None) #file containing the mutations
	muts = muts[[col for col in muts.columns if "Unnamed" not in str(col) and "column" not in str(col)]].sample(frac=1) #randomly reorder rows and remove "Unnamed" columns
	codon = pd.read_csv(codonfile, sep="\t", header=0, index_col=0) #file containing the genomic codon data
	codonpair = pd.read_csv(codonpairfile, sep="\t", header=0, index_col=0) #file containing the genomic codon pair data
	data = {}
	init_time = time.time()
	last_gene = ""
	codondata = {}
	codonpairdata = {}
	try:
		for j, (k, row) in enumerate(muts.iterrows()):
			if os.path.isfile(seqdir + row[gene] + "_" + row[nm] + ".tsv"):
				row = dict(row)
				if row[gene] + "_" + row[nm] != last_gene:
					seqs = tools.read_fasta(os.path.join(seqdir, row[gene] + "_" + row[nm] + ".fasta"))
					last_gene = row[gene] + "_" + row[nm]
					wt_codondata, wt_codonpairdata = get_codon_counts(seqs["ORF"])
					exon_locs = get_exon_locs(seqs["ORF_aligned"], space=space)

					exons = {}
					u = 0
					for v in range(len(seqs["ORF_aligned"])):
						if seqs["ORF_aligned"][v] != space:
							u += 1
						if v > 0 and seqs["ORF_aligned"][v-1] == space and seqs["ORF_aligned"][v] != space: #new exon starting
							if len(exons) < 1:
								exons[1] = [u]
							else:
								exons[max(exons.keys())+1] = [u]
						if v+1 < len(seqs["ORF_aligned"]) and seqs["ORF_aligned"][v] != space and seqs["ORF_aligned"][v+1] == space: #exon ending
							exons[max(exons.keys())].append(u)
				wtseq = seqs["ORF"]
				mutseq = wtseq[:int(row[pos])-1] + row[nt_mut] + wtseq[int(row[pos]):]
				mut_codondata, mut_codonpairdata = get_codon_counts(mutseq)
				
				#compute codon and codon pair for WT and mutant
				posincodon = (int(row[pos]) -1) % 3
				codonstart = (int(row[pos]) - 1) - posincodon
				c_WT = seqs["ORF"][codonstart:codonstart+3]
				cp1_WT = seqs["ORF"][codonstart-3:codonstart+3]
				cp2_WT = seqs["ORF"][codonstart:codonstart+6]

				c_mut = mutseq[codonstart:codonstart+3]
				cp1_mut = mutseq[codonstart-3:codonstart+3]
				cp2_mut = mutseq[codonstart:codonstart+6]	
		
				#sanity checks
				if len(wtseq) < int(row[pos]):
					warnings.warn("\033[93m" + "Mutation is later than sequence " + row[gene] + " " + str(row[pos]) + ": " + row[nt_wt] + "\033[0m")
					continue
				if wtseq[int(row[pos]) - 1] != row[nt_wt]:
					warnings.warn("\033[93m" + "Nucleotides don't match at " + row[gene] + " " + str(row[pos]) + ": " + row[nt_wt] + " vs " + wtseq[int(row[pos]) - 1] + "\033[0m")
					continue
				'''
				if c_WT != row["Codon (WT)"]:
					warnings.warn("\033[93m" + "Codons don't match at " + row[gene] + " " + str(row[pos]) + ": " + row["Codon (WT)"] + " vs " + c_WT + "\033[0m")
					continue
				'''
				mutid = row[gene] + "_" + row[nt_wt] + str(row[pos]) + row[nt_mut]
				data[mutid] = {key:row[key] for key in row}
				data[mutid]["ORF AA#"] = ((int(row[pos])-1)//3) + 1
				data[mutid]["Position in the codon"] = posincodon + 1
				
				data[mutid].update({"AA (WT)":translate[c_WT], "AA (mut)":translate[c_mut], "Codon (WT)":c_WT, "Codon (mut)":c_mut, "Codon Pair-1 (WT)":cp1_WT, "Codon Pair-2 (WT)":cp2_WT, "Codon Pair-1 (mut)":cp1_mut, "Codon Pair-2 (mut)":cp2_mut})
				if data[mutid]["AA (WT)"] == "P":
					data[mutid]["Proline (in WT)"] = 1
					data[mutid]["Proline Involved"] = 1
				else:
					data[mutid]["Proline (in WT)"] = 0
					if data[mutid]["AA (mut)"] == "P":
						data[mutid]["Proline Involved"] = 1
					else:
						data[mutid]["Proline Involved"] = 0

				if data[mutid]["AA (WT)"] == "C":
					data[mutid]["Cysteine (in WT)"] = 1
					data[mutid]["Cysteine Involved"] = 1
				else:
					data[mutid]["Cysteine (in WT)"] = 0
					if data[mutid]["AA (mut)"] == "C":
						data[mutid]["Cysteine Involved"] = 1
					else:
						data[mutid]["Cysteine Involved"] = 0

				
				for column in ["RSCU", "W"]:
					data[mutid][column + " WT (gene)"] = wt_codondata[column][c_WT]
					data[mutid][column + " mut (gene)"] = mut_codondata[column][c_mut]
					data[mutid]["Δ " + column + " (gene)"] = mut_codondata[column][c_mut] - wt_codondata[column][c_WT]
				for column in ["RSCPU", "W_CP", "noln CPS", "CPS"]:
					data[mutid][column + " WT 1 (gene)"] = wt_codonpairdata[column][cp1_WT]
					data[mutid][column + " mut 1 (gene)"] = mut_codonpairdata[column][cp1_mut]
					data[mutid]["Δ " + column + " 1 (gene)"] = mut_codonpairdata[column][cp1_mut] - wt_codonpairdata[column][cp1_WT]
					data[mutid][column + " WT 2 (gene)"] = wt_codonpairdata[column][cp2_WT]
					data[mutid][column + " mut 2 (gene)"] = mut_codonpairdata[column][cp2_mut]
					data[mutid]["Δ " + column + " 2 (gene)"] = mut_codonpairdata[column][cp2_mut] - wt_codonpairdata[column][cp2_WT]

				exons = {}
				u = 0
				for v in range(len(seqs["ORF_aligned"])):
					if seqs["ORF_aligned"][v] != space:
						u += 1
					if v > 0 and seqs["ORF_aligned"][v-1] == space and seqs["ORF_aligned"][v] != space: #new exon starting
						if len(exons) < 1:
							exons[1] = [u]
						else:
							exons[max(exons.keys())+1] = [u]
					if v+1 < len(seqs["ORF_aligned"]) and seqs["ORF_aligned"][v] != space and seqs["ORF_aligned"][v+1] == space: #exon ending
						exons[max(exons.keys())].append(u)

				for exon, r in exons.items():
					if int(row[pos]) >= r[0] and int(row[pos]) <= r[1]:
						data[mutid]["Exon#"] = exon
						if abs(int(row[pos]) - r[0]) < 7 or abs(int(row[pos]) - r[1]) < 7:
							data[mutid]["Near Splice Junction (< 7 nt)? (1-Yes,  0-No)"] = 1
						else:
							data[mutid]["Near Splice Junction (< 7 nt)? (1-Yes,  0-No)"] = 0
				
				for column in codon:
					data[mutid][column + " WT"] = codon.at[c_WT, column]
					data[mutid][column + " mut"] = codon.at[c_mut, column]
					data[mutid]["Δ " + column] = codon.at[c_mut, column] - codon.at[c_WT, column]
				for column in codonpair:
					try: 
						data[mutid][column + " WT 1"] = codonpair.at[cp1_WT, column]
						data[mutid][column + " mut 1"] = codonpair.at[cp1_mut, column]
						data[mutid]["Δ " + column + " 1"] = codonpair.at[cp1_mut, column] - codonpair.at[cp1_WT, column]
					except:
						data[mutid][column + " WT 1"] = ""
						data[mutid][column + " mut 1"] = ""
						data[mutid]["Δ " + column + " 1"] = ""
					try:
						data[mutid][column + " WT 2"] = codonpair.at[cp2_WT, column]
						data[mutid][column + " mut 2"] = codonpair.at[cp2_mut, column]
						data[mutid]["Δ " + column + " 2"] = codonpair.at[cp2_mut, column] - codonpair.at[cp2_WT, column]
					except:
						data[mutid][column + " WT 2"] = ""
						data[mutid][column + " mut 2"] = ""
						data[mutid]["Δ " + column + " 2"] = ""
				
				pos2, error = convert_position(seqs["ORF_aligned"], seqs["transcript_aligned"], int(row[pos])-1)
				mRNAseq_mut = tools.update_str(seqs["transcript"], row[nt_mut], pos2)
				mRNAseq_wt = seqs["transcript"][max(pos2-75, 0):min(pos2+76, len(seqs["transcript"])+1)]
				mRNAseq_mut = mRNAseq_mut[max(pos2-75, 0):min(pos2+76, len(seqs["transcript"])+1)]
				try:
					data[mutid]["WT ΔG (KineFold)"] = run_kinefold(mRNAseq_wt, seed=0)
					data[mutid]["mut ΔG (KineFold)"] = run_kinefold(mRNAseq_mut, seed=0)
					data[mutid]["ΔΔG (KineFold)"] = data[mutid]["mut ΔG (KineFold)"] - data[mutid]["WT ΔG (KineFold)"]
				except Exception as e:
					data[mutid]["WT ΔG (KineFold)"] = ""
					data[mutid]["mut ΔG (KineFold)"] = ""
					data[mutid]["ΔΔG (KineFold)"] = ""
				
				try:
					data[mutid]["WT ΔG (mFold)"] = run_mfold(mRNAseq_wt)
					data[mutid]["mut ΔG (mFold)"] = run_mfold(mRNAseq_mut)
					data[mutid]["ΔΔG (mFold)"] = data[mutid]["mut ΔG (mFold)"] - data[mutid]["WT ΔG (mFold)"]
				except Exception as e:
					data[mutid]["WT ΔG (mFold)"] = ""
					data[mutid]["mut ΔG (mFold)"] = ""
					data[mutid]["ΔΔG (mFold)"] = ""
				
				try:
					data[mutid]["WT ΔG (remuRNA)"], data[mutid]["mut ΔG (remuRNA)"], data[mutid]["H(wt||mu)"], data[mutid]["GCratio"] = run_remuRNA(mRNAseq_wt, row[nt_wt] + "76" + row[nt_mut])
					data[mutid]["ΔΔG (remuRNA)"] = data[mutid]["mut ΔG (remuRNA)"] - data[mutid]["WT ΔG (remuRNA)"]
				except ValueError as e:
					data[mutid]["WT ΔG (remuRNA)"] = data[mutid]["mut ΔG (remuRNA)"] = data[mutid]["H(wt||mu)"] = data[mutid]["GCratio"] = ""
					data[mutid]["ΔΔG (remuRNA)"] = ""
				
				try:
					data[mutid]["WT ΔG (NUPACK)"] = run_nupack(mRNAseq_wt)
					data[mutid]["mut ΔG (NUPACK)"] = run_nupack(mRNAseq_mut)
					data[mutid]["ΔΔG (NUPACK)"] = data[mutid]["mut ΔG (mFold)"] - data[mutid]["WT ΔG (mFold)"]
				except Exception as e:
					data[mutid]["WT ΔG (NUPACK)"] = ""
					data[mutid]["mut ΔG (NUPACK)"] = ""
					data[mutid]["ΔΔG (NUPACK)"] = ""
				
				genomicpos, error = convert_position(seqs["ORF_aligned"], seqs["genomic_aligned"], row[pos])

				if seqs["genomic"][genomicpos-1] != row[nt_wt]:
					warnings.warn("\033[93m" + "Nucleotides don't match at " + row[gene] + " " + str(genomicpos) + ": " + row[nt_wt] + " vs " + seqs["genomic"][genomicpos-1] + "\033[0m")
			
				for hexfile in ["ESR.tsv", "Z_EI.tsv", "Z_WS.tsv"]:
					hexdat = pd.read_csv(os.path.join(os.path.split(os.path.dirname(seqdir))[0], "sources", hexfile), sep="\t", header=None, index_col=0)
					hexval = 0
					data[mutid][hexfile.split(".")[0] + " WT"] = 0
					for i in range(-5, 1):
						subseq = seqs["genomic"][genomicpos+i-1:genomicpos+i+5]
						try: 
							data[mutid][hexfile.split(".")[0] + " WT"] += float(hexdat.at[subseq, 1])
						except KeyError:
							pass

					data[mutid][hexfile.split(".")[0] + " mut"] = 0
					for i in range(-5, 1):
						subseq = seqs["genomic"][genomicpos+i-1:genomicpos-1] + row[nt_mut] + seqs["genomic"][genomicpos:genomicpos+i+5]
						try: 
							data[mutid][hexfile.split(".")[0] + " mut"] += float(hexdat.at[subseq, 1])
						except KeyError:
							pass

					data[mutid]["Δ " + hexfile.split(".")[0]] = data[mutid][hexfile.split(".")[0] + " mut"] - data[mutid][hexfile.split(".")[0] + " WT"]

				
				try:
					if fas_ess(seqs["genomic"]) != fas_ess(seqs["genomic"][:genomicpos-1] + row[nt_mut] + seqs["genomic"][genomicpos:]):
						data[mutid]["FAS-ESS difference"] = 1
					else:
						data[mutid]["FAS-ESS difference"] = 0
				except Exception as e:
					pass

				try:
					if get_exonscan(seqs["genomic"]) != get_exonscan(seqs["genomic"][:genomicpos-1] + row[nt_mut] + seqs["genomic"][genomicpos:]):
						data[mutid]["Exonscan difference"] = 1
					else:
						data[mutid]["Exonscan difference"] = 0
				except Exception as e:
					pass
				
				for key in aa_data["A"].keys():
					try:
						data[mutid][key + " WT"] = aa_data[translate[c_WT]][key]
						data[mutid][key + " mut"] = aa_data[translate[c_mut]][key]
					except KeyError as e:
						pass

				data[mutid]["Transition (0) vs transversion (1)"] = ""
				if (row[nt_wt] == "A" and row[nt_mut] == "G") or (row[nt_wt] == "G" and row[nt_mut] == "A") or (row[nt_wt] == "C" and row[nt_mut] == "T") or (row[nt_wt] == "T" and row[nt_mut] == "C"):
					data[mutid]["Transition (0) vs transversion (1)"] = 0
				else:
					data[mutid]["Transition (0) vs transversion (1)"] = 1
			
				genedata = pd.read_csv(seqdir + row[gene] + "_" + row[nm] + ".tsv", sep="\t", header=0, index_col=0)
				for col in genedata:
					if col in ntfeats:
						data[mutid][col] = genedata.loc[row[pos], col]
					elif col in aafeats:
						data[mutid][col] = genedata.loc[data[mutid]["ORF AA#"], col]
					elif "EVmutation(nt)" in col and " "+row[nt_mut] in col:
						data[mutid]["EVmutation(nt)"] = genedata.loc[row[pos], col]
					elif "EVmutation" in col and " "+translate[c_mut] in col:
						data[mutid]["EVmutation"] = genedata.loc[data[mutid]["ORF AA#"], col]
				
				#compute VEP things
				d = get_vep(int(row[pos]), (row[nt_wt], row[nt_mut]), row[gene])
				data[mutid].update(d)
				
				tools.update_time(j, len(muts.index), init_time)
				if j % 20 == 1:
					compile_data(data, "./", findstr = "mutdata\_computed")
	except Exception as e:
		not_done = list(muts.index)
		for k,v in data.items():
			genename = "_".join(k.split("_")[:-1])
			var = tools.get_mutation_data(k.split("_")[-1])
			for i, row in muts.iterrows():
				if str(row[gene]) == str(genename) and str(row[pos]) == str(var[0]) and row[nt_wt] == var[1][0] and row[nt_mut] == var[1][1]:
					not_done.remove(i)
		if len(re.findall('\d+$', os.path.splitext(os.path.basename(mutfile))[0])) > 0:
			num = int(max(re.findall('\d+$', os.path.splitext(os.path.basename(mutfile))[0]), key=len))+1
		else:
			num = 2
		newfile = os.path.join(os.path.dirname(mutfile), os.path.splitext(os.path.basename(mutfile))[0].replace(str(num-1), "") + str(num) + ".tsv")
		muts.loc[not_done].to_csv(newfile, sep="\t")
		compile_data(data, "./", findstr = "mutdata\_computed")
		print(e)
		traceback.print_exc()
		data.update(compute_features(newfile, codonfile, codonpairfile, seqdir, space=space))
		return(data)
	pd.DataFrame(data).T.to_csv("./mutdata_computed.tsv", sep="\t")
	return(data)

def compile_data(data, directory, findstr = "mutdata\_computed"):
	now = datetime.datetime.now()
	today = str(now.year) + "_" + str(now.month) + "_" + str(now.day)
	tmpdf = copy.deepcopy(data)
	for f in os.listdir(directory):
		if re.match(findstr, f):
			tmpdf.update(pd.read_csv(f, sep="\t", index_col=0).T.to_dict())
	pd.DataFrame(tmpdf).T.to_csv("mutdata_computed_" + today + ".tsv", sep="\t")

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='This script reads a mutation database and computes multiple features for the mutation.', add_help=True)
	parser.add_argument('mutations', type=str, metavar='PATH', help="path to spreadsheet containing mutations")
	parser.add_argument('codondata', type=str, metavar='PATH', help="path to spreadsheet containing codon data (RSCU, tRNA W, etc.)")
	parser.add_argument('codonpairdata', type=str, metavar='PATH', help="path to spreadsheet containing codon pair data (RSCPU, W, etc.)")
	parser.add_argument('--seqdir', type=str, metavar="DIR", default="./gene_data/", help="path to directory containing fasta files for all sequences")
	parser.add_argument('--tmpdir', type=str, metavar="DIR", default=ram_disk, help="path to directory to store temporary files (preferably in a memory-mounted drive)")
	parser.add_argument('--idfile', type=str, default="", metavar="PATH", help="path to file containing gene name, accession id, and PDB ids") 	
	parser.add_argument('--skip_pre', action="store_true", help="whether or not to compute features for each gene (very time consuming)")
	parser.add_argument('--skip_done', action="store_true", help="whether to skip precomputation for genes that already have spreadsheets")
	parser.add_argument('--quiet', '-v', action="store_true", help="whether to print intermittent updates")
	parser.add_argument('--jobs', type=int, default=1, help="maximum number of jobs to run concurrently")
	parser.add_argument('--taxid', type=int, default=9606, help="taxonomic id of the organism")
	clargs = parser.parse_args()

	genes = []
	for f in os.listdir(clargs.seqdir):
		if f.endswith(".fasta") or f.endswith(".fa") and "aligned" not in f:
			seqs = tools.read_fasta(os.path.join(clargs.seqdir, os.path.splitext(f)[0] + ".fasta"))
			genes.append((os.path.basename(os.path.splitext(f)[0]), len(min(seqs.values(), key=len))))
	genes = sorted(genes, key=lambda kv:kv[1])
	genes = [gene[0] for gene in genes]
	if not clargs.skip_pre:

		if clargs.idfile != None and len(clargs.idfile) > 0: 
			ids = pd.read_csv(clargs.idfile, sep="\t", index_col=0).T.to_dict()
		else:
			ids = {}
		q = Queue()
		p = []

		for genename in genes:
			q.put(genename)

		data = Manager().dict()
		for i in range(clargs.jobs):
			p.append(Process(target=work, args=(q, data, ids, clargs.seqdir, clargs.tmpdir, clargs.seqdir, clargs.skip_done, None, "-", clargs.quiet)))
			p[-1].start()

		for i in range(len(p)):
			p[i].join()

	data = compute_features(clargs.mutations, clargs.codondata, clargs.codonpairdata, clargs.seqdir)
	
