# David Holcomb
# 2018

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

import argparse, prody, gzip, subprocess, os.path
import numpy as np
import re
from itertools import product

from codons import codon_to_aa
from polymer import Polymer
from substructures import find_substructures, calc_co
from calc_consensus_contacts import *
from calc_codon_usage import read_fasta
from clean_duplicates import *
from clean_fasta import *
from tools import write_fasta

#read only the accession key from a string
def get_accession(key):
	return re.findall("\w+\_?\d+\.?\d*", key)[0]

#compare aligned nucleotide sequences to aligned amino acid sequences
def get_quality(aligned_nts, aligned_seq):
	aa_seq = aaseq(aligned_nts)
	return sum(1 if aa_seq[i] == aligned_seq[i] else 0 for i in range(len(aa_seq)))/(len(aa_seq)+0.0)

#produces a dictionary of position:[codon, amino acid] entries for a nucleotide sequence
def aaseq_posdict(nt_seq):
	codons = [nt_seq[3*i:3*(i+1)] for i in range(len(nt_seq)//3)]
	aas = [codon_to_aa[c] for c in codons if codon_to_aa[c] != 'Stop']
	d = {i:[codons[i], aas[i]] for i in range(len(aas))}
	return d
				

def align_codons(fasta, output="aligned_nts.fasta", idthreshold=0.5, lenthreshold=0.2, verbose=False, write=True):
	if isinstance(fasta, str) and os.path.isfile(fasta):
		nt_seqs = read_fasta(fasta)
	else:
		nt_seqs = fasta

	#convert NT sequences to amino acid sequences
	pos_dict = {} #dictionary of dictionaries of pos:[nt, aa] values where key is gi
	aa_seqs = {}  #dictionary of amino acid sequences where key is gi
	for gi in nt_seqs.keys():
		pos_dict[gi] = aaseq_posdict(nt_seqs[gi]) #create a dictionary of pos:[nt, aa] values for gi
		aa_seqs[gi] = ''.join(pos_dict[gi][pos][1] for pos in pos_dict[gi].keys()) #create amino acid sequence for gi
	
	#align the newly created amino acid sequences
	aligned_seqs = align_sequences(nt_seqs)	#align the aa sequences
	'''
	new_aligned_seqs = []
	percid = get_percent_id(aligned_seqs)
	for i, gi in enumerate(list(nt_seqs)):
		if i == 0:
			gi0 = gi
			new_aligned_seqs.append(aligned_seqs[i])
		elif percid[i-1] < idthreshold or len(aa_seqs[gi]) < (1-lenthreshold)*len(aa_seqs[gi0]) or len(aa_seqs[gi]) > (1+lenthreshold)*len(aa_seqs[gi0]):
			del nt_seqs[gi]
			del aa_seqs[gi]
			del pos_dict[gi]
		else:
			new_aligned_seqs.append(aligned_seqs[i])
	aligned_seqs = new_aligned_seqs
	'''
	#convert the aligned amino acid sequences to aligned NT sequences
	aligned_nts = {}
	for i, gi in enumerate(nt_seqs.keys()):
		aligned_nt_seq = ""
		aa_count = 0
		for j in range(len(aligned_seqs[i])):
			if aligned_seqs[i][j] is '-':
				aligned_nt_seq = aligned_nt_seq + '---'
			else:
				if pos_dict[gi][aa_count][1] is not aligned_seqs[i][j]:
					raise Exception("Frame shift occurred somewhere.")
				aligned_nt_seq = aligned_nt_seq + pos_dict[gi][aa_count][0]
				aa_count = aa_count + 1
		aligned_nts[gi] = aligned_nt_seq
	if write:
		with open(output, "w") as f:
			for key in aligned_nts.keys():
				f.write(">" + key + "\n")
				f.write(aligned_nts[key] + "\n")
		f.close()
		'''
		failures = 0
		for i, key in enumerate(nt_seqs.keys()):
			qual = get_quality(aligned_nts[key], aligned_seqs[i])
			if qual < 0.50:
				failures = failures + 1
				print(get_accession(key) + "\t\t" + str(qual))
		print((failures+0.0)/(len(nt_seqs.keys())))
		'''
		#if there are duplicate sequence names, keep only the longest sequence
		clean_duplicates(output, output, aformat="FIRST", verbose=verbose)
	else:
		return(aligned_nts)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='This script aligns multiple nucleotide sequences by converting them to amino acid sequences, aligning those, and converting back to nucleotide sequences.', add_help=True)
	parser.add_argument('fasta', type=str, help="path to fasta sequence")
	parser.add_argument('--output', type=str, metavar='PATH', default='aligned_nts.fasta', \
						help="path to output file")
	parser.add_argument("--idthreshold", type=float, default=0.50, help="minimum percent match to focus sequence")
	parser.add_argument("--lenthreshold", type=float, default=0.20, help="maximum percent length difference relative to focus sequence")
	args = parser.parse_args()

	align_codons(args.fasta, args.output, args.idthreshold, args.lenthreshold)
