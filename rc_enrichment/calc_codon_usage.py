# Copyright (C) 2017 William M. Jacobs

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

import argparse, sys, os, random, math, gzip, pickle, os.path
from collections import defaultdict
from itertools import combinations, combinations_with_replacement
import numpy as np
import scipy.stats

from codons import codon_to_aa, codon_table

#read a fasta file
def read_fasta(path):
	seqs = []
	keys = []
	with open(path, 'r') as f:
		for line in f:
			if len(line) > 1 and '#' not in line:
				if '>' in line:
					seqs.append('')
					keys.append(line[1:].strip())
				else:
					seqs[-1] += line.strip()
	return {keys[i] : seqs[i] for i in range(len(keys))}

aa_codons = {aa : [c for c in codon_to_aa if codon_to_aa[c] == aa] \
			 for aa in codon_to_aa.values() if aa != 'Stop'}

#determines whether or not a codon is rare
def israre(codon_usage, rare_model, rare_threshold, c):
	if rare_model == 'no_norm':
		if codon_usage[c] <= rare_threshold:
			return True
		else:
			return False
	elif rare_model == 'cmax_norm':
		if codon_usage[c] / max(codon_usage[cc] for cc in aa_codons[codon_to_aa[c]]) <= rare_threshold:
			return True
		else:
			return False

def calc_codon_usage(fasta, abundances=None, output="", rare_model='no_norm', rare_threshold=0.1, max_len_diff=0.2, group_dpercentile=10, wt_gi='gi|556503834|ref|NC_000913.3|', gi_index=None, verbose=False):
	#read fasta files
	seqs = {}
	if isinstance(fasta, str):
		gene = "".join(os.path.basename(fasta).split(".")[:-1])
		seqs[gene] = read_fasta(fasta)
	elif isinstance(fasta, (list, tuple)):
		for path in fasta:
			gene = "".join(os.path.basename(path).split(".")[:-1])
			seqs[gene] = read_fasta(path)
	if verbose:
		print("Loaded sequences for %d genes" % len(seqs))
	gis = sorted(set(gi for gene in seqs for gi in seqs[gene].keys()))

	#read abundance files
	try:
		with open(abundances, 'r') as f:
			abundances = {line.split()[0] : float(line.split()[1]) for line in f if len(line) > 1 and line[0] != '#'}
	except Exception as e:
		abundances = {}
	
	'''	
	if gi_index != None:
		with open(gi_index, 'r') as f:
			gi_index = {line.split()[0] : ' '.join(line.split()[1:]) \
						for line in f if len(line) > 1 and line[0] != '#'}
		print("GIs:")
		for gi in gis:
			print("%32s: %s" % (gi, gi_index[gi]))
	'''

	#delete the sequences whose length differs from the WT too much
	nonwt_gis = [gi for gi in gis if gi != wt_gi]
	for gene in seqs:
		if wt_gi in seqs[gene]:
			wtlen = len(seqs[gene][wt_gi]) - seqs[gene][wt_gi].count('-')
			for gi in nonwt_gis:
				if gi in seqs[gene]:
					gilen = len(seqs[gene][gi]) - seqs[gene][gi].count('-')
					if abs(1. - gilen / wtlen) > max_len_diff:
						del seqs[gene][gi]

	rerun_flag = False
	try: # split sequences into deciles based on rare codon usage (calculated from first run)
		with gzip.open(os.path.join(output, 'input_codon_usage.p.gz'), 'rb') as f:
			input_relative_usage = pickle.load(f)['overall_codon_usage']
		def get_frac_rare(seq):
			return np.mean([sum(1 if israre(input_relative_usage[gi], rare_model, \
											rare_threshold, seq[gi][3*i:3*(i + 1)]) else 0 \
								for i in range(len(seq[gi]) // 3) if seq[gi][3*i:3*(i + 1)] != '---' \
								and codon_to_aa[seq[gi][3*i:3*(i + 1)]] != 'Stop') / \
							sum(1 for i in range(len(seq[gi]) // 3) if seq[gi][3*i:3*(i + 1)] != '---' \
								and codon_to_aa[seq[gi][3*i:3*(i + 1)]] != 'Stop') \
							for gi in gis if gi in seq])
		frac_rare = {gene : get_frac_rare(seq) for gene,seq in seqs.items() if len(seq) > 0}
		groups = ['ND'] + [(np.percentile(list(frac_rare.values()), percentile), \
							np.percentile(list(frac_rare.values()), percentile \
										  + group_dpercentile)) \
						   for percentile in range(0, 100, group_dpercentile)][::-1]
		def get_gene_group(gene):
			if len(seqs[gene]) == 0:
				return 0
			else:
				x = get_frac_rare(seqs[gene])
				for i in range(1, len(groups)):
					if x >= groups[i][0] and x <= groups[i][1]:
						return i
		gene_group_labels = ['%05.3f:%05.3f' % (groups[i][0], groups[i][1]) \
							 if i > 0 else 'ND' for i in range(len(groups))]
	except IOError: #this is the first run, get general usage info
		rerun_flag = True
		groups = ['all']
		def get_gene_group(gene):
			return 0
		gene_group_labels = ['all']
	except KeyError: #code was run in the same output directory, but on a different set of inputs (input_codon_usage.p.gz isn't correct)
		os.remove(os.path.join(output, 'input_codon_usage.p.gz'))
		rerun_flag = True
		groups = ['all']
		def get_gene_group(gene):
			return 0
		gene_group_labels = ['all']
	gene_groups = {gene : get_gene_group(gene) for gene in seqs}
	if verbose:
		print("Gene groups:")
		for i in range(len(gene_group_labels)):
			print("%11s: n = %3d" % (gene_group_labels[i], \
									 sum(1 for gene in seqs if gene_groups[gene] == i)))
	

	#compute codon usage 
	computed_codon_usage = {}
	computed_codon_usage_unw = {}
	computed_codon_usage_groupw = {}
	absolute_usage = {}
	relative_usage = {}
	relative_usage_unw = {}
	relative_usage_groupw = {}
	for gi in gis:
		computed_codon_usage[gi] = defaultdict(int)
		computed_codon_usage_unw[gi] = defaultdict(int)
		computed_codon_usage_groupw[gi] = [defaultdict(int) for i in range(len(groups))]
		for gene,gene_seqs in seqs.items():
			if gi in gene_seqs:
				seq = gene_seqs[gi]
				for i in range(len(seq) // 3):
					c = seq[3*i:3*(i + 1)]
					if c != '---' and codon_to_aa[c] != 'Stop':
						if gene in abundances:
							computed_codon_usage[gi][c] += abundances[gene]
						else:
							computed_codon_usage[gi][c] += 1
						computed_codon_usage_unw[gi][c] += 1
						computed_codon_usage_groupw[gi][gene_groups[gene]][c] += 1
		codons_total_gi = sum(computed_codon_usage[gi].values())
		absolute_usage[gi] = {c : x / codons_total_gi for c,x in computed_codon_usage[gi].items()}
		relative_usage[gi] = {}
		relative_usage_unw[gi] = {}
		relative_usage_groupw[gi] = {i : {} for i in range(len(groups))}
		for aa in aa_codons:
			aa_total_gi = 0
			aa_total_unw_gi = 0
			for c in list(codon_to_aa):
				if codon_to_aa[c] == aa:
					aa_total_gi = aa_total_gi + computed_codon_usage[gi][c]
					aa_total_unw_gi = aa_total_unw_gi + computed_codon_usage_unw[gi][c]
			for c in aa_codons[aa]:
				try:
					relative_usage[gi][c] = computed_codon_usage[gi][c] / aa_total_gi
					relative_usage_unw[gi][c] = computed_codon_usage_unw[gi][c] / aa_total_unw_gi
				except:
					relative_usage[gi][c] = 1.0/len([c in aa_codons[aa]])
					relative_usage_unw[gi][c] = 1.0/len([c in aa_codons[aa]])
			for i in range(len(groups)):
				aa_total_groupw_gi_i = sum(computed_codon_usage_groupw[gi][i][c] for c in aa_codons[aa])
				for c in aa_codons[aa]:
					if aa_total_groupw_gi_i > 0:
						relative_usage_groupw[gi][i][c] \
							= computed_codon_usage_groupw[gi][i][c] / aa_total_groupw_gi_i
					else:
						relative_usage_groupw[gi][i][c] = 0

	if rerun_flag: #first run through, print general codon usage data
		if verbose:
			print("Writing input_codon_usage.p.gz")
		with gzip.open(os.path.join(output, 'input_codon_usage.p.gz'), 'wb') as f:
			pickle.dump({'groups' : groups,
						 'gene_groups' : gene_groups,
						 'overall_codon_usage' : relative_usage,
						 'unweighted_codon_usage' : relative_usage_unw,
						 'gene_group_codon_usage' : relative_usage_groupw}, f)
			if verbose:
				print("WARNING: Rerun analysis to compute frac-rare groups")
	else: #second run through, print group codon usage data
		codon_list = sorted(c for c in codon_to_aa if codon_to_aa[c] != 'Stop')
		rare_codons = {}
		all_rare_codons = defaultdict(int)
		for gi in gis:
			rare_codons[gi] = sorted(c for c in codon_list \
									 if israre(relative_usage[gi], rare_model, \
											   rare_threshold, c))
			for c in rare_codons[gi]:
				all_rare_codons[c] += 1
		if verbose:
			print("Always common codons:", ' '.join(c for c in sorted(codon_list) \
													if c not in all_rare_codons))
			print("Rare codons:")
			for c in sorted(all_rare_codons, key=lambda y: (-all_rare_codons[y], y)):
				print("%s %s %d" % (c, codon_to_aa[c], all_rare_codons[c]))
			print("Writing rare_codons.dat")
		with open(os.path.join(output, 'rare_codons.dat'), 'w') as f:
			for gi in gis:
				f.write("%s %s\n" % (gi, ','.join("%s:%5.3f" % (c, relative_usage_unw[gi][c]) \
												  for c in sorted(rare_codons[gi]))))
		codon_list_aa_sorted = sorted(codon_list, \
									  key=lambda y: (codon_to_aa[y], \
											relative_usage_groupw[wt_gi][len(groups)-1][y]))
		if verbose:
			print("Writing codon_usage.dat")
		with open(os.path.join(output, 'codon_usage.dat'), 'w') as f:
			f.write("# GI gene_group_index gene_group codon_index "
					"amino_acid codon israre relative_usage\n")
			for gi in gis:
				for c in codon_list_aa_sorted:
					if c in rare_codons[gi]:
						israrecodon = 1
					else:
						israrecodon = 0
					for i in range(len(gene_group_labels)):
						f.write("%32s %2d %s %2d %s %s %d %6.4f\n" % \
								(gi, i, gene_group_labels[i], codon_list_aa_sorted.index(c), \
								 codon_to_aa[c], c, israrecodon, relative_usage_groupw[gi][i][c]))
					f.write("\n")
				f.write("\n")
		if verbose:
			print("Writing codon_usage_wt.dat")
		with open(os.path.join(output, 'codon_usage_wt.dat'), 'w') as f:
			f.write("# GI gene_group_index gene_group codon_index "
					"amino_acid codon israre relative_usage\n")
			for c in codon_list_aa_sorted:
				if c in rare_codons[wt_gi]:
					israrecodon = 1
				else:
					israrecodon = 0
				for i in range(len(gene_group_labels)):
					f.write("%32s %2d %s %2d %s %s %d %6.4f\n" % \
							(wt_gi, i, gene_group_labels[i], codon_list_aa_sorted.index(c), \
							 codon_to_aa[c], c, israrecodon, relative_usage_groupw[wt_gi][i][c]))
				f.write("\n")
		if verbose:
			print("Writing codon_usage.p.gz")
		with gzip.open(os.path.join(output, 'codon_usage.p.gz'), 'wb') as f:
			pickle.dump({'groups' : groups,
						 'gene_groups' : gene_groups,
						 'overall_codon_usage' : relative_usage,
						 'unweighted_codon_usage' : relative_usage_unw,
						 'gene_group_codon_usage' : relative_usage_groupw}, f)

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('fasta', type=str, nargs='+', help="path to input MSA fasta file(s)")
	parser.add_argument('--abundances', default=None, type=str, help="path to protein abundance file")
	parser.add_argument('--output', type=str, metavar='PATH', default="", help="path to the directory into which the output should be written")
	parser.add_argument('--rare-model', choices={'no_norm', 'cmax_norm'}, default='no_norm', \
						help="normalization mode for defining rare codons ['no_norm']")
	parser.add_argument('--rare-threshold', type=float, default=0.1, \
						help="threshold for codon rarity [0.1]")
	parser.add_argument('--max-len-diff', type=float, default=0.2, metavar='DIFF', \
						help="maximum relative sequence-length difference compared to the WT [0.2]")
	parser.add_argument('--group-dpercentile', type=int, default=10, metavar='D', \
						help="percentile width for gene-group calculations [10]")
	parser.add_argument('--wt-gi', type=str, default='gi|556503834|ref|NC_000913.3|', \
						help="GI for WT sequence")
	parser.add_argument('--gi-index', type=str, default=None, \
						help="path to index of GIs versus subject titles [None]")
	args = parser.parse_args()

	calc_codon_usage(args.fasta, args.abundances, args.output, args.rare_model, args.rare_threshold, args.max_len_diff, args.group_dpercentile, args.wt_gi, args.gi_index)
