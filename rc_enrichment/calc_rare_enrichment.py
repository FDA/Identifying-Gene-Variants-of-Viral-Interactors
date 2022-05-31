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

import scipy.misc as spmisc
import argparse, sys, math, gzip, pickle, os.path
import numpy as np
import scipy.optimize
from itertools import combinations
from functools import reduce
import operator
import psutil

from codons import codon_to_aa

sense_codons = set(c for c in codon_to_aa if (codon_to_aa[c] != 'Stop' and codon_to_aa[c] != '-'))
aa_codons = {aa : [c for c in codon_to_aa if codon_to_aa[c] == aa] \
			 for aa in codon_to_aa.values() if (aa != 'Stop' and aa != '-')}

def read_fasta(path):
	seqs = []
	keys = []
	with open(path, 'r') as f:
		for line in f:
			if len(line) > 1 and line[0] != '#':
				if '>' in line:
					seqs.append('')
					keys.append(line[1:].strip())
				else:
					seqs[-1] += line.strip()
	return {keys[i] : seqs[i] for i in range(len(keys))}

#replaces all codons after a stop codon with '---'
def clean_sequences(seqs):
	for gi in seqs:
		for i in range(len(seqs[gi]) // 3):
			c = seqs[gi][i*3:(i+1)*3]
			if c != '---' and codon_to_aa[c] == 'Stop':
				#print("# Stop codon encountered at : " + gi + " " + str(i) + c)
				seqs[gi] = seqs[gi][:i*3] + '-' * (len(seqs[gi]) - i * 3)
				break
	return seqs

#get all sequence ids (sorted)
def sorted_gis(seqs, wtgi):
	return [wtgi] + sorted(gi for gi in seqs.keys() if gi != wtgi)

#get codons from a sequence
def get_codons(seq, exclude_stop=False):
	return [seq[3*i:3*(i+1)] for i in range(len(seq)//3) \
			if not exclude_stop or codon_to_aa[seq[3*i:3*(i+1)]] != 'Stop']

def align_sequences(msa_codons, wtgi):
	wt_msa_codons = msa_codons[wtgi]
	wt_ncodons = sum(1 for i in range(len(wt_msa_codons)) if wt_msa_codons[i] != '---')
	msa_indices = []
	j = -1
	for i in range(len(wt_msa_codons)):
		if '-' in wt_msa_codons[i] and wt_msa_codons[i] != '---':
			raise Exception("Frame shift in MSA WT sequence:", wt_msa_codons[i])
		elif j >= wt_ncodons - 1:
			j = wt_ncodons - 1
		elif '-' not in wt_msa_codons[i]:
			j += 1
		if j == -1:
			msa_indices.append(0)
		else:
			msa_indices.append(j)
	return {i : [j for j in range(len(msa_indices)) if msa_indices[j] == i \
				 and any(c[j] != '---' for c in msa_codons.values())] \
			for i in range(wt_ncodons)}

#calculates the percentage of each sequence that matches the WT
def aa_identity(msa_codons, wtgi):
	gis = sorted_gis(msa_codons, wtgi)
	msa_wt_aa_prob = {}
	for gi in gis:
		nalike = sum(1 if msa_codons[gi][i] != '---' and msa_codons[wtgi][i] != '---' \
					 and codon_to_aa[msa_codons[gi][i]] == codon_to_aa[msa_codons[wtgi][i]] \
					 else 0 for i in range(len(msa_codons[wtgi])))
		ntotal = sum(1 if msa_codons[gi][i] != '---' and msa_codons[wtgi][i] != '---' \
					 else 0 for i in range(len(msa_codons[wtgi])))
		if ntotal > 0:
			msa_wt_aa_prob[gi] = nalike / ntotal
		else:
			msa_wt_aa_prob[gi] = np.nan
	return msa_wt_aa_prob

#computes probability of using a rare codon for each amino acid in each sequence (S4)
def gene_avg_codon_probabilities(rare_codons, grp, group_codon_usage, msa_codons, \
								 nstart=30, verbose=False):
	gene_codon_usage = {}
	rare_codon_prob = {}
	def solve_gi(gi): #finds the lambda scaling factor (S5, S6)
		ncodons_gi = sum(1 for c in msa_codons[gi][nstart:] if c != '---' and codon_to_aa[c] != 'Stop')
		overall_frac = sum(1 if c in rare_codons[gi] else 0 for c in msa_codons[gi][nstart:] \
						   if c != '---' and codon_to_aa[c] != 'Stop') / ncodons_gi
		if overall_frac == 0:
			overall_frac = 0.5 / ncodons_gi
		p0 = {aa : sum(group_codon_usage[gi][grp][c] for c in aa_codons[aa] if c in rare_codons[gi]) \
			  for aa in aa_codons}
		x0 = {aa : p0[aa] / (1. - p0[aa]) for aa in p0}
		def diff(lmbda, sln=False):
			p = {aa : (lmbda * x0[aa]) / (lmbda * x0[aa] + 1) for aa in aa_codons}
			expected_frac = sum(p[codon_to_aa[c]] for c in msa_codons[gi][nstart:] \
								if c != '---' and codon_to_aa[c] != 'Stop') / ncodons_gi
			if not sln:
				return expected_frac - overall_frac
			else:
				return p
		def codon_usage(lmbda):
			p = diff(lmbda, sln=True)
			return {aa : {c : group_codon_usage[gi][grp][c] * (p[aa] / p0[aa]) \
						  if c in rare_codons[gi] else \
						  group_codon_usage[gi][grp][c] * ((1. - p[aa]) / (1. - p0[aa])) \
						  for c in aa_codons[aa]} for aa in aa_codons}
		x = scipy.optimize.brentq(diff, 0, 10)
		if verbose:
			print("lambda =", x)
		return codon_usage(x), diff(x, sln=True)
	for gi in msa_codons.keys():
		gene_codon_usage[gi], rare_codon_prob[gi] = solve_gi(gi)
		if any(pc < 0 or pc >= 1 for pc in rare_codon_prob[gi].values()):
			raise Exception("rare codon prob =", rare_codon_prob[gi])
		if any(math.fabs(1. - sum(gene_codon_usage[gi][aa].values())) > 1.e-9 for aa in aa_codons):
			raise Exception("sum(codon prob) != 1", gene_codon_usage[gi])
		if verbose:
			for aa in sorted(aa_codons):
				print("%s:" % aa, ' '.join("%s:%5.3f,%5.3f" % (c, group_codon_usage[gi][grp][c], \
															   gene_codon_usage[gi][aa][c]) \
										   for c in sorted(aa_codons[aa])))
	return gene_codon_usage, rare_codon_prob

#computes the probability of at least n Bernoulli trials being "yes" while the rest are "no"
#p is a list of probability values, corresponding to each Bernoulli trial
def prob_ntuple(p, n):
	entries = set(i for i in range(len(p)))
	q = 1. - p
	if n < len(p) / 2:
		pn = 0.
		for m in range(n):
			for selected in combinations(entries, m):
				notselected = set(entries) - set(selected)
				pn += reduce(operator.mul, [p[j] for j in selected], 1) * \
					  reduce(operator.mul, [q[j] for j in notselected], 1)
		return 1. - pn
	else:
		pn = 0.
		for m in range(n, len(p) + 1):
			for selected in combinations(entries, m):
				notselected = set(entries) - set(selected)
				pn += reduce(operator.mul, [p[j] for j in selected], 1) * \
					  reduce(operator.mul, [q[j] for j in notselected], 1)
		return pn

#same as prob_ntuple, but uses convolutions to avoid redundant calculations
#much faster than prob_ntuple, but introduces (small) rounding errors
def prob_ntuple_conv(p, n, ld=False, verbose=False):
	p_len = len(p)
	pmf = np.zeros(p_len + 1, dtype=np.longdouble)
	if ld:
		p = np.longdouble(p)
	pmf[0] = 1.0
	for i in range(p_len):
		success = pmf[:i + 1] * p[i]
		pmf[:i + 2] *= (1 - p[i])  # failure
		pmf[1:i + 2] += success
	pn = 0.0
	for j in range(n, len(pmf)):
		pn += pmf[j]
	return pn

def msa_rare_codon_analysis_wtalign_nseq(msa_codons, wtgi, msa_index_dict, \
										 rare_codons, rare_codon_prob, L=10, zsig=1., verbose=False):
	gis = sorted_gis(msa_codons, wtgi)
	wt_ncodons = sum(1 for i in range(len(msa_codons[wtgi])) \
					 if msa_codons[wtgi][i] != '---' and codon_to_aa[msa_codons[wtgi][i]] != 'Stop')
	f_enriched_avg = {gi : (sum(1 for i in range(wt_ncodons) for j in msa_index_dict[i] \
							  if msa_codons[gi][j] != '---' \
							  and msa_codons[gi][j] in rare_codons[gi]) / \
						  sum(1 for i in range(wt_ncodons) for j in msa_index_dict[i] \
							  if msa_codons[gi][j] != '---')) for gi in gis}
	f_gi_avg = np.mean(list(f_enriched_avg.values())) #S2
	def prob_poisson(l, n):
		return l**n * math.exp(-l) / math.factorial(n)
	def min_n_poisson_cum(l, z): #returns the enrichment threshold
		pz = math.erf(z)
		s = 0
		n = -1
		while s < pz:
			n += 1
			s += prob_poisson(l, n)
		return n
	def max_n_poisson_cum(l, z): #returns the depletion threshold
		pz = math.erfc(z)
		s = 0
		n = -1
		while s < pz:
			n += 1
			s += prob_poisson(l, n)
		return n
	p_enriched = {gi : {} for gi in gis} #neutral model probability of rare-codon enrichment for gi at center (S7)
	p_depleted = {gi : {} for gi in gis} #neutral model probability of rare-codon depletion for gi at center 
	f_enriched = {gi : {} for gi in gis} #fraction of rare codons for each gi and center
	n_rare = {gi : {} for gi in gis} #number of rare codons
	nseq_enriched = {} #number of enriched sequences at center
	nseq_depleted = {} #number of depleted sequences at center
	fseq_enriched = {} #fraction of sequences enriched with rare codons
	fseq_depleted = {} #fraction of sequences depleted of rare codons
	p_nseq_enriched = {} #neutral model probability of rare-codon enrichment (S8)
	p_nseq_depleted = {} #neutral model probability of rare-codon depletion
	for i in range(wt_ncodons - L + 1):
		if verbose:
			print("# Percent finished: " + str(100.0*i/(wt_ncodons-L+1)), end="\r", flush=True)
		center = i + L // 2
		nseq_enriched[center] = nseq_depleted[center] = 0
		nseq_possible_enriched_center = nseq_possible_depleted_center = 0
		for gi in gis:
			all_indices = sorted(k for j in range(i, i + L) for k in msa_index_dict[j] \
								 if msa_codons[gi][k] != '---') #find all non-insertion/deletion positions of sequence within window
			indices = [j for j in all_indices \
					   if rare_codon_prob[gi][codon_to_aa[msa_codons[gi][j]]] > 0] #find all positions with a nontrivial probability of using a rare codon
			n_rare[gi][center] = sum(1 for j in indices if msa_codons[gi][j] in rare_codons[gi]) #find the number of rare codons of sequence within window
			if len(all_indices) > 0: #avoid a divide-by-zero error
				f_enriched[gi][center] = n_rare[gi][center] / len(all_indices)
			else:
				f_enriched[gi][center] = np.nan
			p_rc = np.array([rare_codon_prob[gi][codon_to_aa[msa_codons[gi][j]]] for j in indices]) #probability of using a rare codon for each gi and amino acid
			nmin_enriched = min_n_poisson_cum(f_gi_avg * len(all_indices), zsig) #rare codon enrichment threshold (S3)
			p_enriched[gi][center] = prob_ntuple(p_rc, nmin_enriched) #use original prob_ntuple to avoid rounding errors
			if n_rare[gi][center] >= nmin_enriched:
				nseq_enriched[center] += 1
			if len(p_rc) >= nmin_enriched:
				nseq_possible_enriched_center += 1
			if len(all_indices) > 0:
				nseq_possible_depleted_center += 1
			nmax_depleted = max_n_poisson_cum(f_gi_avg * len(all_indices), zsig)
			p_depleted[gi][center] = prob_ntuple(1. - p_rc, len(p_rc) - nmax_depleted) #use original prob_ntuple to avoid rounding errors
			if n_rare[gi][center] <= nmax_depleted:
				nseq_depleted[center] += 1
		p_enriched_center = np.array([p_enriched[gi][center] for gi in gis])
		p_nseq_enriched[center] = prob_ntuple_conv(p_enriched_center, nseq_enriched[center]) #use prob_ntuple_conv to speed up
		p_depleted_center = np.array([p_depleted[gi][center] for gi in gis])
		p_nseq_depleted[center] = prob_ntuple_conv(p_depleted_center, nseq_depleted[center]) #use prob_ntuple_conv to speed up
		if nseq_possible_enriched_center > 0:
			fseq_enriched[center] = nseq_enriched[center] / nseq_possible_enriched_center
		else:
			fseq_enriched[center] = np.nan
		if nseq_possible_depleted_center > 0:
			fseq_depleted[center] = nseq_depleted[center] / nseq_possible_depleted_center
		else:
			fseq_depleted[center] = np.nan
		if verbose:
			print("> %4d %2d %2d %6.3f %6.3f" % \
				  (center, nseq_enriched[center], nseq_depleted[center], \
				   np.mean(list(f_enriched[gi][center] for gi in gis)), f_gi_avg))
			print(' p =', ' '.join("%5.3f" % x for x in p_enriched_center))
			print('\a')
	return {'nmin_enriched' : min_n_poisson_cum(f_gi_avg * L, zsig),
			'nmax_depleted' : max_n_poisson_cum(f_gi_avg * L, zsig),
			'p_enriched' : p_enriched,
			'f_enriched' : f_enriched,
			'f_enriched_avg' : f_enriched_avg,
			'n_rare' : n_rare,
			'nseq_enriched' : nseq_enriched,
			'fseq_enriched' : fseq_enriched,
			'p_nseq_enriched' : p_nseq_enriched,
			'nseq_depleted' : nseq_depleted,
			'fseq_depleted' : fseq_depleted,
			'p_nseq_depleted' : p_nseq_depleted}

#loads the null model
def load_null_model(msa_codons, gis, cl_usage, cl_rare_model, cl_use_wt_rare_codons, \
					cl_rare_threshold, cl_null_model, cl_gene, verbose=False):
	with gzip.open(cl_usage, 'rb') as f:
		usage_data = pickle.load(f)
	#relative synonymous codon usage, without n scaling factor, RSCU/(n_i) from Sharp and Li
	if cl_rare_model == 'no_norm': # number of c / number of a <= threshold
		if not cl_use_wt_rare_codons:
			rare_codons = {gi : [c for c in sense_codons \
								 if usage_data['overall_codon_usage'][gi][c] <= cl_rare_threshold] \
						   for gi in gis}
		else:
			rare_codons = {gi : [c for c in sense_codons \
								 if usage_data['overall_codon_usage'][cl_wt_gi][c] \
								 <= cl_rare_threshold] for gi in gis}
	#relative adaptiveness of a codon, w from Sharp and Li
	elif cl_rare_model == 'cmax_norm': # number of c / number of c_max <= threshold
		if not cl_use_wt_rare_codons:
			rare_codons = {gi : [c for c in sense_codons \
								 if usage_data['overall_codon_usage'][gi][c] \
								 / max(usage_data['overall_codon_usage'][gi][c] \
									   for c in aa_codons[codon_to_aa[c]]) <= cl_rare_threshold] \
						   for gi in gis}
		else:
			rare_codons = {gi : [c for c in sense_codons \
								 if usage_data['overall_codon_usage'][cl_wt_gi][c] \
								 / max(usage_data['overall_codon_usage'][cl_wt_gi][c] \
									   for c in aa_codons[codon_to_aa[c]]) <= cl_rare_threshold] \
						   for gi in gis}
	else:
		raise Exception("Unknown rare-codon model")
	if cl_null_model == 'eq': #all codons are equally likely for any amino acid
		gene_codon_usage = {gi : {aa : {c : 1. / len(aa_codons[aa]) for c in aa_codons[aa]} \
								  for aa in aa_codons} for gi in gis}
		rare_codon_prob = {gi : {aa : sum(1. / len(aa_codons[aa]) \
										  for c in aa_codons[aa] if c in rare_codons[gi]) \
								 for aa in aa_codons} for gi in gis}
		relative_usage_giavg = {c : 1. / len(aa_codons[codon_to_aa[c]]) for c in codon_to_aa \
								if codon_to_aa[c] != 'Stop'}
	elif cl_null_model == 'genome': #not weighted by abundance
		gene_codon_usage = {gi : {aa : {c : usage_data['unweighted_codon_usage'][gi][c] \
										for c in aa_codons[aa]} \
								  for aa in aa_codons} for gi in gis}
		rare_codon_prob = {gi : {aa : sum(usage_data['unweighted_codon_usage'][gi][c] \
										  for c in aa_codons[aa] if c in rare_codons[gi]) \
								 for aa in aa_codons} for gi in gis}
		relative_usage_giavg = {c : sum(usage_data['unweighted_codon_usage'][gi][c] for gi in gis) \
								/ len(gis) for c in codon_to_aa if codon_to_aa[c] != 'Stop'}
	elif cl_null_model == 'groups': #uses protein abundance data
		grp = usage_data['gene_groups'][cl_gene]
		gene_codon_usage, rare_codon_prob = gene_avg_codon_probabilities(rare_codons, grp, \
												usage_data['gene_group_codon_usage'], msa_codons, \
																		 verbose=verbose)
		relative_usage_giavg = {c : sum(usage_data['gene_group_codon_usage'][gi][grp][c] \
										for gi in gis) / len(gis) \
								for c in codon_to_aa if codon_to_aa[c] != 'Stop'}
	else:
		raise Exception("Unknown null model: %s" % cl_null_model)
	return rare_codons, gene_codon_usage, rare_codon_prob, relative_usage_giavg

def calc_rare_enrichment(gene, msa, usage, output="./", rare_model="no_norm", rare_threshold=0.1, max_len_diff=0.2, min_aa_iden=0.5, null_model="genome", L=15, wt_gi='gi|556503834|ref|NC_000913.3|', use_wt_rare_codons=False, verbose=False, wt_only=False):
	# Load multiple sequence alignment

	insufficient_overlap = []
	insufficient_aa_ident = []

	seqs = read_fasta(msa)
	seqs = clean_sequences(seqs)
	all_gis = sorted(gi for gi in seqs.keys())
	if len(all_gis) == 0:
		raise Exception("No sequences loaded; gene:", gene)
	if wt_gi not in all_gis:
		raise Exception("WT GI not in GIs; gene:", gene)
	wt_len = len(seqs[wt_gi]) - seqs[wt_gi].count('-')
	for gi in all_gis:
		gi_len = len(seqs[gi]) - seqs[gi].count('-')
		wt_gi_overlap = sum(1 for i in range(len(seqs[wt_gi])) \
							if seqs[wt_gi][i] != '-' and seqs[gi][i] != '-')
		if abs(1. - gi_len / wt_len) > max_len_diff or \
		   1. - wt_gi_overlap / wt_len > max_len_diff:
			insufficient_overlap.append([gi, str(round(wt_gi_overlap / wt_len,5))])
			del seqs[gi]
	if wt_only:
		seqs = {wt_gi : seqs[wt_gi]}

	if len(insufficient_overlap) > 0 and verbose:
		print("Ignoring GIs due to insufficient overlap:")
		for gi in insufficient_overlap:
			print(gi[0] + "\t" + gi[1])

	#compare sequences to WT, remove those that don't match WT well enough
	gis = sorted_gis(seqs, wt_gi)
	msa_codons = {gi : get_codons(seq) for gi,seq in seqs.items()}
	aa_perc_id = aa_identity(msa_codons, wt_gi)
	for gi in gis:
		if aa_perc_id[gi] < min_aa_iden:
			insufficient_aa_ident.append([gi, str(round(aa_perc_id[gi],5))])
			del seqs[gi]
	if len(seqs) == 1:	
		sys.stderr.write("WARNING: Only one usable sequence in alignment; gene: %s\n" % gene)

	if len(insufficient_aa_ident) > 0 and verbose:
		print("Ignoring GIs due to insufficient amino acid identity:")
		for gi in insufficient_aa_ident:
			print(gi[0] + "\t" + gi[1])

	gis = sorted_gis(seqs, wt_gi)
	msa_codons = {gi : get_codons(seq) for gi,seq in seqs.items()}
	msa_index_dict = align_sequences(msa_codons, wt_gi)
	if verbose:
		print("# Loaded sequence alignment with %d sequences" % len(seqs))
		print("# Maximum allowed length difference (relative to WT) =", max_len_diff)
		print("# Minimum allowed AA percent identity (relative to WT) =", min_aa_iden)

	# Load codon-usage data and define null model
	rare_codons, gene_codon_usage, rare_codon_prob, relative_usage_giavg \
		= load_null_model(msa_codons, gis, usage, rare_model, use_wt_rare_codons, \
						  rare_threshold, null_model, gene, verbose=verbose)
	if verbose:
		print("WT rare-codon probabilities:")
		for aa in sorted(aa_codons):
			if rare_codon_prob[wt_gi][aa] > 0:
				print("%s %5.3f" % (aa, rare_codon_prob[wt_gi][aa]))

	# Rare-codon calculations
	if verbose:
		print("# Window width =", L)
	rc_analysis = msa_rare_codon_analysis_wtalign_nseq(msa_codons, wt_gi, \
													   msa_index_dict, rare_codons, rare_codon_prob, \
													   L=L, verbose=verbose)

	with open(output + "%s_final_MSA.fasta" % gene, "w") as f:
		for gi in seqs:
			f.write(">%s\n" % gi)
			f.write(seqs[gi] + "\n")

	with open(output + '%s_rc_profile.dat' % gene, 'w') as f:
		f.write("# n_msa = %d\n" % len(msa_codons))
		f.write("# L = %d\n" % L)
		f.write("# i p_nseq_enriched p_nseq_depleted f_enriched_wt f_enriched_avg "
				"f_enriched_mean f_enriched_stddev frac_seq_enriched frac_seq_depleted "
				"nmin_enriched nmin_depleted p_enriched\n")
		f_enriched_avg = np.array([rc_analysis['f_enriched_avg'][gi] for gi in gis])
		for i in sorted(rc_analysis['p_nseq_enriched'].keys()):
			f_enriched = np.array([rc_analysis['f_enriched'][gi][i] for gi in gis])
			f.write("%d %g %g %g %g %g %g %g %g %g %g %g\n" % \
					(i, rc_analysis['p_nseq_enriched'][i], rc_analysis['p_nseq_depleted'][i], \
					 rc_analysis['f_enriched'][wt_gi][i], np.nanmean(f_enriched_avg), \
					 np.nanmean(f_enriched), np.nanstd(f_enriched), \
					 rc_analysis['fseq_enriched'][i], rc_analysis['fseq_depleted'][i], \
					 rc_analysis['nmin_enriched'], rc_analysis['nmax_depleted'], \
					 rc_analysis['p_enriched'][wt_gi][i]))

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('gene', type=str, help="gene name")
	parser.add_argument('msa', type=str, help="path to MSA fasta file")
	parser.add_argument('usage', type=str, help="path to gzip'd codon usage data")
	parser.add_argument('--output', type=str, default='./', help="prefix for output files ['./']")
	parser.add_argument('--rare-model', choices={'no_norm', 'cmax_norm'}, default='no_norm', \
						help="normalization mode for defining rare codons ['no_norm']")
	parser.add_argument('--rare-threshold', type=float, default=0.1, \
						help="threshold for codon rarity [0.1]")
	parser.add_argument('--max-len-diff', type=float, default=0.2, \
						help="maximum relative sequence-length difference compared to the WT [0.2]")
	parser.add_argument('--min-aa-iden', type=float, default=0.5, \
						help="minimum amino-acid percent identity compared to the WT [0.5]")
	parser.add_argument('--null-model', choices={'groups', 'genome', 'eq'}, \
						default='genome', help="codon-usage null model ['groups']")
	parser.add_argument('--L', type=int, default=15, \
						help="window width for local rare-codon concentration calculation, "
						"in codons [15]")
	parser.add_argument('--wt-gi', type=str, default='gi|556503834|ref|NC_000913.3|', \
						help="GI for WT sequence [=Escherichia coli str. K-12 substr. MG1655]")
	parser.add_argument('--use-wt-rare-codons', action='store_true', \
						help="use WT rare codons for all GIs [False]")
	parser.add_argument('--verbose', action='store_true', help="print more information [False]")
	parser.add_argument('--wt-only', action='store_true', help="examine WT only instead of MSA [False]")
	args = parser.parse_args()

	calc_rare_enrichment(args.gene, args.msa, args.usage, args.output, args.rare_model, args.rare_threshold, args.max_len_diff, args.min_aa_iden, args.null_model, args.L, args.wt_gi, args.use_wt_rare_codons, args.verbose, args.wt_only)
