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

import argparse, math, random
import networkx as nx
import numpy as np

excluded_regions = {}

def read_rc_profile(path):
	data = {}
	N = 0
	with open(path, 'r') as f:
		for line in f:
			if '#' in line and '=' in line:
				data[line.split()[1]] = int(line.split()[-1])
			elif '#' in line:
				header = line.split()[1:]
				for field in header:
					data[field] = {}
			else:
				N += 1
				center = int(line.split()[0])
				for i in range(1, len(header)):
					data[header[i]][center] = float(line.split()[i])
		data['N'] = N + data['L']
	return data

def join_overlapping_regions(regions):
	G = nx.Graph()
	for i in range(len(regions)):
		G.add_node(regions[i])
		for j in range(i + 1, len(regions)):
			if (regions[j][0] > regions[i][0] and regions[j][0] < regions[i][1]) or \
			   (regions[j][1] > regions[i][0] and regions[j][1] < regions[i][1]):
				G.add_edge(regions[i], regions[j])
	ccs = sorted(nx.connected_components(G))
	return [(min(x for node in cc for x in node), max(x for node in cc for x in node)) for cc in ccs]

def find_conserved_enriched_regions(profile, pval_max, frac_seq_min):
	conserved_enriched = []
	for center in sorted(profile['p_nseq_enriched']):
		if profile['p_nseq_enriched'][center] <= pval_max and \
		   profile['frac_seq_enriched'][center] >= frac_seq_min:
			conserved_enriched.append(center)
	conserved_enriched_indices = sorted(set(i for center in conserved_enriched \
											for i in range(center - profile['L'] // 2, \
														   center + profile['L'] // 2 + 1)))
	ccs = sorted(nx.connected_components(nx.Graph([(i, i + 1) for i in conserved_enriched_indices \
												   if i + 1 in conserved_enriched_indices])))
	regions = [(min(cc), max(cc)) for cc in ccs]
	return regions

def read_elongation_profile(path):
	with open(path, 'r') as f:
		data = [line.split() for line in f if '#' not in line]
	offset = int(data[0][0])
	nvertices = np.array([0]*offset + [int(line[2]) for line in data])
	fe = np.array([0]*offset + [float(line[3]) for line in data])
	return offset, nvertices, fe

def find_fe_drops(nvertices, F, offset=0, L_min_separation=15, min_nvertices=30, min_nvfrac=0.5):
	def isdrop(i, femin):
		if F[i] - femin < 0 \
		   and nvertices[i] >= min_nvertices \
		   and nvertices[i] / i >= min_nvfrac * nvfrac:
			return True
		else:
			return False
	nvfrac = nvertices[-1] / len(nvertices)
	femin = F[0]
	drops = [[]]
	for i in range(len(F)):
		if isdrop(i, femin) and (i == len(F) - 1 or isdrop(i + 1, femin)):
			femin = F[i]
			drops[-1].append((i, femin))
		elif len(drops[-1]) != 0 and i - drops[-1][-1][0] > L_min_separation:
			drops.append([])
	if len(drops[-1]) == 0:
		del drops[-1]
	drops = [((x[0][0], x[-1][0]), x[-1][1] - F[x[0][0] - 1], \
			  max(nvertices[i] / i for i in range(x[0][0], x[-1][0] + 1))) for x in drops]
	return drops

def calc_normalized_conserved_regions(profiles, pval_max, frac_seq_min, binsize=1, centeronly=False):
	normalized_conserved_regions = [0] * (100 // binsize + 1)
	for gene in profiles:
		N = profiles[gene]['N']
		rc = find_conserved_enriched_regions(profiles[gene], pval_max, frac_seq_min)
		if centeronly:
			for i in rc:
				x = round((i[0] + i[1]) / (2 * N) * (len(normalized_conserved_regions) - 1))
				normalized_conserved_regions[x] += 1
		else:
			s = set()
			for i in rc:
				for j in range(i[0], i[1] + 1):
					s.add(round(j / N * (len(normalized_conserved_regions) - 1)))
			for x in s:
				normalized_conserved_regions[x] += 1
	norm = len(profiles)
	return [x / norm for x in normalized_conserved_regions]

def calc_normalized_fe_drops(profiles, binsize=4):
	normalized_drops = [0] * int(math.ceil(100 / binsize) + 1)
	for gene in profiles:
		if 'fe' in profiles[gene]:
			N = len(profiles[gene]['fe'])
		else:
			N = profiles[gene]['N']
		for i in profiles[gene]['fe_drop_boundaries']:
			x = round(i / N * (100 / binsize))
			normalized_drops[x] += 1
	norm = len(profiles)
	return [x / norm for x in normalized_drops]

def calc_normalized_conserved_regions_trunc(profiles, pval_max, frac_seq_min):
	normalized_conserved_regions = [0] * 101
	for gene in profiles:
		N = len(profiles[gene]['p_nseq_enriched'])
		L = profiles[gene]['L']
		rc = find_conserved_enriched_regions(profiles[gene], pval_max, frac_seq_min)
		for i in rc:
			x = round((i[0] + i[1] - 2 * L) / (2 * N) * (len(normalized_conserved_regions) - 1))
			normalized_conserved_regions[x] += 1
	norm = len(profiles)
	return [x / norm for x in normalized_conserved_regions]

def randomize_rc_profiles_from_distribution(data, pval_max, frac_seq_min, verbose=False):
	dist = calc_normalized_conserved_regions_trunc(data, pval_max, frac_seq_min)
	Ls = [r[1] - r[0] - data[gene]['L'] for gene in data \
		  for r in find_conserved_enriched_regions(data[gene], pval_max, frac_seq_min)]
	data_rnd = {}
	for gene in data:
		data_rnd[gene] = data[gene].copy()
		N = len(data[gene]['p_nseq_enriched'])
		L = data[gene]['L']
		Lm = L + np.mean(Ls)
		rc = [(i, max(0.5, 0.5 + np.percentile(Ls, random.random() * 100) / 2)) \
			  for i in data[gene]['p_nseq_enriched'].keys() \
			  if random.random() <= dist[round((i - L) / N * 100)] * (100. / N) * (1. + 2 * Lm / N)]
		if verbose:
			print("randomize:", gene, N, \
				  find_conserved_enriched_regions(data[gene], pval_max, frac_seq_min), '->', rc)
		data_rnd[gene]['p_nseq_enriched'] = {i : 0 for i in data[gene]['p_nseq_enriched'].keys()}
		data_rnd[gene]['frac_seq_enriched'] = {i : 0 for i in data[gene]['frac_seq_enriched'].keys()}
		for r in rc:
			for x in range(int(r[0] - r[1]), int(math.ceil(r[0] + r[1]))):
				data_rnd[gene]['p_nseq_enriched'][x] = pval_max
				data_rnd[gene]['frac_seq_enriched'][x] = frac_seq_min
	return data_rnd

def randomize_rc_profiles_poisson(data, pval_max, frac_seq_min, min_nterm_offset=50):
	n_obs = sum(1 for gene in data \
				for r in find_conserved_enriched_regions(data[gene], pval_max, frac_seq_min) \
				if r[0] >= min_nterm_offset)
	Ltot = sum(max(0, data[gene]['N'] - min_nterm_offset - data[gene]['L']) for gene in data)
	f_avg = n_obs / Ltot
	Ls = [r[1] - r[0] - data[gene]['L'] for gene in data \
		  for r in find_conserved_enriched_regions(data[gene], pval_max, frac_seq_min)]
	data_rnd = {}
	for gene in data:
		data_rnd[gene] = data[gene].copy()
		N = data[gene]['N']
		L = data[gene]['L']
		Lm = L + np.median(Ls)
		def pref(i):
			if i < min_nterm_offset + L:
				return 1. + (i - min_nterm_offset - L // 2) / (L // 2)
			elif i > N - L:
				return 1. + (L // 2 - (N - i)) / (L // 2)
			else:
				return 1.
		rc = [(i, max(0.5, 0.5 + np.percentile(Ls, random.random() * 100) / 2)) \
			  for i in range(min_nterm_offset + L // 2, N - L // 2 + 1) \
			  if random.random() <= f_avg * pref(i) * (1. + 2 * Lm / N)]
		data_rnd[gene]['p_nseq_enriched'] = {i : 0 for i in data[gene]['p_nseq_enriched'].keys()}
		data_rnd[gene]['frac_seq_enriched'] = {i : 0 for i in data[gene]['frac_seq_enriched'].keys()}
		for r in rc:
			for x in range(int(r[0] - r[1]), int(math.ceil(r[0] + r[1]))):
				data_rnd[gene]['p_nseq_enriched'][x] = pval_max
				data_rnd[gene]['frac_seq_enriched'][x] = frac_seq_min
	return data_rnd

def test_rc_explanation(data, pval_max=0.001, frac_seq_min=0.75, \
						max_upstream=20, max_downstream=60, min_nterm_offset=80, verbose=False):
	nfound, ncorrect, ncorrect_end = 0, 0, 0
	for gene in data.keys():
		all_rc_regions = find_conserved_enriched_regions(data[gene], pval_max, frac_seq_min)
		rc_regions = [r for r in all_rc_regions \
					  if not any(rr >= r[0] and rr <= r[1] for rr in excluded_regions.get(gene, []))]
		if verbose and len(all_rc_regions) != len(rc_regions):
			print("excluded RC region:", gene, set(all_rc_regions) - set(rc_regions))
		fe_drop_boundaries = data[gene]['fe_drop_boundaries']
		for r in rc_regions:
			if r[0] >= min_nterm_offset:
				rcenter = (r[0] + r[1]) / 2.
				nfound += 1
				if any(x - y >= max_upstream and x - y <= max_downstream \
					   for y in fe_drop_boundaries for x in r):
					ncorrect += 1
					if data[gene]['N'] - rcenter < max_upstream:
						ncorrect_end += 1
	if nfound == 0:
		frac_explained, frac_explained_err = np.nan, np.nan
	else:
		frac_explained = ncorrect / nfound
		frac_explained_err = math.sqrt(frac_explained * (1. - frac_explained) / nfound)
	return {'number_found' : nfound,
			'frac_explained' : frac_explained,
			'frac_explained_err' : frac_explained_err}

def test_fe_prediction(data, pval_max=1., frac_seq_min=0.75, \
					   max_upstream=20, max_downstream=60, min_nterm_offset=50, verbose=False):
	nboundaries, ncorrect_predict = 0, 0
	nenriched_expected, nregions = 0, 0
	expected_enriched_area, total_area = 0, 0
	for gene in data.keys():
		all_rc_regions = find_conserved_enriched_regions(data[gene], pval_max, frac_seq_min)
		rc_regions = [r for r in all_rc_regions \
					  if not any(rr >= r[0] and rr <= r[1] for rr in excluded_regions.get(gene, []))]
		if verbose and len(all_rc_regions) != len(rc_regions):
			print("excluded RC region:", gene, set(all_rc_regions) - set(rc_regions))
		fe_drop_boundaries = data[gene]['fe_drop_boundaries']
		for y in fe_drop_boundaries:
			nboundaries += 1
			if any(x - y >= max_upstream and x - y <= max_downstream \
				   for r in rc_regions for x in r if r[0] >= min_nterm_offset):
				ncorrect_predict += 1
		expected_enriched_regions = join_overlapping_regions([(y + max_upstream, min(data[gene]['N'], \
										y + max_downstream)) for y in fe_drop_boundaries \
										if y + max_upstream < data[gene]['N']])
		for r in rc_regions:
			if r[0] >= min_nterm_offset:
				rcenter = (r[0] + r[1]) / 2.
				nregions += 1
				if any(rcenter >= reg[0] and rcenter <= reg[1] for reg in expected_enriched_regions):
					nenriched_expected += 1
		expected_enriched_area += sum(reg[1] - reg[0] + 1 for reg in expected_enriched_regions)
		total_area += data[gene]['N'] - min_nterm_offset
	frac_confirmed = ncorrect_predict / nboundaries
	frac_confirmed_err = math.sqrt(frac_confirmed * (1. - frac_confirmed) / nboundaries)
	if nregions == 0:
		frac_predicted, frac_predicted_err = np.nan, np.nan
	else:
		frac_predicted = nenriched_expected / nregions
		frac_predicted_err = math.sqrt(frac_predicted * (1. - frac_predicted) / nregions)
	pred_norm_area = expected_enriched_area / total_area
	return {'number_of_regions' : nregions,
			'frac_confirmed' : frac_confirmed,
			'frac_confirmed_err' : frac_confirmed_err,
			'frac_predicted' : frac_predicted,
			'frac_predicted_err' : frac_predicted_err,
			'frac_predicted_norm_area' : frac_predicted / pred_norm_area,
			'frac_predicted_norm_area_err' : frac_predicted_err / pred_norm_area,
			'frac_notpredicted_norm_area' : (1 - frac_predicted) / (1 - pred_norm_area),
			'frac_notpredicted_norm_area_err' : frac_predicted_err / (1 - pred_norm_area)}

def calc_profile_comparison(genes, structure_path="./", output="./", max_upstream=20, max_downstream=60, Lmax=1e9, blacklist=None):
	excluded_genes = set()
	if blacklist != None:
		with open(lacklist, 'r') as f:
			blacklist = [(line.split()[0], int(line.split()[1]), int(line.split()[-1])) for line in f \
						 if len(line) > 0 and line[0] != '#']
		print("Excluded:", blacklist)
		for d in blacklist:
			if d[-1] == 0:
				if d[0] not in excluded_regions:
					excluded_regions[d[0]] = [d[1]]
				else:
					excluded_regions[d[0]].append(d[1])
			else:
				excluded_genes.add(d[0])

	# Load profile data
	data = {}
	for gene in genes:
		if gene[:2] == 'rp':
			print("WARNING: skipping gene %s ('rp')" % gene)
			continue
		if gene in excluded_genes:
			print("WARNING: skipping excluded gene %s (blacklist)" % gene)
			continue
		try:
			rc_data = read_rc_profile(structure_path + '%s_rc_profile.dat' % gene)
		except IOError:
			rc_data = None
			print("WARNING: could not find rare-codon profile data for gene %s" % gene)
			continue
		try:
			offset, nvertices, fe = read_elongation_profile(structure_path + '%s_elongation_profile.dat' % gene)
			if len(fe) > rc_data['N']:
				print("WARNING: skipping gene %s with inconsistent profile lengths" % gene)
				continue
		except IOError:
			continue
		if rc_data['N'] > Lmax:
			print("WARNING: skipping gene %s with L > Lmax = %d" % (gene, Lmax))
			continue
		data[gene] = rc_data
		data[gene]['f_enriched_avg'] = rc_data['f_enriched_avg'][10]
		data[gene]['nvertices'] = nvertices
		data[gene]['fe'] = fe

	comparison_genes = set(gene for gene in data if data[gene]['n_msa'] >= 8)
	for gene in comparison_genes:
		data[gene]['fe_drop_boundaries'] \
			= [drop[0][0] for drop in find_fe_drops(data[gene]['nvertices'], data[gene]['fe'])]
	comparison_data = {gene : data[gene] for gene in comparison_genes}
	datasets = {'model' : comparison_data}

	# Perform comparisons
	def run_comparison(test_fcn, data, pval_max, frac_seq_min, \
					   randomize_fcns={'c1' : randomize_rc_profiles_from_distribution,
									   'c2' : randomize_rc_profiles_poisson}, ntrials=100):
		test_results = test_fcn(data, pval_max=pval_max, frac_seq_min=frac_seq_min)
		control_results = {k : {key : [] for key in test_results.keys()} for k in randomize_fcns.keys()}
		for k in randomize_fcns:
			randomize_fcn = randomize_fcns[k]
			for trial in range(ntrials):
				randomized_data = randomize_fcn(data, pval_max, frac_seq_min)
				trial_results = test_fcn(randomized_data, pval_max=pval_max, frac_seq_min=frac_seq_min)
				for key in test_results.keys():
					control_results[k][key].append(trial_results[key])
		return test_results, control_results

	for dataset_name, data in datasets.items():
		print("Writing %sgenes_%s.lst" % (output, dataset_name))
		with open('%sgenes_%s.lst' % (output, dataset_name), 'w') as f:
			for gene in sorted(data):
				f.write("%s\n" % gene)

		output_keys = ['frac_explained']
		print("Writing %scompare_profiles_%s_explain.dat" % (output, dataset_name))
		with open('%scompare_profiles_%s_explain.dat' % \
				  (output, dataset_name), 'w') as f:
			f.write("# frac_seq_min pval_max n")
			
			for key in output_keys:
				f.write(" %s %s %s" % (key, key + '_err', \
									   ' '.join("%s %s" % (key + '_control_%s_mean', \
														   key + '_control_%s_std'))))
			f.write("\n")
			
			for frac_seq_min in [0.5, 0.625, 0.75, 0.875, 1.]:
				for pval_max in [0.1**p for p in range(7)]:
					try:
						test_results, control_results = run_comparison(test_rc_explanation, \
																	   data, pval_max, frac_seq_min)
						f.write("%5.3f %8.6f %4d" % \
								(frac_seq_min, pval_max, test_results['number_found']))
						for key in output_keys:
							f.write(" %8.6f %8.6f %s" % \
									(test_results[key], test_results[key + '_err'], \
									 ' '.join("%8.6f %8.6f" % \
											  (np.mean(control_results[k][key]), \
											   np.std(control_results[k][key])) \
											  for k in sorted(control_results.keys()))))
						f.write("\n")
					except ZeroDivisionError:
						continue
				f.write("\n")

		output_keys = ['frac_predicted', 'frac_confirmed', 'frac_predicted_norm_area', \
					   'frac_notpredicted_norm_area']
		print("Writing %scompare_profiles_%s_predict.dat" % (output, dataset_name))
		with open('%scompare_profiles_%s_predict.dat' % \
				  (output, dataset_name), 'w') as f:
			f.write("# frac_seq_min pval_max number_of_regions")
			
			for key in output_keys:
				f.write(" %s %s %s" % (key + '', key + '_err', \
									   ' '.join("%s %s" % (key + '_control_%s_mean', \
														   key + '_control_%s_std') )))
			f.write("\n")
			
			for frac_seq_min in [0.5, 0.625, 0.75, 0.875, 1.]:
				for pval_max in [0.1**p for p in range(7)]:
					try:
						test_results, control_results = run_comparison(test_fe_prediction, \
																	   data, pval_max, frac_seq_min)
						f.write("%5.3f %8.6f %4d" %
								(frac_seq_min, pval_max, \
								 test_results['number_of_regions']))
						for key in output_keys:
							f.write(" %8.6f %8.6f %s" % \
									(test_results[key], test_results[key + '_err'], \
									 ' '.join("%8.6f %8.6f" % \
											  (np.mean(control_results[k][key]), \
											   np.std(control_results[k][key])) \
											  for k in sorted(control_results.keys()))))
						f.write("\n")
					except ZeroDivisionError:
						continue
				f.write("\n")

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('genes', type=str, nargs='+', help="space-separated list of gene name(s)")
	parser.add_argument('--structure-path', type=str, metavar='PATH', default='./', \
						help="path to input files [./]")
	parser.add_argument('--output', type=str, default='./', help="prefix for output files [./]")
	parser.add_argument('--max-upstream', type=int, default=20, \
						help="maximum upstream offset for rare-codon/FE-drop association [20]")
	parser.add_argument('--max-downstream', type=int, default=60, \
						help="maximum downstream offset for rare-codon/FE-drop association [60]")
	parser.add_argument('--Lmax', type=int, default=1e9, help="maximum sequence length allowed [1e9]")
	parser.add_argument('--blacklist', type=str, default=None, help="path to excluded regions [None]")

	args = parser.parse_args()

	calc_profile_comparison(args.genes, args.structure_path, args.output, args.max_upstream, args.max_downstream, args.Lmax, args.blacklist)
