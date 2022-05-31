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

import math, argparse, time
import numpy as np
import networkx as nx
from itertools import combinations
from tools import update_time
import os.path

from polymer import Polymer
from substructures import find_substructures, calc_ghost_edges, \
	calc_independent_configs, get_polymer_segments, get_polymer_from_config, get_expected_segments, \
	get_config_vertices
from fe import avg_contacts, subpolymer_epsilon, minimize_fe, segment_end_probability, find_loops


#unused
def calc_all_configs(polymer, substructures, epsilon, mu, ndim, kuhnlength, verbose=False):
	configs = calc_independent_configs(polymer, substructures)
	polymer_segments = get_polymer_segments(polymer, substructures)
	distances = polymer.distances()
	results = {}
	for L in range(polymer.number_of_residues()):
		allowed_substructures = elongation_substructures(substructures, L)
		config = tuple(sorted(c for c in allowed_substructures))
		if set(config) in configs and config not in results:
			subpolymer = get_polymer_from_config(polymer, polymer_segments, \
												 substructures, config, kuhnlength)
			subpolymer_segments = get_expected_segments(subpolymer, polymer_segments, \
														substructures, config)
			epsilon_config = subpolymer_epsilon(epsilon, subpolymer.nodes())
			rho_config, fe_config, warning = minimize_fe(subpolymer, epsilon_config, distances, \
														 subpolymer_segments, \
														 mu=mu, d=ndim, kuhnlength=kuhnlength)
			avg_seg_end_prob = segment_end_probability(epsilon_config, rho_config, subpolymer_segments)
			results[config] = {'vertices' : list(subpolymer.nodes()),
							   'rho' : rho_config,
							   'fe' : fe_config,
							   'avg_seg_end_prob' : avg_seg_end_prob,
							   'warning' : warning}
	return results

#returns indices of substructures containing only vertices in the first L residues
def elongation_substructures(substructures, L):
	allowed_substructures = [i for i in range(len(substructures))
							 if all(v <= L for b in substructures[i] for v in b)]
	return allowed_substructures

#returns indices of substructures containing any of the first L residues
def get_allowed_partial_substructures(substructures, L):
	allowed_substructures = [i for i in range(len(substructures))
							 if any(all(v <= L for v in b) for b in substructures[i])]
	return allowed_substructures

#returns induced subpolymer containing given vertices
def get_subpolymer(polymer, vertices):
	return Polymer([edge for edge in polymer.edges() if all(u in vertices for u in edge)], \
				   nresidues=polymer.number_of_residues())

#calculates energy for each nascent length
#unused
def calc_fe_L(polymer, distances, substructures, epsilon, mu, ndim, kuhnlength):
	polymer_segments = get_polymer_segments(polymer, substructures)
	def get_allowed_segments(vertices):
		return {c : [[u for u in seg if u in vertices] for seg in polymer_segments[c]] \
				for c in polymer_segments if all(any(u in vertices for u in seg) \
												 for seg in polymer_segments[c])}
	results = {}
	for L in range(min(polymer.nodes()), polymer.number_of_residues() + 1):
		allowed_substructures = get_allowed_partial_substructures(substructures, L)
		allowed_vertices = [u for u in polymer.nodes() if u <= L and any(v <= L for v in polymer[u])]
		allowed_segments = get_allowed_segments(allowed_vertices)
		polymer_L = get_subpolymer(polymer, allowed_vertices)
		epsilon_L = subpolymer_epsilon(epsilon, allowed_vertices)
		rho_L, fe_L, warning = minimize_fe(polymer_L, epsilon_L, distances, allowed_segments, \
										   mu=mu, d=ndim, kuhnlength=kuhnlength)
		avg_seg_end_prob = segment_end_probability(epsilon_L, rho_L, allowed_segments)
		results[L] = {'vertices' : allowed_vertices,
					  'substructures' : allowed_substructures,
					  'rho' : rho_L,
					  'fe' : fe_L,
					  'loops' : find_loops(polymer_L, allowed_segments),
					  'avg_seg_end_prob' : avg_seg_end_prob,
					  'warning' : warning}
	return results

#calculates energy for each nascent length and segment
def calc_fe_Ls(polymer, distances, substructures, epsilon, mu, ndim, kuhnlength):
	polymer_segments = get_polymer_segments(polymer, substructures)
	def get_allowed_segments(vertices):
		return {c : [[u for u in seg if u in vertices] for seg in polymer_segments[c]] \
				for c in polymer_segments if all(any(u in vertices for u in seg) \
												 for seg in polymer_segments[c])}
	results = {}
	r = polymer.number_of_residues() + 1 - min(polymer.nodes())

	#this is the time bottleneck
	###################################################################################
	init_time = time.time()
	for i, L in enumerate(range(min(polymer.nodes()), polymer.number_of_residues() + 1)):
		update_time(i, polymer.number_of_residues()+1 - min(polymer.nodes()), init_time)
		allowed_substructures = get_allowed_partial_substructures(substructures, L)
		ustart = {}
		for Ls in [0] + sorted(min(seg[0] for seg in polymer_segments[c]) \
							   for c in allowed_substructures):
			allowed_vertices = tuple(sorted(u for u in polymer.nodes() if (u <= L and u >= Ls) \
											and any(v <= L and v >= Ls for v in polymer[u])))
			if allowed_vertices not in ustart:
				ustart[allowed_vertices] = Ls
		for allowed_vertices,Ls in ustart.items():
			allowed_segments = get_allowed_segments(allowed_vertices)
			polymer_Ls = get_subpolymer(polymer, allowed_vertices)
			epsilon_Ls = subpolymer_epsilon(epsilon, allowed_vertices)
			rho_Ls, fe_Ls, warning = minimize_fe(polymer_Ls, epsilon_Ls, distances, allowed_segments, \
												 mu=mu, d=ndim, kuhnlength=kuhnlength)
			avg_seg_end_prob = segment_end_probability(epsilon_Ls, rho_Ls, allowed_segments)
			results[(Ls, L)] = {'vertices' : allowed_vertices,
								'substructures' : allowed_substructures,
								'rho' : rho_Ls,
								'fe' : fe_Ls,
								'loops' : find_loops(polymer_Ls, allowed_segments),
								'avg_seg_end_prob' : avg_seg_end_prob,
								'warning' : warning}
	##########################################################################################
	return results

#find the minimum energy set of substructures for each nascent chain length
def min_fe_L_Ls(results_Ls, multiregion=True):
	N = max(k[1] for k in results_Ls)
	results = {}
	stable_regions = []
	def update_stable_regions(Ls, L): #delete all stable regions with Ls <= region <= L and append (Ls, L) to stable regions
		todelete = []
		for i in range(len(stable_regions)):
			if Ls <= stable_regions[i][0] and L >= stable_regions[i][1]:
				todelete.append(i)
		for i in todelete[::-1]:
			del stable_regions[i]
		stable_regions.append((Ls, L))
		stable_regions.sort()
	def min_fe_stable_regions():
		G = nx.Graph()
		for i in range(len(stable_regions)):
			G.add_node(stable_regions[i])
			for j in range(len(stable_regions)):
				if stable_regions[j][0] >= stable_regions[i][0] \
				   and stable_regions[j][0] <= stable_regions[i][1]:
					G.add_edge(stable_regions[i], stable_regions[j])
		ccs = list(nx.connected_components(G))
		selected = []
		for cc in ccs:
			if len(cc) == 1:
				for y in cc:
					selected.append(y)
			else:
				ccmin = (None, 0)
				for n in range(1, len(cc)):
					for x in combinations(cc, n):
						if not any(x[j][0] >= x[i][0] and x[j][0] <= x[i][1] \
								   for i in range(len(x)) for j in range(len(x)) if i != j):
							fe = sum(results_Ls[x[i]]['fe'] for i in range(len(x)))
							if fe < ccmin[1]:
								ccmin = (x, fe)
				for y in ccmin[0]:
					selected.append(y)
		return selected
	def combine_region_results(regions):
		combined_result = {}
		combined_result['vertices'] = sorted(set(u for r in regions for u in results_Ls[r]['vertices']))
		combined_result['substructures'] = sorted(set(c for r in regions \
													  for c in results_Ls[r]['substructures']))
		combined_result['rho'] = np.sum(results_Ls[r]['rho'] for r in regions)
		combined_result['fe'] = sum(results_Ls[r]['fe'] for r in regions)
		combined_result['loops'] = sorted(l for r in regions for l in results_Ls[r]['loops'])
		combined_result['avg_seg_end_prob'] = {k : v for r in regions \
											   for k,v in results_Ls[r]['avg_seg_end_prob'].items()}
		combined_result['warning'] = any(results_Ls[r]['warning'] for r in regions)
		combined_result['Ls'] = max(r[0] for r in regions)
		combined_result['regions'] = sorted(regions)
		return combined_result
	for L in range(N + 1):
		print("# Percent finished: " + str(100.0*L/(N+1)), end='\r', flush=True)
		keys = set(k for k in results_Ls if k[1] == L)
		if len(keys) > 0:
			minLs = min((k for k in keys), key=lambda y: results_Ls[y]['fe']) #find the minimum energy set of substructures
			if multiregion:
				if results_Ls[minLs]['fe'] < 0:
					update_stable_regions(*minLs)
				regions = min_fe_stable_regions()
				if len(regions) > 0:
					for i in range(len(regions), 0, -1):
						if minLs[0] <= regions[i - 1][1]:
							del regions[i - 1]
				regions.append(minLs)
				valid_regions = [r for r in regions if results_Ls[r]['rho'].sum() >= 1. or r[1] == L]
				results[L] = combine_region_results(valid_regions)
			else:
				results[L] = results_Ls[minLs]
				results[L]['Ls'] = minLs[0]
	return results

def calc_elongation(gene, polymer, output = "./", ndim=3.0, mu=-2.0, kuhnlength=2.0):
	with open(polymer, 'r') as f:
		nresidues = int(f.readline().split()[-1])
		energies = {(int(line.split()[0]), int(line.split()[1])) : -float(line.split()[2]) \
					for line in f}
		for edge in list(energies.keys()):
			energies[(edge[1], edge[0])] = energies[edge]
	polymer = Polymer([edge for edge in energies.keys()], nresidues=nresidues)
	distances = polymer.distances()

	print("# nresidues =", polymer.number_of_residues())

	for edge in calc_ghost_edges(polymer):
		polymer.add_edge(*edge)
	epsilon = np.zeros((polymer.number_of_residues(), polymer.number_of_residues()))
	for u,v in energies.keys():
		epsilon[u,v] = epsilon[v,u] = energies[(u,v)]

	native_L = polymer.number_of_residues()
	substructures = find_substructures(polymer)
	results_Ls = calc_fe_Ls(polymer, distances, substructures, epsilon, \
							mu, ndim, kuhnlength)
	results = min_fe_L_Ls(results_Ls)

	print("# Writing elongation_profile.dat")
	with open(os.path.join(output, gene + '_elongation_profile.dat'), 'w') as f:
		f.write("# L Ls nvertices F_L F_L_constrained rho_sum nloops regions\n")
		minfe = 0
		for L in results:
			Ls = results[L].get('Ls', 0)
			minfe = min(minfe, results[L]['fe'])
			regions = results[L].get('regions', [(Ls, L)])
			f.write("%3d %3d %3d %12.6g %12.6g %12.6g %2d %s\n" % \
					(L, Ls, len(results[L]['vertices']), minfe, results[L]['fe'], \
					 results[L]['rho'].sum(), len(results[L]['loops']), str(regions)))

if __name__ == '__main__':

	np.seterr(divide='ignore')

	parser = argparse.ArgumentParser()
	parser.add_argument('gene', type=str, help="gene name")
	parser.add_argument('polymer', type=str, help="path to input absolute-energy polymer file")
	parser.add_argument('--output', type=str, metavar='PATH', default = "./", help="path to the directory to write output to")
	parser.add_argument('--ndim', type=float, default=3., help="number of spatial dimensions [3]")
	parser.add_argument('--mu', type=float, default=-2., help="residue ordering penalty / kT [-2.]")
	parser.add_argument('--kuhnlength', type=float, default=2., metavar='B', help="Kuhn length [2]")
	args = parser.parse_args()

	calc_elongation(args.gene, args.polymer, args.output, args.ndim, args.mu, args.kuhnlength)
