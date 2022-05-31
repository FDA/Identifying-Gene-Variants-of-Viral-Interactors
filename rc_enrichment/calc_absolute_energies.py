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

import argparse, os.path, time
import numpy as np
import networkx as nx
import scipy.optimize
from itertools import product
from tools import update_time

from polymer import Polymer
from substructures import find_substructures, get_polymer_segments, calc_ghost_edges
from fe import minimize_fe, subpolymer_epsilon
from calc_consensus_contacts import write_polymer

#creates a list of (list of ranges connected by edges)
def native_connected_components(polymer, ranges):
	G = nx.Graph()
	for ri in ranges:
		for rj in ranges:
			if ri not in G or (ri in G and rj not in G[ri]):
				for u,v in polymer.edges():
					if u >= ri[0] and u < ri[1] and v >= rj[0] and v < rj[1]: #there is an edge connecting both ranges
						G.add_edge(ri, rj)
						break
	#G is a graph where vertices represent ranges, and there is an edge between ranges if there is an edge in the polymer connecting them
	return sorted(sorted(cc) for cc in nx.connected_components(G))

#creates a subpolymer containing only amino acids between umin and umax
def get_subpolymer(polymer, ranges):
	umin = min(u for r in ranges for u in r)
	umax = max(u for r in ranges for u in r)
	subpolymer = Polymer([edge for edge in polymer.edges() \
						  if all(u >= umin and u < umax for u in edge)], \
						 nresidues=polymer.number_of_residues())
	return subpolymer

#returns indices of substructures entirely between Lmin and Lmax
def get_allowed_partial_substructures(substructures, Lmin, Lmax):
	allowed_substructures = [i for i in range(len(substructures))
							 if any(all(v >= Lmin and v <= Lmax for v in b) for b in substructures[i])]
	return allowed_substructures

#returns the minimum free energy over all start and end vertices
def minimize_fe_Lse(polymer, epsilon, distances, segments, mu, d, kuhnlength, min_len_frac=0.8):
	Lmin = min(u for u in range(polymer.number_of_residues()) \
			   if any(epsilon[u,v] != 0 for v in range(epsilon.shape[1]))) #first amino acid involved in a native contact
	Lmax = max(u for u in range(polymer.number_of_residues()) \
			   if any(epsilon[u,v] != 0 for v in range(epsilon.shape[1]))) + 1 #last amino acid involved in a native contact
	substructures = find_substructures(polymer)
	allowed_substructures = get_allowed_partial_substructures(substructures, Lmin, Lmax)
	polymer_segments = get_polymer_segments(polymer, substructures)
	def get_allowed_segments(vertices): 
		return {c : [[u for u in seg if u in vertices] for seg in polymer_segments[c]] \
				for c in polymer_segments if all(any(u in vertices for u in seg) \
												 for seg in polymer_segments[c])}
	def get_allowed_subpolymer(vertices):
		return Polymer([edge for edge in polymer.edges() if all(u in vertices for u in edge)], \
				   nresidues=polymer.number_of_residues())
	ustart = set([Lmin] + [min(seg[0] for seg in polymer_segments[c]) for c in allowed_substructures]) #set of earliest vertices for all segments in all substructures
	uend = set([max(seg[-1] for seg in polymer_segments[c]) for c in allowed_substructures] + [Lmax]) #set of latest vertices for all segments in all substructures
	minfe = (None, None, np.inf, None)
	for Ls,Le in product(sorted(ustart), sorted(uend)): #consider every possible subset
		if Le - Ls >= min_len_frac * (Lmax - Lmin):
			allowed_vertices = tuple(sorted(u for u in polymer.nodes() if (u >= Ls and u <= Le) \
											and any(v >= Ls and v <= Le for v in polymer[u])))
			allowed_segments = get_allowed_segments(allowed_vertices)
			polymer_Lse = get_allowed_subpolymer(allowed_vertices)
			epsilon_Lse = subpolymer_epsilon(epsilon, allowed_vertices)
			rho_Lse, fe_Lse, warning = minimize_fe(polymer_Lse, epsilon_Lse, distances, \
												   allowed_segments, mu=mu, d=d, kuhnlength=kuhnlength)
			if fe_Lse < minfe[2]:
				minfe = (Ls, Le), rho_Lse, fe_Lse, warning
	return minfe

#calculates estimate for k_B*T, the Boltzmann constant x temperature?
def calc_beta(polymer, epsilon, distances, segments, mu, ndim, kuhnlength, dF=0, xtol=1.e-4):
	def fe_diff(beta, check=False): 
		rho, fe, warning = minimize_fe(polymer, beta * epsilon, distances, segments, \
									   mu=mu, d=ndim, kuhnlength=kuhnlength)
		if check:
			return warning
		else:
			return fe - dF
	def fe_diff_fm(beta, check=False): 
		r, rho, fe, warning = minimize_fe_Lse(polymer, beta * epsilon, distances, segments, \
											  mu, ndim, kuhnlength)
		if check:
			return warning
		else:
			return fe - dF
	def minimize_beta(): #find beta so that minimize_fe(polymer, beta * epsilon, ...) = dF
		minbeta = 1.
		for maxbeta in range(4, 7):
			try:
				beta = scipy.optimize.brentq(fe_diff, minbeta, maxbeta, xtol=xtol)
				if beta == 0:
					minbeta /= 2.
					continue
				return beta
			except ValueError:
				continue
		return np.nan
	def minimize_beta_fm(guess): #find beta so that minimize_fe_Lse(polymer, beta * episilon, ...) = dF
		return scipy.optimize.brentq(fe_diff_fm, 0.95 * guess, 1.05 * guess, xtol=xtol)
	beta = minimize_beta()
	return minimize_beta_fm(beta)

def calc_absolute_energies(polymer, output=None, ndim=3.0, mu=-2.0, kuhnlength=2.0, dFdN=-0.075):
	
	with open(polymer, 'r') as f:
		nresidues = int(f.readline().split()[-1])
		energies = {(int(line.split()[0]), int(line.split()[1])) : -float(line.split()[2]) \
					for line in f}
	polymer = Polymer([edge for edge in energies], nresidues=nresidues)
	for edge in calc_ghost_edges(polymer):
		polymer.add_edge(*edge)
	for u,v in list(energies.keys()):
		energies[(v,u)] = energies[(u,v)]
	epsilon = np.array([[energies[(u,v)] if (u,v) in energies else 0 for u in range(nresidues)] \
						for v in range(nresidues)])
	dF = polymer.number_of_residues() * dFdN
	distances = polymer.distances()
	ranges_cc = native_connected_components(polymer, [(0, polymer.number_of_residues())])

	print("# Polymer:", polymer)
	print("# mu =", mu)
	print("# ndim =", ndim)
	print("# kuhnlength =", kuhnlength)
	print("# dFdN =", dFdN)
	print("# dF =", dF)
	print("# Number of residues:", polymer.number_of_residues())
	print("# Number of contacts:", polymer.number_of_edges())

	betaepsilon = np.zeros(epsilon.shape)
	init_time = time.time()
	for i, ranges in enumerate(ranges_cc):
		fullrange = (min(u for r in ranges for u in r), max(u for r in ranges for u in r))
		epsilon_fullrange = subpolymer_epsilon(epsilon, [u for u in range(fullrange[0], fullrange[1])])
		epsilon_ranges = {r : subpolymer_epsilon(epsilon, [u for u in range(r[0], r[1])]) \
						  for r in ranges} #unused
		epsilon_ranges[fullrange] = epsilon_fullrange #unused
		subpolymer = get_subpolymer(polymer, ranges)
		subpolymer_ranges = {r : get_subpolymer(polymer, [r]) for r in ranges} #unused
		substructures = find_substructures(subpolymer)
		subpolymer_segments = get_polymer_segments(subpolymer, substructures)
		segments_ranges = {r : {ss : [[u for u in seg] for seg in segs] \
								for ss,segs in subpolymer_segments.items() \
								if all(u >= r[0] and u < r[1] for seg in segs for u in seg)} \
						   for r in ranges}
		segments_ranges[fullrange] = {ss : [[u for u in seg] for seg in segs] \
									  for ss,segs in subpolymer_segments.items()}
		beta_full = calc_beta(subpolymer, epsilon_fullrange, distances, segments_ranges[fullrange], \
							  mu, ndim, kuhnlength, dF=dF)
		minbeta_range, minbeta = fullrange, beta_full
		betaepsilon[fullrange[0]:fullrange[1],fullrange[0]:fullrange[1]] \
			+= minbeta * epsilon_fullrange[fullrange[0]:fullrange[1],fullrange[0]:fullrange[1]]
		update_time(i, len(ranges_cc), init_time)

	if output == None:
		output = polymer.replace('.dat', '_abs_energies.dat')
	print("Writing %s" % os.path.join(output, "polymer_abs_energies.dat"))
	with open(os.path.join(output, "polymer_abs_energies.dat"), 'w') as f:
		write_polymer(-betaepsilon, f)

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('polymer', type=str, help="path to input polymer-graph file")
	parser.add_argument('--output', type=str, default=None, \
						help="path to output directory")
	parser.add_argument('--ndim', type=float, default=3., help="number of spatial dimensions [3]")
	parser.add_argument('--mu', type=float, default=-2., help="residue ordering penalty / kT [-2.]")
	parser.add_argument('--kuhnlength', type=float, default=2., metavar='B', help="Kuhn length [2]")
	parser.add_argument('--dFdN', type=float, default=-0.075, \
						help="free-energy of the native state per residue [-0.075]")
	args = parser.parse_args()

	calc_absolute_energies(args.polymer, args.output, args.ndim, args.mu, args.kuhnlength, args.dFdN)

