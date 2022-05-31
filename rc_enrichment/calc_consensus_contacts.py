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

#This script will analyze crystal structures (in pdb format) to find native contacts between amino acids. It then computes an energy function based on the native contacts, hydrogen bonding, and alpha helices between each pair of amino acids. This output is used to construct a native contact graph model of the protein.

import argparse, prody, gzip, subprocess, os.path
import numpy as np
from itertools import product
from Bio import pairwise2
import matplotlib.pyplot as plt
import matplotlib.collections as collections
from tools import *

from codons import codon_to_aa
from polymer import Polymer
from substructures import find_substructures, calc_co
from calc_codon_usage import read_fasta

def aaseq(seq):
	codons = [seq[3*i:3*(i+1)] for i in range(len(seq)//3)]
	return ''.join(codon_to_aa[c] for c in codons if codon_to_aa[c] != 'Stop')

def align_sequences(seqs, clustalo_cmd='/usr/bin/clustalo'):
	with subprocess.Popen([clustalo_cmd, '-i', '-'],
						  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE) \
						  as proc:
		clustalo_input = bytes(''.join(">%d\n%s\n" % (i, seqs[i]) for i in range(len(seqs))), \
							   'ascii')
		clustalo_output, clustalo_err = proc.communicate(clustalo_input)
	aligned_seqs = []
	for line in clustalo_output.decode('ascii').split('\n'):
		if len(line) > 0 and line[0] == '>':
			aligned_seqs.append('')
		elif len(line) > 0:
			aligned_seqs[-1] += line.strip()
	return aligned_seqs

def get_alignment_indices(aligned_seqs):
	indices = []
	for i in range(1, len(aligned_seqs)):
		indices.append({})
		c = 0
		for j in range(len(aligned_seqs[i])):
			if aligned_seqs[i][j] != '-':
				indices[-1][c] = j
				c += 1
	return indices

def get_percent_id(aligned_seqs):
	percid = []
	for i in range(1, len(aligned_seqs)):
		percid.append(sum(1 if aligned_seqs[0][c] == aligned_seqs[i][c] else 0 \
						  for c in range(len(aligned_seqs[0])) if aligned_seqs[0][c] != '-') \
					  / len(aligned_seqs[0].replace('-', '')) * 100.)
	return percid

def get_shifted_residues(chain):
	residues = [r for r in chain.iterResidues() if r.isprotein]
	shifted_residues = [[] for i in range(len(residues))]
	for i in range(len(residues)):
		residue = residues[i]
		for atom in residue.iterAtoms():
			if atom.getElement() == 'H':
				continue
			elif atom.getName() == 'N' and i > 0:
				shifted_residues[i - 1].append(atom)
			else:
				shifted_residues[i].append(atom)
	return [res for res in shifted_residues if len(res) > 0]

def build_contact_map_from_residues(n, indices, residues, cutoff=4., min_sequence_distance=3):
	atoms = np.array([[x for x in atom.getCoords()] for residue in residues \
					  for atom in residue if atom.getElement() != 'H'])
	atomranges = []
	j = 0
	for residue in residues:
		start = j
		for atom in residue:
			if atom.getElement() != 'H':
				j += 1
		stop = j
		atomranges.append((start, stop))
	centers = np.array([np.mean([atoms[ai] for ai in range(*atomranges[i])]) \
						for i in range(len(residues))])
	A = np.zeros((n, n))
	for i in range(len(residues)):
		for j in range(i + min_sequence_distance, len(residues)):
			if np.linalg.norm(centers[i] - centers[j]) <= cutoff * 4:
				for ai in range(*atomranges[i]):
					for aj in range(*atomranges[j]):
						if np.linalg.norm(atoms[ai] - atoms[aj]) <= cutoff:
							A[indices[i],indices[j]] += 1
							A[indices[j],indices[i]] += 1
	return A

def identify_backbone_hydrogen_bonds(contacts, indices, residues, cutoff=4.):
	H = np.zeros((len(contacts), len(contacts)))
	for i in indices:
		for j in indices:
			if contacts[indices[i],indices[j]] > 0:
				for atom1, atom2 in product(residues[i], residues[j]):
					if ((atom1.getName() == 'N' and atom2.getName() == 'O') or \
						(atom1.getName() == 'O' and atom2.getName() == 'N')) and \
						np.linalg.norm(atom1.getCoords() - atom2.getCoords()) <= cutoff:
						H[indices[i],indices[j]] = H[indices[j],indices[i]] = 1
						break
	return H

def consensus_contacts(wtseq, allcontacts, allhbonds, indices, consensus_frac=0.25):
	n = list(allcontacts.values())[0].shape[0]
	wtn = len(wtseq.replace('-', ''))
	contacts = np.zeros((wtn, wtn))
	hbonds = np.zeros((wtn, wtn))
	ii = -1
	for i in range(n):
		if wtseq[i] != '-':
			ii += 1
			jj = -1
			for j in range(n):
				if wtseq[j] != '-':
					jj += 1
					Cij = np.array([allcontacts[pdbid][i,j] for pdbid in allcontacts.keys()])
					Hij = np.array([allhbonds[pdbid][i,j] for pdbid in allcontacts.keys()])
					pdbs_covering = sum(1 if (i in indices[pdbid].values() and j in indices[pdbid].values()) else 0 for pdbid in allcontacts.keys()) #find the number of pdbs covering both i and j (excludes pdbs not covering those positions)
					if pdbs_covering > 0 and sum(1 if Cij[k] > 0 else 0 for k in range(len(Cij))) / pdbs_covering \
					   >= consensus_frac:
						contacts[ii,jj] = np.max(Cij)
						hbonds[ii,jj] = np.max(Hij)
	return contacts, hbonds

def bond_energies(polymer, substructures, contacts, hbonds, alpha_helix=0.625, alpha_hbond=16.):
	ishelix = {}
	for i in range(contacts.shape[0]):
		for j in range(contacts.shape[1]):
			ishelix[(i,j)] = ishelix[(j,i)] = 0
			for ss in substructures:
				if (i,j) in ss or (j,i) in ss:
					if calc_co(ss) < 5:
						ishelix[(i,j)] = ishelix[(j,i)] = 1
						break
	U = np.zeros(contacts.shape)
	for i in range(contacts.shape[0]):
		for j in range(contacts.shape[1]):
			U[i,j] = alpha_helix**ishelix[(i,j)] * (contacts[i,j] / alpha_hbond + hbonds[i,j])
	return U

def write_polymer(U, stream):
	stream.write("nresidues = %d\n" % U.shape[0])
	for i in range(U.shape[0]):
		for j in range(i, U.shape[1]):
			if U[i,j] > 0:
				stream.write("%d %d %g\n" % (i, j, U[i,j]))

def select_chains(fasta, pdbids, structure_path, focus="", excluded_chains="", all_chains=False):

	if excluded_chains is not "" and excluded_chains is not None:
		#process excluded chains
		chain_list = []
		excluded_chains = (excluded_chains.upper()).split(",")
		for chains in excluded_chains:
			pdbid, chainids = chains.split("_")
			for chainid in chainids.strip():
				chain_list.append(pdbid + "_" + chainid)
		excluded_chains = chain_list
	else:
		excluded_chains = []

	seqs = read_fasta(fasta)
	if focus == None or len(focus) < 1:
		focus = list(seqs.keys())[0]
	wtseq = aaseq(seqs[focus])

	chains = {}
	for pdbid in pdbids.split(','):
		pdbid = pdbid.lower()
		try:
			pdbfile = structure_path + '%s.pdb.gz' % pdbid
			f = gzip.open(pdbfile, 'rt')
		except:
			pdbfile = structure_path + '%s.pdb' % pdbid
			f = open(pdbfile, 'rt')
		try:
			mol, header = prody.parsePDBStream(f, header=True)
		except ValueError as e:
			print("PDB Error: %s: %s" % (pdbid, e))
			continue
		''' removed due to possibility of NMR and EM structures
		if 'experiment' not in header or 'X-RAY' not in header['experiment']:
			print("%s is not an X-ray structure; skipping" % pdbfile)
			continue
		'''
		if header != None and 'resolution' in header:
			resolution = float(header['resolution'])
		else:
			resolution = 100.
		pdbchains = {}
		for chain in mol.iterChains():
			#determine whether or not chain is excluded
			chainid = pdbid + "_" + str(chain).split()[-1]
			if chainid.upper() in excluded_chains:
				continue
			
			#score the chain against the sequence of interest
			chainseq = chain.getSequence()
			score = pairwise2.align.globalms(wtseq, chainseq, 2, -1, -0.25, -.1, score_only=True)
			residues = get_shifted_residues(chain)
			if len(residues) != len(chainseq):
				raise Exception("mismatched lengths")
			if isinstance(score, (float, int)):
				pdbchains[str(chain).split()[-1]] = (chainseq, residues, score)

		bestchain = max([k for k in pdbchains.keys() if is_numeric(pdbchains[k][2])], key=lambda y: pdbchains[y][2])
		chains[pdbid + "_" + bestchain] = pdbchains[bestchain] #collect only the best matching chain from the pdb

		if all_chains: #collect all the chains with a positive score and sort them out later
			for chain, value in pdbchains.items():
				if chain is not bestchain and value[2] >= 0:
					chains[pdbid + "_" + chain] = value
		f.close()
	if len(chains) > 100:
		chains = {pdbid : chains[pdbid] for pdbid in sorted(chains, key=lambda y: chains[y][2], reverse=True)[:100]}

	pdbids = sorted(chains.keys())
	if len(pdbids) == 0:
		print("No chains found.")
		raise SystemExit
	print([wtseq] + [chains[pdbid][0] for pdbid in pdbids])
	aligned_chains = align_sequences([wtseq] + [chains[pdbid][0] for pdbid in pdbids])
	aligned_indices = get_alignment_indices(aligned_chains)

	return (wtseq, pdbids, chains, aligned_chains, aligned_indices)

#returns the ranges where the pdbs match the sequence
def get_alignment_profile(aligned_seqs):
	return_ranges = []
	for i in range(1, len(aligned_seqs)):

		#look for segments where the pdb matches the sequence
		ranges = []
		segment = []
		for c in range(len(aligned_seqs[0])):
			if aligned_seqs[0][c] == aligned_seqs[i][c] and aligned_seqs[0][c] != '-': #segment matching the sequence
				if segment == []: #start of matching segment
					segment = [c]
			else:	#segment doesn't match the sequence
				if segment != []: #end the matching segment
					segment.append(c-1)
					ranges.append(segment)
				segment = []

		#if the sequence ended, but a new matching segment was started, finish it
		if segment != []:
			segment.append(len(aligned_seqs[0]))
			ranges.append(segment)
		return_ranges.append(ranges)

	return return_ranges

#determines whether i is in any of the ranges
def contained(i, ranges):
	for range_vals in ranges:
		if i >= range_vals[0] and i <= range_vals[1]:
			return 1
	return 0

#determines whether or not r1 and r2 have any overlap
def overlap(r1, r2):
	if (r1[0] >= r2[0] and r1[0] <= r2[1]) or (r1[1] >= r2[0] and r1[1] <= r2[1]) or (r2[1] >= r1[0] and r2[1] <= r1[1]) or (r2[1] >= r1[0] and r2[1] <= r1[1]):
		return([min(r1[0], r1[0]), max(r1[1], r2[1])])
	else: 
		return None

#useful for combining alignment profiles for different chains of the same protein
def combine_ranges(set1, set2):
	ranges = set1 + set2
	for i, r1 in enumerate(sorted(set1, key=lambda tup: tup[0])):
		for j, r2 in enumerate(sorted(set2, key=lambda tup: tup[0])):
			new_r = overlap(r1, r2)
			if new_r is not None:
				try:
					ranges.remove(r1)
					ranges.remove(r2)
				except:
					pass
				ranges.append(new_r)
	return(ranges)

def plot_pdbs(chainids=[], aligned_chains=[], output='fig_pdbs'):
	
	#find the ranges where each pdb matches the sequence
	alignment_profile = get_alignment_profile(aligned_chains)

	#combine identity ranges for different chains of the same pdb
	combined_chains = {}
	for i, chainid in enumerate(chainids):
		pdb_name = chainid.split('_')[0]
		if pdb_name not in combined_chains:
			combined_chains[pdb_name] = alignment_profile[i]
		else:
			combined_chains[pdb_name] = combine_ranges(combined_chains[pdb_name], alignment_profile[i])

	#plot the data
	fig = plt.figure(figsize = (len(aligned_chains[0])/6.0, len(chainids)))
	x = np.arange(0,len(aligned_chains[0])-1,1)
	
	for i, pdb in enumerate(combined_chains):
		#compute the places where the pdb matches the sequence
		y = [contained(x_i,combined_chains[pdb]) for x_i in x]
		booly = [y_i > 0 for y_i in y]
		y = [0 for x_i in x]
	
		#add the plot for the current pdb
		ax = fig.add_subplot(len(combined_chains.keys()),1,i+1)
		ax.set_title(pdb, fontsize=50, verticalalignment='top', horizontalalignment='left')
		ax.plot(x, y, color='white')
		ax.axhline(-1, color='black', lw=2)
		ax.axhline(-1, color='black', lw=2)
		ax.fill_between(x, y1=-1, y2=1, where=booly, facecolor='blue', alpha=0.5)

	plt.xticks(x, aligned_chains[0], rotation='horizontal')
	plt.savefig(output, bbox_inches="tight")

def calc_consensus_contacts(fasta, pdbids, output, structure_path, focus="", plot=True, min_percent_id=95, consensus_frac=0.25, excluded_chains="", all_chains=False):
	wtseq, pdbids, chains, aligned_chains, aligned_indices = select_chains(fasta, pdbids, structure_path, focus=focus, excluded_chains=excluded_chains, all_chains=all_chains)

	if plot:
		plot_pdbs(pdbids, aligned_chains, output=os.path.join(output, "fig_pdb"))

	contact_maps = {}
	hbonds = {}
	for i in range(len(pdbids)):
		indices = aligned_indices[i]
		residues = chains[pdbids[i]][1]
		contact_maps[pdbids[i]] = build_contact_map_from_residues(len(aligned_chains[0]), indices, residues)
		hbonds[pdbids[i]] = identify_backbone_hydrogen_bonds(contact_maps[pdbids[i]], indices, residues)
	percid = get_percent_id(aligned_chains)
	for i in range(len(pdbids)):
		print("%s %6.2f" % (pdbids[i], percid[i]))
	selected_pdbids = [pdbids[i] for i in range(len(pdbids)) if percid[i] >= min_percent_id]
	if len(selected_pdbids) == 0:
		print("No chains selected.")
		raise SystemExit
	print("Using %d chains with percentid >= %g." % (len(selected_pdbids), min_percent_id))

	contacts, hbonds = consensus_contacts(aligned_chains[0], \
										  {pdbid : contact_maps[pdbid] for pdbid in selected_pdbids}, \
										  {pdbid : hbonds[pdbid] for pdbid in selected_pdbids}, \
										  {pdbids[i]: aligned_indices[i] for i in range(len(pdbids))}, \
										  consensus_frac=consensus_frac)
	polymer = Polymer([(i,j) for i in range(contacts.shape[0]) for j in range(contacts.shape[1]) \
					   if contacts[i,j] > 0])
	substructures = find_substructures(polymer)
	U = bond_energies(polymer, substructures, contacts, hbonds)
	print("Number of residues:", contacts.shape[0])
	print("Minimum consensus fraction:", consensus_frac)
	print("Total bond energy:", 0.5 * U.sum())

	print("Writing %s" % os.path.join(output, "polymer.dat"))
	with open(os.path.join(output, "polymer.dat"), 'w') as f:
		write_polymer(U, f)

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('fasta', type=str, help="path to fasta sequence")
	parser.add_argument('pdbids', type=str, help="comma-separated list of pdbids")
	parser.add_argument('--output', type=str, metavar='PATH', default='./', \
						help="path to output directory [./]")
	parser.add_argument('--structure-path', type=str, metavar='PATH', default='./', \
						help="path to gzip'ed input pdb files [./]")
	parser.add_argument('--excluded-chains', type=str, default="", help="comma-separated list of chains to exclude (formatted as PDBID_CHAIN, e.g. 2wpm_S)")
	parser.add_argument('--min-percent-id', type=float, default=95., \
						help="minimum percent sequence identity [95.]")
	parser.add_argument('--consensus-frac', type=float, default=0.25, \
						help="minimum contact consensus fraction [0.25]")
	parser.add_argument('--plot_pdbs', action="store_true", help="whether or not to plot the coverage of the pdbs")
	parser.add_argument('--all-chains', action="store_true", help="whether or not to collect all chains in the pdb, or only the best chain")
	args = parser.parse_args()

	calc_consensus_contacts(args.fasta, args.pdbids, args.output, args.structure_path, plot=args.plot_pdbs, min_percent_id=args.min_percent_id, consensus_frac=args.consensus_frac, excluded_chains=args.excluded_chains, all_chains=args.all_chains)
