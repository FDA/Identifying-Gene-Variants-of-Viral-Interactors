import random, re, os, sys, os.path, copy, warnings, math, time, scipy, warnings,  pickle, statistics, copy, tools, itertools, Bio
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, DistanceMatrix
import sklearn
import numpy as np
import pandas as pd
import align_codons

#Unfinished
def process_clust(n_points, fit_cluster):
	ii = itertools.count(n_points)
	clusters = [{'node_id': next(ii), 'left': x[0], 'right':x[1]} for x in fit_cluster.children_]
	members = {i:[i] for i in range(n_points)}
	for cluster in clusters:
		node_id = cluster["node_id"]
		members[node_id] = copy.deepcopy(members[cluster["left"]])
		members[node_id].extend(copy.deepcopy(members[cluster["right"]]))

	on_split = {c["node_id"]: [c["left"], c["right"]] for c in clusters}
	up_merge = {c["left"]: {"into": c["node_id"], "with": c["right"]} for c in clusters}
	up_merge.update({c["right"]: {"into": c["node_id"], "with": c["left"]} for c in clusters})

	return(members, on_split)

#Input: multiple sequence alignment (seqs), a focus sequence (genename), and a substitution matrix (submatrix),
#Output:a phylogeny tree based on pairwise distances
#Note: by default, uses Levenshtein distance for the distance matrix
def gen_tree(seqs, genename, submatrix=None, space="-", ram_disk="/media/temp/"):
	if submatrix != None and isinstance(submatrix, dict):
		piecelen = len(list(submatrix.keys())[0])
		maxpenalty = min([submatrix[k1][k2] for k1 in submatrix for k2 in submatrix[k1]])
	else:
		piecelen = 1
		maxpenalty = 0
	if genename == None or len(genename) < 1:
		genename = list(seqs.keys())[-1]
	if isinstance(seqs, (dict,)):
		names = list(seqs.keys())
	elif isinstance(seqs, (list, tuple)):
		names = [str(i) for i in range(len(seqs))]
	seqs = [seqs[name] for name in names]
	n = len(seqs)
	try: # assume user provided a substitution matrix that can be used by Biopython's Phylo library
		calculator = Bio.Phylo.TreeConstruction.DistanceCalculator(submatrix)
		tools.write_fasta({names[i]:seqs[i] for i in range(len(seqs))}, os.path.join(ram_disk, genename + "_tmp.fasta"))
		aln = Bio.AlignIO.read(os.path.join(ram_disk, genename + "_tmp.fasta"), "fasta")
		distance = calculator.get_distance(aln)

		combined = "".join(seqs)
		alphabet = sorted(set([combined[i] for i in range(len(combined))]))
		submatrix = {m:{n:np.inf for n in alphabet} for m in alphabet}
		for m in alphabet:
			for n in alphabet:
				tools.write_fasta({"m":m, "n":n}, os.path.join(ram_disk, "tmp_distance.fasta"))
				aln = Bio.AlignIO.read(os.path.join(ram_disk, "tmp_distance.fasta"), "fasta")
				tmp_dist = calculator.get_distance(aln)
				submatrix[m][n] = tmp_dist[1][0]
	except Exception as e: 
		distance = [[0 for j in range(len(seqs))] for i in range(len(seqs))]
		init_time = time.time()
		for i in range(len(distance)):
			seq1 = seqs[i]
			for j in range(len(distance[i])):
				seq2 = seqs[j]	
				if isinstance(submatrix, (dict, pd.DataFrame)): #submatrix is usable as a substitution matrix
					for k in range(len(seq1)//piecelen):
						if piecelen*(k+1) <= len(seq1) and seq1[piecelen*k:piecelen*(k+1)] in submatrix and seq2[piecelen*k:piecelen*(k+1)] in submatrix[seq1[piecelen*k:piecelen*(k+1)]]:
							distance[i][j] -= submatrix[seq1[piecelen*k:piecelen*(k+1)]][seq2[piecelen*k:piecelen*(k+1)]]
						elif seq1[piecelen*k:piecelen*(k+1)] == seq2[piecelen*k:piecelen*(k+1)]:
							pass
						else:
							distance[i][j] -= maxpenalty
				else: #default to Levenshtein distance
					distance[i][j] = tools.levenshtein(seq1.replace(space, ""), seq2.replace(space, ""))
				tools.update_time(i*n+j, n**2, init_time)
	distance = [[(distance[i][j] + distance[j][i])/2 for j in range(i+1)] for i in range(len(seqs))]
	minval = min([distance[i][j] for i in range(len(distance)) for j in range(len(distance[i]))])
	maxval = max([distance[i][j] for i in range(len(distance)) for j in range(len(distance[i]))])
	distance = [[(distance[i][j] - minval)/(maxval - minval) for j in range(len(distance[i]))] for i in range(len(distance))]
	t = DistanceTreeConstructor().nj(DistanceMatrix(names, matrix=distance))
	return(t, submatrix)
	
#unfinished
#attempts to recreate rate4site for general purpose sequences (including NTs and codons)
def rate4site2(seqs, submat, focus=None, alphabet=None, tree=None):
	if focus == None or focus not in seqs.keys():
		focus = list(seqs.keys())[-1]
	if alphabet == None:
		alphabet = tools.infer_alphabet("".join(list(seqs.values())))
	if isinstance(submat, str):
		try:
			submat = pd.read_csv(submat, sep="\t", index_col=0).T.to_dict()
		except Exception as e:
			pass
	if isinstance(submat, pd.DataFrame):
		submat = submat.T.to_dict()
	if isinstance(submat, dict):
		piecelen = len(list(submat.keys())[0])
	else:
		piecelen = 1
	if any([True if len(seqs[k]) != len(seqs[focus]) else False for k in seqs.keys()]):
		seqs = tools.align_sequences(seqs)
	if tree == None or not isinstance(submat, dict):
		tree, submat = gen_tree(seqs, focus, submatrix=submat)
	distance = {m:{n:tree.distance(m, n) for n in seqs.keys()} for m in seqs.keys()}
	minval = min([distance[k][k2] for k in seqs.keys() for k2 in seqs.keys()])
	maxval = max([distance[k][k2] for k in seqs.keys() for k2 in seqs.keys()])
	distance = pd.DataFrame({(m,n):{"distance":(distance[m][n] - minval)/(maxval - minval)} for n in seqs.keys() for m in seqs.keys()}).T
	names = list(seqs.keys())
	if isinstance(submat, dict):
		minval = min([submat[k][k2] for k in submat.keys() for k2 in submat[k].keys()])
		maxval = max([submat[k][k2] for k in submat.keys() for k2 in submat[k].keys()])
		for k,v in submat.items():
			for k2, v2 in v.items():
				submat[k][k2] = (submat[k][k2] - minval)/(maxval - minval)
	jrange = [i for i in range(len(seqs[focus])//piecelen) if seqs[focus][piecelen*i:piecelen*(i+1)] in submat]
	seqmat = {(m,n):{i:1 for i in jrange} for n in seqs.keys() for m in seqs.keys()}
	for i, j in enumerate(jrange):
			for m in seqs.keys():
				for n in seqs.keys():
					try:
						if seqs[m][piecelen*j:piecelen*(j+1)] != seqs[n][piecelen*j:piecelen*(j+1)]:
							seqmat[(m, n)][j] = submat[seqs[m][piecelen*j:piecelen*(j+1)]][seqs[n][piecelen*j:piecelen*(j+1)]]
						else:
							seqmat[(m, n)][j] = 0
					except KeyError as e:
						seqmat[(m, n)][j] = 1
	seqmat = pd.DataFrame(seqmat).T
	r = scipy.linalg.lstsq(seqmat, distance)[0]
	r = [abs(xi[0]) for xi in r]
	return(r)
