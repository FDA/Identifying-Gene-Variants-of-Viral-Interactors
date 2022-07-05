import os, sys, re, subprocess, math, time, statistics, pickle, requests, distance
import random
import numpy as np
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.cluster

#Input: DNA sequence
#Output: complement sequence
def complement(seq):
	d = {"C":"G", "G":"C", "A":"T", "T":"A"}
	if all([True if seq[i] in d.keys() else False for i in range(len(seq))]):
		return("".join([d[seq[i]] for i in range(len(seq))]))

#Input: a list of values, a key
#Output: percentile rank of key relative to list l
def pctile(l, k):
	l = [float(li) for li in list(l) if is_numeric(li) and not np.isnan(float(li))]
	if len(l) >= 0 and is_numeric(k) and not np.isnan(float(k)):
		return(len([li for li in l if li <= float(k)])/len(l))

#Input: string representing HTML/XML formatted site
#Output: possible elements matching the given key ("tr" by default)
#finds an HTML/XML element when regex doesn't work
def get_element(s, k="tr"):#find HTML tagged areas (regex doesn't work)
	elements = []
	while True:
		start = s.find('<' + k)
		if start < 0:
			break
		end = s.find("/" + k + ">")
		if end < 0:
			break
		elements.append(s[start:end+len(k)+2])
		s = s[end+len(k)+2:]
	return([x for x in elements if len(x) > 0])

#Input: a complex dictionary structure d, a key k
#Output: any value matching the key k in d
def seek_in_struct(d, k):
	return_d = {}
	if isinstance(d, dict):
		if any([True if k[i] in d.keys() else False for i in range(len(k))]):
			return({k[i]:d[k[i]] for i in range(len(k)) if k[i] in d.keys()})
		else:
			for k2,v2 in d.items():
				return_d.update(seek_in_struct(v2, k))	
	elif isinstance(d, (tuple, list)):
		return_d = {}
		for i in range(len(d)):
			return_d.update(seek_in_struct(d[i], k))
	return(return_d)


#Input: strings s1, s2
#Output: the Levenshtein distance between s1 and s2 (edit distance)
def levenshtein(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]

#Input: dictionary of k:count or frequency values, number of samples to return
#Output: a sample of keys from dictionary, sampled proportionally to their values
def custom_sample(d, k=None):
	for k2, v2 in d.items():
		if np.isnan(v2):
			d[k2] = 0.0
	total = sum(d.values())
	if total != 1.0:
		d = {k2:v2/total for k2,v2 in d.items()}
	try:
		return(np.random.choice(list(d.keys()), p=list(d.values()), size=k))
	except Exception as e:
		return(None)

#Input: a function, any positional and keyword arguments of function
#Output: outcome of function with arguments, retried multiple times if issues are reached
def retry_func(func, positional_arguments=[], keyword_arguments={}, lim=10, wait=2):
	for i in range(lim):
		try:
			return(func(*positional_arguments, **keyword_arguments))
		except Exception as e:
			time.sleep(wait*i)
	warnings.warn("\033[93m" + str(func.__name__) + " failed after " + str(lim) + " tries." + "\033[0m") 
	return(None)

#Input: sequence, position in sequence, and length of desired subsequence
#Output: subsequence of length l, centered on position pos (l/2 characters on either side of position)
def subseq(seq, pos, l):
	return(seq[max(0, pos - math.floor(l/2)):min(len(seq), pos + math.ceil(l/2))])

#Input: multiple sequence alignment, and name of a focus gene
#Output: uniqueness measure of sequences
def seq_unique(seqs, genename=None, space="-"):
	if genename == None:
		genename = list(seqs.keys())[-1]
	names = list(seqs.keys())
	uniqueness = {}
	identmat = [[-1 for seq1 in seqs] for seq2 in seqs]
	init_time = time.time()
	n = len(seqs)
	for i, (name1, seq1) in enumerate(seqs.items()):
		l = len(seq1.replace(space, ""))
		unique = len(seq1)
		for j, (name2, seq2) in enumerate(seqs.items()):
			if name1 == name2 or j == i:
				continue
			identmat[i][j] = len([i for i in range(len(seq1)) if (seq1[i] == seq2[i] and seq1[i] != space)])
			if j < i and identmat[j][i] >= 0:
				identmat[i][j] = identmat[j][i]
			update_time(i * n + j, n ** 2, init_time)
		uniqueness[name1] = 1-max([x/l for x in identmat[i] if x >= 0])
	identmat = [[identmat[i][j]/len(seqs[names[i]].replace(space, "")) for i in range(len(seqs))] for j in range(len(seqs))]
	toremove = [k for k in seqs.keys() if k != genename]
	for i in range(len(names)):
		if not any([True for j in range(len(identmat[i])) if identmat[i][j] == 1 and i != j]):
			try:
				toremove.remove(names[i])
			except ValueError as e:
				pass
		else:
			longest = max([names[j] for j in range(len(identmat[i])) if i != j and identmat[i][j] == 1], key=lambda kv: len(seqs[kv].replace(space, "")))
			try:
				toremove.remove(longest)
			except ValueError as e:
				pass
	seqs_wo_dups = {k:seqs[k] for k in names if k not in toremove}
	total = sum(uniqueness.values())
	uniqueness = {k:v/total*len(seqs) for k,v in uniqueness.items()}
	identmat = [[identmat[i][j] for i in range(len(names)) if names[i] not in toremove] for j in range(len(names)) if names[j] not in toremove]
	return(uniqueness, seqs_wo_dups, identmat)

#Inputs: multiple sequence alignment, name of focus gene
#Optional input: distance function between sequences, default is Levenshtein
#Output: similarity between sequences
def seqs_compute_similarity(seqs, genename=None, space="-", func=lambda k: levenshtein(k[0], k[1])):
	if isinstance(seqs, (dict,)):
		names = list(seqs.keys())
	elif isinstance(seqs, (list, tuple)):
		names = [i for i in range(len(seqs))]
	seqs = [seqs[name] for name in names]
	similarity = -1*np.array([[func((seq1.replace(space, ""), seq2.replace(space, ""))) for seq1 in seqs] for seq2 in seqs])
	return(similarity)

#Input: similarity matrix of sequences, multiple sequence alignment
#Output: cluster centers when sequences are clustered by similarity using affinity propagation
def cluster_seqs(similarity, seqs, genename=None, space="-"):
	if isinstance(seqs, (dict,)):
		names = list(seqs.keys())
	elif isinstance(seqs, (list, tuple)):
		names = [i for i in range(len(seqs))]
	seqs = [seqs[name] for name in names]
	affprop = sklearn.cluster.AffinityPropagation(affinity="precomputed", damping=0.5)
	affprop.fit(similarity)
	newseqs = {}
	for cluster_id in np.unique(affprop.labels_):
		exemplar = seqs[affprop.cluster_centers_indices_[cluster_id]]
		cluster = np.unique([names[i] for i in range(len(names)) if affprop.labels_[i]==cluster_id])
		cluster_str = ", ".join(cluster)
		#print(" - *%s:* %s" % (exemplar, cluster_str))
		newseqs[names[affprop.cluster_centers_indices_[cluster_id]]] = exemplar
	newseqs[genename] = seqs[names.index(genename)]
	return(newseqs)

#Input: similarity matrix of sequences, multiple sequence alignment
#Output: cluster centers when sequences are clustered by similarity using affinity propagation, also cluster labels:sequences in cluster dictionary
def cluster_seqs2(seqs, genename=None, space="-", frac=False, mode="identity"):
	if genename == None or len(genename) < 1:
		genename = list(seqs.keys())[-1]
	if frac: #compute uniqueness and remove duplicate sequences
		uniqueness, seqs, identmat = seq_unique(seqs, genename=genename, space=space)
		if not isinstance(frac, (bool,)) and isinstance(frac, (int, float)) and mode.lower() in ["levenshtein", "leven", "ls"]:
			if 0 <= frac and frac <= 1:
				frac = 100 * frac
			elif frac > 100 or frac < 0:
				raise Exception("Can't interpret percentile value " + str(frac))
			pctile = np.percentile(list(uniqueness.values()), 100 - frac)
			seqs = {k:v for k,v in seqs.items() if k == genename or (k in uniqueness and uniqueness[k] > pctile)}
	if isinstance(seqs, (dict,)):
		names = list(seqs.keys())
	elif isinstance(seqs, (list, tuple)):
		names = [i for i in range(len(seqs))]
	seqs = [seqs[name] for name in names]
	lev_similarity = np.array([[1 for seq1 in seqs] for seq2 in seqs])
	init_time = time.time()
	n = len(seqs)
	for i, seq1 in enumerate(seqs):
		for j, seq2 in enumerate(seqs):
			if mode.lower() not in ["levenshtein", "leven", "ls"] and 'identmat' in locals():
				lev_similarity[i][j] = 1-identmat[i][j]
			else:
				if j < i and lev_similarity[j][i] < 0:
					lev_similarity[i][j] = lev_similarity[j][i]
				else:
					lev_similarity[i][j] = -levenshtein(seq1.replace(space, ""), seq2.replace(space, ""))
			update_time(i * n + j, n ** 2, init_time)
	affprop = sklearn.cluster.AffinityPropagation(affinity="precomputed", damping=0.9)
	affprop.fit(lev_similarity)
	clusters = {i:[names[j] for j in range(len(names)) if affprop.labels_[j] == i] for i in set(affprop.labels_)}
	newseqs = {}
	for cluster_id in np.unique(affprop.labels_):
		exemplar = seqs[affprop.cluster_centers_indices_[cluster_id]]
		cluster = np.unique([names[i] for i in range(len(names)) if affprop.labels_[i]==cluster_id])
		cluster_str = ", ".join(cluster)
		#print(" - *%s:* %s" % (exemplar, cluster_str))
		newseqs[names[affprop.cluster_centers_indices_[cluster_id]]] = exemplar
	newseqs[genename] = seqs[names.index(genename)]
	return(newseqs, clusters)

#Input: similarity matrix, clustering from sklearn that can be fitted on matrix
#Output: cluster centers for clustering computed from similarity matrix
def cluster_centers(similarity, cluster, names=[]):
	if names == None or len(names) < 1:
		names = [i for i in range(len(similariy))]
	cluster.fit(-np.array(similarity))
	labels = cluster.predict(-np.array(similarity))
	scores = cluster.transform(-np.array(similarity))
	distances = {i:statistics.mean([(1-similarity[i][j])**2 for j in range(len(labels)) if labels[i] == labels[j]]) for i in range(len(labels))}
	clusters = {names[min([i for i in range(len(labels)) if labels[i] == labels[j]], key=lambda kv: distances[kv])]:[names[i] for i in range(len(labels)) if labels[i] == labels[j]] for j in range(len(labels))}
	return(clusters)

#Input: dictionary of sequences, name of focus sequences
#Output: clustering sequences by similarity
def cluster_seqs_n(seqs, focus=None, n=150):
	if focus == None or len(focus) < 1 or focus not in seqs:
		focus = list(seqs.keys())[-1]
	uniqueness, seqs_wo_dups, identmat = seq_unique(seqs, focus)
	names = list(seqs_wo_dups.keys())
	for i, k in enumerate(names):
		if seqs_wo_dups[k].replace("-", "") == seqs[focus].replace("-", ""):
			names[i] = focus
			break
	kmeans = sklearn.cluster.KMeans(n_clusters=n)
	clusters = cluster_centers(identmat, kmeans, names=names)
	return(clusters)
	
#Input: multiple sequence alignment, name of focus gene, weighting of sequences
#Output: weighted distribution of characters (NT, AA) at each position of focus sequence.
def dist_with_weight(seqs, genename, weight={}, space="-", alphabet=""):
	if set(weight.keys()) != set(seqs.keys()):
		weight, seqs, identmat = seq_unique(seqs, genename=genename, space=space)
	names = list(seqs.keys())
	joined_seqs = "".join(seqs.values()).replace(space, "")
	alphabet = sorted(set([alphabet[i] for i in range(len(alphabet))] + [joined_seqs[i] for i in range(len(joined_seqs))]))
	alphabetsize = len(alphabet)
	#print(seqs)
	inverted_seqs = {i:"".join([seqs[name][i] for name in names]) for i in range(len(seqs[genename]))}
	
	location = 1
	distribution = {}
	for i in range(len(seqs[genename])):
		if seqs[genename][i] != space:
			distribution[location] = {alpha:0 for alpha in alphabet}
			for j in range(len(inverted_seqs[i])):
				if inverted_seqs[i][j] != space:
					distribution[location][inverted_seqs[i][j]] += weight[names[i]]
			total = sum(distribution[location].values())
			try:
				distribution[location] = {k:v/float(total) for k,v in distribution[location].items()}
			except ZeroDivisionError as e:
				pass
			location += 1
	return(distribution)

#Input: RefSeq accession ID, genetic variant, identifer ("c." for coding sequence variant, "g." for genomic variant), assembly to use
#Output: any other coordinates of variant from Mutalyzer
def get_mutalyzer(refid, var, c="c.", assembly="GRCh37"):
	var = get_mutation_data(var)
	version = int(refid.split(".")[-1])
	refid = "".join(refid.split(".")[:-1])
	ids = []
	for i in range(version, 0, -1): #retry until you get the right version
		address = "https://mutalyzer.nl/position-converter?assembly_name_or_alias=" + assembly + "&description=" + refid + "." + str(i) + "%3A" + c + str(var[0]) + var[1][0] + "%3E" + var[1][1]
		r = retry_request(requests.get, [address])
		ids += re.findall("[A-Za-z]+\_\d+\.?\d*\:[a-z]\.[\-|\+]?\d+[\-|\+]?\d*[A-Za-z]\&gt\;[A-Za-z]", r.text)+re.findall("[A-Za-z]+\_\d+\.?\d*\:[a-z]\.\-?\d+\-?\d*[A-Za-z]\>[A-Za-z]", r.text)
		ids = list(set([i.replace("&gt;", ">") for i in ids]))
	return(ids)

#Input: Refseq accession ID, genetic variant, identifer ("c." for coding sequence variant, "g." for genomic variant)
#Output: Returns anay errors or warnings from Mutalyzers
def check_variant(refid, var, c="c."):
	var = get_mutation_data(var)
	address = "https://mutalyzer.nl/name-checker?description=" + str(refid) + "%3A" + str(c) + str(var[0]) + str(var[1][0]) + "%3E" + str(var[1][1])
	r = retry_request(requests.get, [address])
	errors = re.findall("(\d+) Error", r.text)
	warnings = re.findall("(\d+) Warnings", r.text)
	return([int(x) for x in errors], [int(x) for x in warnings])

#Input: Requests function, positional and keyword arguments
#Output: output of function applied on arguments
#retries multiple times if errors are encountered
def retry_request(func, positional_arguments=[], keyword_arguments={}, lim=10, wait=2):
	for i in range(lim):
		try:
			response = func(*positional_arguments, **keyword_arguments)
			response.raise_for_status()
			return(response)
		except:
			time.sleep(wait*i)
	warnings.warn("\033[93m" + str(func.__name__) + " failed after " + str(lim) + " tries." + "\033[0m") 
	return(None)

#Input: sequence, possible alphabet, cuboff fraction
#Output: alphabet of input (NT or amino acid)
def infer_alphabet(seq, alphabet="", cutoff=0.95):
	options = {"nt":["A", "C", "G", "T", "U", "-", "."], "aa":["G", "P", "A", "V", "L", "I", "M", "C", "F", "Y", "W", "H", "K", "R", "Q", "N", "E", "D", "S", "T", "*", "-", "."]}
	combined = []
	for k,v in options.items():
		combined += v
	alphabet = list(set([(seq+alphabet)[i] for i in range(len(seq)+len(alphabet))] + combined))
	frequencies = {alphabet[i]:float(len([j for j in range(len(seq)) if seq[j] == alphabet[i]]))/len(seq) for i in range(len(alphabet))}
	for alpha, letters in sorted(options.items(), key=lambda kv: len(kv[1])):
		if sum([frequencies[letters[i]] for i in range(len(letters))]) >= cutoff:
			return(alpha)
	return(None)

#Input: a list of values
#Output: whether or not list is mostly numeric (also excludes NaN values)
#determines if a list is mostly numeric (useful for lists and dataframes)
def infer_numeric(l, cutoff=0.9, forbidden=["", np.nan, np.inf, None]):
	l = [l[i] for i in range(len(l)) if l[i] not in forbidden]
	if np.nan in forbidden:
		l = [l[i] for i in range(len(l)) if not is_numeric(l[i]) or not np.isnan(float(l[i]))]
	if np.inf in forbidden:
		l = [l[i] for i in range(len(l)) if not is_numeric(l[i]) or not np.isinf(abs(float(l[i])))]
	l = [l[i] if not is_numeric(l[i]) else float(l[i]) for i in range(len(l))]
	if len(l) < 1 or not len([k for k in l if is_numeric(k) and not np.isnan(float(k))])/len(l) > cutoff:
		return(False)
	else:
		return(True)

#one-to-three letter amino acid codes
residue_abbreviations = \
{
    'G' : 'gly', # nonpolar
    'A' : 'ala',
    'V' : 'val',
    'L' : 'leu',
    'I' : 'ile',
    'M' : 'met',
    'F' : 'phe',
    'W' : 'trp',
    'P' : 'pro',
    'S' : 'ser', # polar
    'T' : 'thr',
    'C' : 'cys',
    'Y' : 'tyr',
    'N' : 'asn',
    'Q' : 'gln',
    'D' : 'asp', # negative
    'E' : 'glu',
    'K' : 'lys', # positive
    'R' : 'arg',
    'H' : 'his'
}
residue_abbreviations_inv = {v : k for k,v in residue_abbreviations.items()}

#standard genetic code reverse-translation and translation
codon_table = {'I' : ['ATT', 'ATC', 'ATA'],
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
               'Stop' : ['TAA', 'TAG', 'TGA'],
				'-' : ['---']}
codon_to_aa = {codon : aa for aa,v in codon_table.items() for codon in v}

#list of codons
codons = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"]

#Input: string
#Output: longest numeric substring
def get_num(s):
	num = float(max(re.findall("\d+\.?\d*", s), key=len))
	if num - int(num) == 0.0:
		return(int(num))
	else:
		return(num)

#Input: string
#Output: web-formatted string
#useful when sending data online which must be formatted
def slugify(value):
	value = str(re.sub('[^\w\s-]', '', value).strip().lower())
	value = str(re.sub('[-\s]+', '-', value))
	return(value)

#Input: path to pickle file
#Output: the memory object corresponding to the file
def pickleload(f):
	with open(f, "rb") as pklfile:
		d = pickle.load(pklfile)
	return(d)

#Input: any memory object (variable), path to pickle file
#Output: writes object to disk
def pickledump(d, f):
	with open(f, "wb") as pklfile:
		pickle.dump(d,pklfile)

#Input: string, character, position in string
#Output: mutant sequence from replacing string at position pos with character c
def update_str(s, c, pos):
	return(s[:pos] + c + s[pos+1:])

#Input: iteration count i, total number of iterations n, and initial time
#Output: prints progress and estimated remaining time
def update_time(i, n, init_time, func_name="", linelen=80):
	if i+1 >= n:
		print(" " * (70 + len(func_name)), end="\r", flush=True)
	else:
		elapsed = time.time() - init_time
		estimated = int(elapsed/(i+1) * (n-(i+1)))
		s = padstring(str(round(100*(i+1)/n, 2)), n=5, space="0") + "% done"
		if func_name != None and len(func_name) < 1:
			s += ". "
		else:
			s += " with " + func_name + ". "
		s += "Estimated time remaining " + str(estimated // 60) + " minutes, " + str(estimated % 60) + " seconds."
		s += " " * (linelen - len(s))
		print(s, end="\r", flush=True)

#Input: iteration count i, total number of iterations n, and initial time
#Output: prints progress (using a progress bar) and estimated remaining time
def update_time_2(i, n, init_time, l=304, linelen=80):
	eigth = {1:'\u258F', 2:'\u258E', 3:'\u258D', 4:'\u258C', 5:'\u258B', 6:'\u258A', 7:'\u2589'}
	if i+1 >= n:
		print(" " * (70 + int(l/8)), end="\r", flush=True)
	else:
		elapsed = time.time() - init_time
		estimated = int(elapsed/(i+1) * (n-(i+1)))
		pct = math.floor((l + 0.0)*i/n)
		s = '\u2588' * math.floor(pct/8.0)
		if pct % 8 in eigth:
			s += eigth[pct % 8]
		s += " " * (int(l/8) - len(s)) + '\u258F' + " "
		s += str(estimated // 60) + " minutes, " + str(estimated % 60) + " seconds." + (" " * 5)
		print(s, end="\r", flush=True)

#Input: two dimensional array (matrix)
#Output: transpose matrix
def trans_mat(mat):
	return([[mat[i][j] for i in range(len(mat))] for j in range(len(mat[0]))])

#Input: list a, list b
#Output: list a without any elements of b
def remove_from_list(a,b):
	for bi in b:
		if bi in a:
			a.remove(bi)
	return(a)

#Input: string s, list of strings l:
#Output: if any strings from li in s
def check_stringmatch(s, l, case_dependent=False):
	for li in l:
		if li in s or (not case_dependent and li.lower() in s.lower()):
			return(True)
	return(False)

#Input: hue, saturation, and value for a color
#Output: RGB color system values
def hsv_to_rgb(h, s, v):
	if s == 0.0:
		return v, v, v
	i = int(h*6.0) # XXX assume int() truncates!
	f = (h*6.0) - i
	p = v*(1.0 - s)
	q = v*(1.0 - s*f)
	t = v*(1.0 - s*(1.0-f))
	i = i%6
	if i == 0:
		return v, t, p
	if i == 1:
		return q, v, p
	if i == 2:
		return p, v, t
	if i == 3:
		return p, q, v
	if i == 4:
		return t, p, v
	if i == 5:
		return v, p, q

#Input: iteration i, total number of iterations n, cycling factor r
#Output: RGB color value
#Designed to alternate colors to avoid having adjacent similar colors
def gen_color(i, n, r=2, s=1.0, v=1.0):
	n = n+1
	h = i/(n+0.0) + (i%r)*(1.0/r)
	if h > 1:
		h = h-1
	return(tuple(hsv_to_rgb(h, s, v)))

#Input: hex string
#Output: RGB color values
def hex_to_rgb(s):
	if "#" in s:
		s = s.split("#")[1:]
	color = []
	for i in range(len(s)//2):
		color.append(int(s[2*i:2*i+2],16)/16**2)
	return(color)

#Input: nucleotide sequence
#Output: dictionary of position:[codon, amino acid] for sequence
def aaseq_posdict(nt_seq):
	codons = [nt_seq[3*i:3*(i+1)] for i in range(len(nt_seq)//3)]
	aas = [codon_to_aa[c] for c in codons if codon_to_aa[c] != 'Stop']
	d = {i:[codons[i], aas[i]] for i in range(len(aas))}
	return d

#Input: path to CSV/TSV file
#Output: matrix corresponding to CSV/TSV file
def read_csv(inf, sep="\t"):
	data = []
	with open(inf, "r") as f:
		row = []
		for line in f.readlines():
			data.append((line.strip()).split(sep))
	return(data)

#Input: path to CSV/TSV file
#Optional inputs: index_column, header
#Output: matrix corresponding to CSV/TSV file
def read_csv_2(inf, sep="\t", index_col=None, header=0):
	data = {}
	headers = {}
	with open(inf, "r") as f:
		for i, line in enumerate(f.readlines()):
			line = line.strip()
			if i == header:
				headers = line.split(sep)
				if index_col != None and index_col in headers:
					index_col = headers.index(index_col)
			else:
				if not is_numeric(index_col):
					index = i
				else:
					index = line.split(sep)[index_col]
				data[index] = {headers[i]:line.split(sep)[i] for i in range(len(headers))}
	return(data)

#Input: number, number of significant digits, color to indicate formatting of significant values
#Output: rounded number string, potentially colored to indicate significance
def significance(num, digits=3, color="red"):
	if round(num,digits) < math.pow(10, -digits):
		returnstr = str(math.pow(10, -digits))
	else:
		returnstr = str(round(num,digits))
	if num <= 0.05 and color:
		returnstr = "\033[91m" + returnstr + "\033[00m"
	return(returnstr)

#Input: dictionary of k:count
#Output: dictionary of k:frequency
def convert_to_freqs(d):
	tot = sum(d.values())
	l = len(list(d.keys())[0])
	d = {k:d[k]/tot*(10**l) for k in d}
	return(d)

#Input: string
#Output: string padded with spaces to fill a number of characters
def padstring(s, n=10, space=" "):
	if len(s) < n:
		return(s+(space*(n-len(s))))
	return(s)

#Input: path to file, two-dimensional matrix, labels of rows and columns
#Output writes a CSV/TSV file
def write_to_csv(filename, data, labelsX=[], labelsY = []):

	if labelsX == None or len(labelsX) < 1:
		labelsX = list(range(dim(data)[0]))
	if labelsY == None or len(labelsY) < 1:
		labelsY = list(range(dim(data)[1]))
	if	dim(data)[0] != len(labelsX) and dim(data)[1] != len(labelsY):
		raise Exception("Dimension of labels and data do not match")
	with open(filename, "w") as outf:
		for i in range(-1, len(data)):
			for j in range(-1, len(data[0])):
				if i == -1 and j == -1:
					outf.write("")
				elif i == -1:
					outf.write(labelsY[j])
				elif j == -1:
					outf.write(labelsX[i])
				else:
					outf.write(data[i][j])
				outf.write("\t")
			outf.write("\n") 

#Input: CSV/TSV file with no header and only two columns
#Output: reads CSV/TSV and computes an index:value dictionary
def two_column_csv(infile, sep="\t"):
	data = {}
	with open(infile, "r") as inf:
		for line in inf:
			pieces = line.strip().split(sep)
			if len(pieces) >= 2:
				if is_numeric(pieces[1]):
					data[pieces[0]] = convert_to_num(pieces[1])
				else:
					data[pieces[0]] = pieces[1]
	return(data)

#Input: dataframe, names of x and y columns
#Output: generates a scatterplot of x and y data in dataframe
def scatterplot(df, col1, col2):
	x = list(df[col1])
	y = list(df[col2])
	try:
		plt.scatter(x, y, c="black", alpha=0.1)
		plt.xlabel(col1)
		plt.ylabel(col2)
		plt.savefig("./out/" + col1 + "_" + col2 + ".png")
		plt.cla()
		plt.clf()
	except Exception as e:
		print(col1 + " " + col2)
		print(e)

#Input: Dictionary of sequences
#Output: name and sequence for sequence of median length
def median_seqs(seqs):
	seqs = [(name, seq) for name,seq in seqs.items()]
	lens = [len(pair[1]) for pair in seqs]
	middle = statistics.median(lens)
	index = lens.index(middle)
	return(seqs[index])

#Input: path to codon usage table file
#Output: dictionary of codon:count
def parse_cut(infile):
	codon_data = {}
	with open(infile, "r") as codonfile:
		for line in codonfile:
			codons = re.findall("[TCAGUtcagu]{3}", line)
			numbers = re.findall("\(\s*\d+\.?\d*\s*\)", line)
			for i,codon in enumerate(codons):
				if i < len(numbers):
					number = max(re.findall("\d+\.?\d*", numbers[i]))
					codon_data[codon.upper().replace("U", "T")] = float(number)
				else:
					raise Exception("Incorrect datafile")
	return(codon_data)

#Input: input, function, default value
#Output: output of function on input, or default if errors are encountered
def func_or_default(l, func=math.log, default=0.0, verbose=False):
	try:
		return(func(l))
	except Exception as e:
		if verbose:
			print(e)
		return(default)


#Input: dictionary of name:seq sequences, path to file
#Output: writes a multiple-sequence fasta file
def write_fasta(seqs, outf):
	with open(outf, "w") as outf:
		for name, seq in seqs.items():
			outf.write(">" + name + "\n")
			outf.write(seq + "\n")

#Input: path to file
#Output: dictionary of name:sequence for sequences
def read_fasta(inf, aformat="FIRST", duplicate="replace"):
	data = {}
	with open(inf, "r") as fa:
		name = ""
		for line in fa.readlines():
			if "#" in line:
				continue
			if ">" in line:
				if aformat.upper() == "NCBI":
					name = re.search(">[a-zA-Z]+_?\d+(\.\d+)*", line).group(0)
				elif aformat.upper() in ["FIRST", "WORD"]:
					name = line.split()[0]
				else:
					name = line.strip()
				name = name[1:].strip()
				if name in data.keys():
					if duplicate.lower() in ["append", "a"]: #simply add to existing sequence
						pass
					elif duplicate.lower() in ["replace", "r"]: #reset sequence to empty
						data[name] = ""
					elif duplicate.lower() in ["separate", "s"]: #add underscore+number to end of sequence name
						matches = re.findall("/_\d+$/", name)
						if matches != None and len(matches) > 0:
							num = int(max(matches)[1:])
							name = name[:-len(str(num))] + str(num+1)
							data[name] = ""
						else:
							name = name + "_2"
							data[name] = ""
				else:
					data[name] = ""
			else:
				data[name] = data[name] + line.strip()
	return(data)

#Input: sequences output from BLAST
#Output: sequences combined if they have the same RefSeq/GenBank accession ID
#combines different hits from the same sequence from BLAST
def agglomerate_seqs(seqs):
	newseqs = {seqid.split(":")[0]:{name.split(":")[1].split("-")[0]:seq for name, seq in seqs.items() if name.split(":")[0] == seqid.split(":")[0]} for seqid in seqs}
	newseqs = {seqid:sorted(seqs.items(), key=lambda kv: kv[0]) for seqid, seqs in newseqs.items()}
	newseqs = {seqid:[kv[1] for kv in seqs] for seqid, seqs in newseqs.items()}
	newseqs = {seqid:"".join(seqs) for seqid, seqs in newseqs.items()}
	return(newseqs)

#Input: array of values
#Output: mean of values
def calc_mean(arr):
	return(sum(arr)/len(arr))

#Input: a value: a list of length-2 lists [[a,b], [c,d], ..., [e,f]]
#Output: indicates if i is in any range
def contained(i, ranges):
	for range_vals in ranges:
		if i >= range_vals[0] and i <= range_vals[1]:
			return 1
	return 0

#Input: string
#Output: finds a substring resembling an NCBI accession, with version
def get_accession(key):
	return max(re.findall("\w+\_?\d+\.?\d*", key))

#Input: string
#Output: position and characters changed
#multiple use mutation parser, can parse different formats
def get_mutation_data(line):
	aas = re.findall("[a-zA-Z*]", line)
	pos = int(max(re.findall("\d+", line), key=lambda kv: int(kv)))
	return pos, aas

#Input: mutation, reference sequence, index of mutation (usually 1)
#Output: indicates if reference character in mutation matches sequence
def check_mut(mut, seq, index=1):
	if isinstance(mut, str) and not isinstance(mut, (tuple, list)):
		mut = get_mutation_data(mut)
	if not seq[mut[0]-index] == mut[1][0]:
		raise Exception("Sequence does not match mutant: " + str(mut) + " vs " + seq[mut[0]-index])

#Input: any string or object
#Output: indicates if string is numeric
def is_numeric(s):
	try:
		float(s)
		return True
	except Exception as e:
		return False

#Input: string
#Output: string converted to number if is numeric
def convert_to_num(s):
	if is_numeric(s):
		if float(s) == int(float(s)):
			return(int(float(s)))
		else:
			return(float(s))
	else:
		return(s)
		
#finds the AA variant, given a NT variant and the sequence
#checks the NT sequence matches the WT NT
#checks the AA sequence if given
def get_mutant_aa(ntmut, ntseq, aaseq=None, index=0):
	codonpos = (ntmut[0]-index)//3*3
	wtaa = codon_to_aa[ntseq[codonpos:codonpos+3]]
	mutseq = update_str(ntseq, ntmut[1][1], ntmut[0]-index)
	mutaa = codon_to_aa[mutseq[codonpos:codonpos+3]]
	if ntseq[ntmut[0]-index] != ntmut[1][0]:
		warnings.warn("NTs don't match at position " + str(ntmut[0]) + " for mutant " + str(ntmut))
	if aaseq != None:
		if len(aaseq) > (ntmut[0] - index)//3 and aaseq[(ntmut[0] - index)//3] != wtaa:
			warnings.warn("AAs don't match at position " + str(codonpos) + " for mutant " + str(ntmut))
		elif len(aaseq) <= (ntmut[0] - index)//3:
			print("AA sequence isn't long enough to match NT sequence")
	return((codonpos//3+index, (wtaa, mutaa)))

#Input: value
#Output: indicates if value is valid numeric (False)
#nan numbers are numeric but not useful
def is_nan(val):
	if is_numeric(val) and not np.isnan(val):
		return(False)
	else:
		return(True)

#Input: list of values, integer n
#Output: n-dimensional cartesian product
def perm(l, n, mode="list"):
	if l == [] or l == "" or n == 0:
		if mode in ["string", "str", "s"]:
			return([""])
		else:
			return([[]])
	return_val = []
	sub = perm(l, n-1, mode=mode)
	for i in l:
		for i_sub in sub:
			if mode in ["string", "str", "s"]: 
				return_val.append(i + i_sub)
			else:
				return_val.append([i] + i_sub)
	return(return_val)

#Input: sequences
#Output: dictionary of codon:count
#calculates codon counts for a sequence and adds them to dictionary
def calc_cu_seq(seq, d={}):
	for s in perm(["T", "G", "A", "C"], 3):
		if s not in d.keys():
			d[s] = 0
	for i in range(len(seq)//3):
		codon = seq[3*i: 3*(i+1)]
		d[codon] = d[codon] + 1
	return(d)

#Input: dictionary of codon:count
#Output: W (codon adaptation index)
#computes relative adaptiveness from a codon usage dictionary
def calc_w(d):
	w = {}
	for k in d.keys():
		try:
			w[k] = d[k]/max([d[codon] for codon in codon_table[codon_to_aa[k]]])
		except Exception as e:
			w[k] = 1
	return(w)

#Input: nucleotide sequences
#Output: complement sequence
#gives complement of NT sequence
def complement(seq):
	seq = seq.upper()
	seq = (((((seq.replace("A", "X")).replace("T", "A")).replace("X", "T")).replace("C", "Y")).replace("G", "C")).replace("Y", "G")
	return(seq)

#Input: dictionary of sequences
#Output dictionary of multiple sequence alignment from Clustal Omega
#aligns sequences (input as a dictionary of sequence name: sequence)
def align_sequences(seqs, clustalo_cmd='/usr/bin/clustalo'):
	with subprocess.Popen([clustalo_cmd, '-i', '-'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
		clustalo_input = bytes(''.join(">%s\n%s\n" % (str(key), str(seqs[key])) for key in seqs), 'ascii') 
		clustalo_output, clustalo_err = proc.communicate(clustalo_input)
	aligned_seqs = {}
	name = ""
	for line in clustalo_output.decode('ascii').split('\n'):
		try:
			if len(line) > 0 and line[0] == '>':
				name = line[1:]
				aligned_seqs[name] = ""
			elif len(line) > 0:
				aligned_seqs[name] += line.strip()
		except Exception as e:
			print(e)
			print(line.strip())
	return(aligned_seqs)

#Input: sequences
#Output: checksum of sequence (for GCG)
def checksumCGC(seq):
	seq = seq.upper()
	check = 0
	for i in range(len(seq)):
		if re.match("[A-Z]", seq[i]):
			check += ((i % 57) + 1) * ord(seq[i])
	return(check % 10000)

#Input: path to input file in FASTA format
#Output: writes GCG file equivalent to FASTA file
#converts a FASTA file to GCG format
def fasta_to_gcg(infile, alphabet="-ACDEFGHIKLMNPQRSTVWY"):
	seqs = read_fasta(infile)
	with open(os.path.splitext(infile)[0]+".gcg", "w") as outf:
		for name, seq in seqs.items():
			outf.write(name.strip() + "\tLength: " + str(len(seq)) + "\t Check: " + str(checksumCGC(seq)) + " ..\n\n")
			for i in range(len(seq)):
				if i % 50 == 0:
					outf.write(" " + ('{0: ' + str(round(math.log(len(seq), 10), None)+1) + 'd}').format(i+1) + "\t")
				elif i % 50 == 49:
					outf.write(seq[i] + "\n")
				elif i % 10 == 0:
					outf.write(" " + seq[i])
				else:
					outf.write(seq[i])
			outf.write("\n\n")
#input: dictionary of sequences, additional parameters to give
#Output: aligned sequences from Clustal Omega
#aligns sequences with opportunities for additional parameters (input as a dictionary of sequence name: sequence)
def align_sequences_params(seqs, mode="clustalo", clustalo_cmd='/usr/bin/clustalo', params="", tmpdir="/media/temp/"):
	write_fasta(seqs, tmpdir + "in.fasta")
	with subprocess.Popen([clustalo_cmd] + params.split() + ["--outfmt", "fa", '-i', tmpdir+"in.fasta", "-o", tmpdir+"out.fasta"], \
						  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE) \
						  as proc:
		clustalo_output, clustalo_err = proc.communicate()
	return(read_fasta(tmpdir + "out.fasta"))

#Input: name of codon pair (ABCDEF), dictionary of codon pair:count
#Output: codon pair score for input codon pair
#computed codon pair score
def cps(bicodon, d):
	a = b = x = y = 0.0
	ab = xy = 0.0
	for i in range(len(d)):
		if d[i][1] == bicodon[0:3]:
			a = a + 1
		if d[i][1] == bicodon[4:6]:
			b = b + 1
		if d[i][0] == codon_to_aa[bicodon[0:3]]:
			x = x + 1
		if d[i][0] == codon_to_aa[bicodon[4:6]]:
			y = y + 1
		if i+1 in d.keys():
			if bicodon == d[i][1] + d[i+1][1]:
				ab = ab + 1
			if codon_to_aa[bicodon[0:3]]+codon_to_aa[bicodon[4:6]] == d[i][0]+d[i+1][0]:
				xy = xy + 1
	try:
		return(ab/(a*b/(x*y)*xy))
	except Exception as e:
		return(e)

#Input: a numeric list of dictionary of key:value
#Output: normalized data (relative to empirical mean and standard deviation)
#normalizes a list of numbers to z-scores
def normalize(data):
	if isinstance(data, (dict,)):
		u = np.mean(data.values())
		o = np.std(data.values())
		for key in data.keys():
			data[key] = (data[key] - u)/o
	else:
		u = np.mean(data)
		o = np.std(data)
		for i in range(len(data)):
			data[i] = (data[i] - u)/o
	return(data)

#Input: value and a list of values
#Output: empirical p-value (2*fraction of values in l more extreme than val)
#will take whichever side of l would result in a lower p-value
def pval_dist(val, l):
	l = [float(xi) for xi in l if is_numeric(xi) and not np.isnan(float(xi))]
	return(2/len(l)*min([len([xi for xi in l if xi < val]), len([xi for xi in l if xi > val])]))

#Input: a list of object
#Output: dimensions of object
#returns dimension of a multi-dimensional data structure
def dim(l):
	if isinstance(l, (str,)):
		return([])
	try: 
		length = len(l)
	except:
		return([])
	if length > 0:
		return([length] + dim(list(l)[0]))
	else:
		return([length])

#Input: two arrays of same length
#Output: same arrays, but only positions that are numeric and not NaN in both
def clean_arrays(arr1, arr2):
	good = [i for i in range(len(arr1)) if is_numeric(arr1[i]) and not np.isnan(float(arr1[i])) and is_numeric(arr2[i]) and not np.isnan(float(arr2[i]))]
	return([float(arr1[i]) for i in good], [float(arr2[i]) for i in good])

#Input: two arrays of same length
#Output: mean absolute percent error of two arrays
def calc_mape(arr1, arr2):
	if len(arr1) != len(arr2):
		raise Exception("Arrays have different lengths.")
	return(100*np.mean([abs(arr1[i]-arr2[i])/arr1[i] for i in range(len(arr1)) if arr1[i] != 0]))

#Input: two arrays of same length
#Output: mean absolute square error of two arrays
def calc_mse(arr1, arr2):
	if len(arr1) != len(arr2):
		raise Exception("Arrays have different lengths.")
	return(sum([(arr1[i]-arr2[i])**2 for i in range(len(arr1))])/len(arr1))

#Input: string (usually parsed from website)
#Output: decoded string if bytes, else original string
#Note: this is useful when receiving data from internet, as encoding may raise issues
def clean_web_bytes(s):
	if isinstance(s, bytes):
		return(s.decode("utf-8"))
	else:
		return(s)
