import os, sys, re, subprocess, math, time, statistics, pickle, requests, distance
import random
import numpy as np
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.cluster

def complement(seq):
	d = {"C":"G", "G":"C", "A":"T", "T":"A"}
	if all([True if seq[i] in d.keys() else False for i in range(len(seq))]):
		return("".join([d[seq[i]] for i in range(len(seq))]))

def pctile(l, k):
	l = [float(li) for li in list(l) if is_numeric(li) and not np.isnan(float(li))]
	if len(l) >= 0 and is_numeric(k) and not np.isnan(float(k)):
		return(len([li for li in l if li <= float(k)])/len(l))

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


#computes levenshtein distance between strings
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

#retry a given function (useful for networking functions)
def retry_func(func, positional_arguments=[], keyword_arguments={}, lim=10, wait=2):
	for i in range(lim):
		try:
			return(func(*positional_arguments, **keyword_arguments))
		except Exception as e:
			time.sleep(wait*i)
	warnings.warn("\033[93m" + str(func.__name__) + " failed after " + str(lim) + " tries." + "\033[0m") 
	return(None)

#find the subsequence centered on pos with the l flanking characters (l/2 on each side)
def subseq(seq, pos, l):
	return(seq[max(0, pos - math.floor(l/2)):min(len(seq), pos + math.ceil(l/2))])

#computes a sequence uniqueness score
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

#computes similarity between sequences using a custom function (usually levenshtein, but can be changed)
def seqs_compute_similarity(seqs, genename=None, space="-", func=lambda k: levenshtein(k[0], k[1])):
	if isinstance(seqs, (dict,)):
		names = list(seqs.keys())
	elif isinstance(seqs, (list, tuple)):
		names = [i for i in range(len(seqs))]
	seqs = [seqs[name] for name in names]
	similarity = -1*np.array([[func((seq1.replace(space, ""), seq2.replace(space, ""))) for seq1 in seqs] for seq2 in seqs])
	return(similarity)

#cluster sequences by by precomputed similarity matrix (uses Affinity Propagation)
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

#find the cluster centers with a similarity matrix
def cluster_centers(similarity, cluster, names=[]):
	if names == None or len(names) < 1:
		names = [i for i in range(len(similariy))]
	cluster.fit(-np.array(similarity))
	labels = cluster.predict(-np.array(similarity))
	scores = cluster.transform(-np.array(similarity))
	distances = {i:statistics.mean([(1-similarity[i][j])**2 for j in range(len(labels)) if labels[i] == labels[j]]) for i in range(len(labels))}
	clusters = {names[min([i for i in range(len(labels)) if labels[i] == labels[j]], key=lambda kv: distances[kv])]:[names[i] for i in range(len(labels)) if labels[i] == labels[j]] for j in range(len(labels))}
	return(clusters)

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

#converts a mutation coordinates (uses genomic from GRCh37 by default)
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

#check a variant with mutalyzer
def check_variant(refid, var, c="c."):
	var = get_mutation_data(var)
	address = "https://mutalyzer.nl/name-checker?description=" + str(refid) + "%3A" + str(c) + str(var[0]) + str(var[1][0]) + "%3E" + str(var[1][1])
	r = retry_request(requests.get, [address])
	errors = re.findall("(\d+) Error", r.text)
	warnings = re.findall("(\d+) Warnings", r.text)
	return([int(x) for x in errors], [int(x) for x in warnings])

#retries a requests library function
#assumes a response object is returned (which can be raised for HTML status)
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

#infers the alphabet of a sequence/string (either NT or amino acid)
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

codons = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"]

#get the (longest) number from a string
def get_num(s):
	num = float(max(re.findall("\d+\.?\d*", s), key=len))
	if num - int(num) == 0.0:
		return(int(num))
	else:
		return(num)

#useful when sending data online which must be formatted
def slugify(value):
	value = str(re.sub('[^\w\s-]', '', value).strip().lower())
	value = str(re.sub('[-\s]+', '-', value))
	return(value)

#load from pickle (binary data) file
def pickleload(f):
	with open(f, "rb") as pklfile:
		d = pickle.load(pklfile)
	return(d)

#save to pickle (binary data) file
def pickledump(d, f):
	with open(f, "wb") as pklfile:
		pickle.dump(d,pklfile)

#creates a mutant string
def update_str(s, c, pos):
	return(s[:pos] + c + s[pos+1:])

#gives updates, including estimated remaining time, when used in a for loop
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

#also updates on time, but uses a wait bar for time
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

def trans_mat(mat):
	return([[mat[i][j] for i in range(len(mat))] for j in range(len(mat[0]))])

def remove_from_list(a,b):
	for bi in b:
		if bi in a:
			a.remove(bi)
	return(a)

def check_stringmatch(s, l, case_dependent=False):
	for li in l:
		if li in s or (not case_dependent and li.lower() in s.lower()):
			return(True)
	return(False)

#converts colors from HSV to RGB system. Useful when automating color choice
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

#used for plotting to generate rainbow colors, but alternate colors by hue
def gen_color(i, n, r=2, s=1.0, v=1.0):
	n = n+1
	h = i/(n+0.0) + (i%r)*(1.0/r)
	if h > 1:
		h = h-1
	return(tuple(hsv_to_rgb(h, s, v)))

#Converts hex strings to RGB system
def hex_to_rgb(s):
	if "#" in s:
		s = s.split("#")[1:]
	color = []
	for i in range(len(s)//2):
		color.append(int(s[2*i:2*i+2],16)/16**2)
	return(color)

#produces a dictionary of position:[codon, amino acid] entries for a nucleotide sequence
def aaseq_posdict(nt_seq):
	codons = [nt_seq[3*i:3*(i+1)] for i in range(len(nt_seq)//3)]
	aas = [codon_to_aa[c] for c in codons if codon_to_aa[c] != 'Stop']
	d = {i:[codons[i], aas[i]] for i in range(len(aas))}
	return d

#parses a csv file. Use pd.read_csv if possible
def read_csv(inf, sep="\t"):
	data = []
	with open(inf, "r") as f:
		row = []
		for line in f.readlines():
			data.append((line.strip()).split(sep))
	return(data)

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

#creates appropriate string indicating significance and rounds to 3 digits
def significance(num, digits=3, color="red"):
	if round(num,digits) < math.pow(10, -digits):
		returnstr = str(math.pow(10, -digits))
	else:
		returnstr = str(round(num,digits))
	if num <= 0.05 and color:
		returnstr = "\033[91m" + returnstr + "\033[00m"
	return(returnstr)

#using a dictionary of key:number/count values, converts all to frequencies
def convert_to_freqs(d):
	tot = sum(d.values())
	l = len(list(d.keys())[0])
	d = {k:d[k]/tot*(10**l) for k in d}
	return(d)

def padstring(s, n=10, space=" "):
	if len(s) < n:
		return(s+(space*(n-len(s))))
	return(s)

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

#reads a two column csv into a data dictionary
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

#create a scatterplot using data from df of col1 vs col2
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

def median_seqs(seqs):
	seqs = [(name, seq) for name,seq in seqs.items()]
	lens = [len(pair[1]) for pair in seqs]
	middle = statistics.median(lens)
	index = lens.index(middle)
	return(seqs[index])

#read a codon usage table (CUT) file
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

#returns either a function on a value, or some default value if error
def func_or_default(l, func=math.log, default=0.0, verbose=False):
	try:
		return(func(l))
	except Exception as e:
		if verbose:
			print(e)
		return(default)


#writes a fasta file
def write_fasta(seqs, outf):
	with open(outf, "w") as outf:
		for name, seq in seqs.items():
			outf.write(">" + name + "\n")
			outf.write(seq + "\n")

#parses a fasta file
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

#function to combine different hits from the same sequence from BLAST
def agglomerate_seqs(seqs):
	newseqs = {seqid.split(":")[0]:{name.split(":")[1].split("-")[0]:seq for name, seq in seqs.items() if name.split(":")[0] == seqid.split(":")[0]} for seqid in seqs}
	newseqs = {seqid:sorted(seqs.items(), key=lambda kv: kv[0]) for seqid, seqs in newseqs.items()}
	newseqs = {seqid:[kv[1] for kv in seqs] for seqid, seqs in newseqs.items()}
	newseqs = {seqid:"".join(seqs) for seqid, seqs in newseqs.items()}
	return(newseqs)

def calc_mean(arr):
	return(sum(arr)/len(arr))

#determines whether i is in any of the ranges
def contained(i, ranges):
	for range_vals in ranges:
		if i >= range_vals[0] and i <= range_vals[1]:
			return 1
	return 0

#read only the accession key from a string
def get_accession(key):
	return max(re.findall("\w+\_?\d+\.?\d*", key))

#reads the mutation data from a string
def get_mutation_data(line):
	aas = re.findall("[a-zA-Z*]", line)
	pos = int(max(re.findall("\d+", line), key=lambda kv: int(kv)))
	return pos, aas

def check_mut(mut, seq, index=1):
	if isinstance(mut, str) and not isinstance(mut, (tuple, list)):
		mut = get_mutation_data(mut)
	if not seq[mut[0]-index] == mut[1][0]:
		raise Exception("Sequence does not match mutant: " + str(mut) + " vs " + seq[mut[0]-index])

#checks if a string (or object) can be represented as a number
def is_numeric(s):
	try:
		float(s)
		return True
	except Exception as e:
		return False

#converts a string (or variable) to a number if it's numeric
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

def is_nan(val):
	if is_numeric(val) and not np.isnan(val):
		return(False)
	else:
		return(True)

#reads a sheet of data containing ranges (domains or exons)
def read_ranges(wb, ws_name):
	ws = wb[ws_name]

	regions = {}
	for row in ws.iter_rows(min_row=2, min_col=1):
		name = ""
		val = []
		for cell in row:
			if cell.column is 'A': #labels column
				name = str(cell.value)
				i = 1
				while name in regions.keys(): #forces unique region names (dictionary keys)
					i = i+1
					name = cell.value + " " + str(i)
			else: #values column
				val.append(int(cell.value))
		regions[name] = val
	return regions

#computes n-dimensional cartesian product of a string or list
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

#calculates codon counts for a sequence and adds them to dictionary
def calc_cu_seq(seq, d={}):
	for s in perm(["T", "G", "A", "C"], 3):
		if s not in d.keys():
			d[s] = 0
	for i in range(len(seq)//3):
		codon = seq[3*i: 3*(i+1)]
		d[codon] = d[codon] + 1
	return(d)

#computes the codon usage dictionary for an input file
def calc_cut(inf):
	d = {}
	seq = ""
	cdinfo = []
	numseqs = 0
	numcodons = 0
	with open(inf, "r") as inf:
		for line in inf:
			if ">" in line or "#" in line:
				for cds in cdinfo:
					try:
						#look for coding sequence info in header
						cd_range = re.findall("\d+", cds)
						if len(cd_range) >= 2:
							seq = seq[int(cd_range[0])-1:int(cd_range[1])-1]
						else:
							pass
						#if sequences have been read, add their counts to dictionary
						if seq != "":
							d = calc_cu_seq(seq, d)
							numcodons = numcodons + (len(seq)//3)
					except Exception as e: #skip sequence
						print(str(e))
				#read new info
				cdinfo = re.findall("CDS:\d+-\d+", line)
				seq = ""
				numseqs = numseqs + 1
			else:
				seq = seq + line.strip()
	print(str(numseqs) + " sequences read")
	print(str(numcodons) + " codons read")
	return(d)

#computes relative adaptiveness from a codon usage dictionary
def calc_w(d):
	w = {}
	for k in d.keys():
		try:
			w[k] = d[k]/max([d[codon] for codon in codon_table[codon_to_aa[k]]])
		except Exception as e:
			w[k] = 1
	return(w)

#gives complement of NT sequence
def complement(seq):
	seq = seq.upper()
	seq = (((((seq.replace("A", "X")).replace("T", "A")).replace("X", "T")).replace("C", "Y")).replace("G", "C")).replace("Y", "G")
	return(seq)

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

def checksumCGC(seq):
	seq = seq.upper()
	check = 0
	for i in range(len(seq)):
		if re.match("[A-Z]", seq[i]):
			check += ((i % 57) + 1) * ord(seq[i])
	return(check % 10000)

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

#aligns sequences with opportunities for additional parameters (input as a dictionary of sequence name: sequence)
def align_sequences_params(seqs, mode="clustalo", clustalo_cmd='/usr/bin/clustalo', params="", tmpdir="/media/temp/"):
	write_fasta(seqs, tmpdir + "in.fasta")
	with subprocess.Popen([clustalo_cmd] + params.split() + ["--outfmt", "fa", '-i', tmpdir+"in.fasta", "-o", tmpdir+"out.fasta"], \
						  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE) \
						  as proc:
		clustalo_output, clustalo_err = proc.communicate()
	return(read_fasta(tmpdir + "out.fasta"))

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

def cps_positions(x, y, d):
	if x in d.keys() and y in d.keys():
		cps = cps(d[x][1]+d[y][1], d)
		if is_numeric(cps):
			return(float(cps))
	return(0)

def opt_cpb(seq, mode="max"):
	d = aaseq_posdict(seq)
	cpb = 0
	count = 0
	for i in range(len(d)-1):
		cps = cps(d[i][1]+d[i+1][1], d)
		if is_numeric(cps):
			cpb = cpb + cps
			count = count + 1
	cpb = cpb/count

	for i in range(len(d)**2):
		rand_pos = random.randint(0, len(d))
		x,y = random.sample([i for i in range(len(d)) if d[i][0] == d[rand_pos][0]],2)
		orig = 0
		changed = 0
		for p in [x, x-1, y, y-1]:
			orig = orig + cps_positions(p, p+1,  d)
		
		for pair in [[x, y+1], [x-1, y], [y, x+1], [y-1, x]]:
			changed = changed + cps_positions(pair[0], pair[1], d)

		rand = random.random()
		if (changed > orig and mode == "max") or (changed < orig and mode == "min") or rand >= 0.5 + i/2.0/((len(d))**2):	
			tmp = d[x]
			d[x] = d[y]
			d[y] = tmp
			cpb = cpb - orig/count + changed/count
	return(cpb, d)

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

def pval_dist(val, l):
	l = [float(xi) for xi in l if is_numeric(xi) and not np.isnan(float(xi))]
	return(2/len(l)*min([len([xi for xi in l if xi < val]), len([xi for xi in l if xi > val])]))

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

def clean_arrays(arr1, arr2):
	good = [i for i in range(len(arr1)) if is_numeric(arr1[i]) and not np.isnan(float(arr1[i])) and is_numeric(arr2[i]) and not np.isnan(float(arr2[i]))]
	return([float(arr1[i]) for i in good], [float(arr2[i]) for i in good])

def calc_mape(arr1, arr2):
	if len(arr1) != len(arr2):
		raise Exception("Arrays have different lengths.")
	return(100*np.mean([abs(arr1[i]-arr2[i])/arr1[i] for i in range(len(arr1)) if arr1[i] != 0]))

def calc_mse(arr1, arr2):
	if len(arr1) != len(arr2):
		raise Exception("Arrays have different lengths.")
	return(sum([(arr1[i]-arr2[i])**2 for i in range(len(arr1))])/len(arr1))

def clean_web_bytes(s):
	if isinstance(s, bytes):
		return(s.decode("utf-8"))
	else:
		return(s)

def combine_fasta(fasta):
	seqs = read_fasta(fasta)
	seqs = {k.split(":")[0]:{int(k2.split(":")[1].split("-")[0]):v2 for k2,v2 in seqs.items() if k2.split(":")[0] == k.split(":")[0]} for k in seqs}
	seqs = {k:sorted(v.items(), key=lambda kv: kv[0]) for k,v in seqs.items()}
	seqs = {k: "".join([vi[1] for vi in v]) for k,v in seqs.items()}
	return(seqs)
