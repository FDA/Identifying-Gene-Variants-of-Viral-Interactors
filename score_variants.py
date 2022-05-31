import os, re, sys, time, os.path, argparse, math, statistics, pickle, subprocess, traceback, datetime, requests, combine_documents, warnings, openpyxl, bio_tools, copy, tools, stat_tools
from collections import Counter
import pandas as pd
import numpy as np
import scipy.stats, scipy.optimize
import compute_features
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression, LinearRegression
import statsmodels.api as sm
from multiprocessing import Process
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['font.family'] = "arial"

ram_disk="insert here"

#convert clinical significance to 0-1 scale
sigdict = {"uncertain-significance":0.5, "likely-benign":0.6, "benign":0.75, "pathogenic":0.25, "likely-pathogenic":0.4, 'not-provided':0.50, "protective":1.0}
def sigdict(k):
	if k in sigdict:
		return(sigdict[k])
	else:
		print("Warning: " + str(k) + " not recognized")
		return(0.5)

markers = {'o': 'circle', 'v': 'triangle_down', '^': 'triangle_up', 's': 'square', '*': 'star'}

geneset = [('F8', 'F8', './LOVD/00495_F8.tsv'), ('F9', 'F9', './LOVD/02237_F9.tsv'), ('DMD', 'DMD', './LOVD/05324_DMD.tsv'), ('G6PD', 'G6PD', './LOVD/00308_G6PD.tsv'), ('HBB', 'HBB', './LOVD/03506_HBB.tsv'), ('VWF', 'VWF/VWD1', './LOVD/00498_VWF.tsv'), ('VWF', 'VWF/VWD2', './LOVD/03361_VWF.tsv'), ('VWF', 'VWF/VWD3', './LOVD/02126_VWF.tsv'), ('F9', 'F9', './F8_F9/F9_dataset_updated.tsv'), ('VWF', 'VWF', './VWD/VWF_dataset_updated.tsv'), ('ADAMTS13', 'ADAMTS13', './cTTP/ADAMTS13_dataset_updated.tsv'), ("ADAMTS13", "ADAMTS13", "./in_vitro_ADAMTS13/Plaimauer_variants.tsv")]

datasets = ["./LOVD/F8/", "./LOVD/F9/", "./LOVD/DMD/", "./LOVD/HBB/", "./LOVD/VWF/VWD1/", "./LOVD/VWF/VWD2/", "./LOVD/VWF/VWD3/", "./LOVD/G6PD/", "./F8_F9/F8/", "./F8_F9/F9/", "./cTTP/ADAMTS13/", "./VWD/VWF/"]

invert = {"Polyphen": True, "SIFT": True, "PROVEAN": True, "PROVEAN (NT)":True, "GnomAD": True, "Conservation (with PDB)": True, "AA percent identity": True, "AA entropy": False, "AA variance": False, "BLOSUM conservation": False, "Rare codon enrichment": True, "Consurf (wo 3D info)": True, "Accessible surface area (PDBPisa)": False, "Solvation energy (PDBPisa)": False, "NT percent identity": True, "NT entropy": False, "NT variance": False, "Relative surface area (NetsurfP2)": False, "MSA loglikelihood (AA)": True, "MSA loglikelihood (NT)": True, "EVmutation (AA)": True, "EVmutation (NT)": True, "Number of mutants": False}

#finds distribution of characters at each position of alignment
def compute_distribution(seqs, genename, space="-", alphabet="", freq=True):
	joined_seqs = "".join(seqs.values()).replace(space, "")
	alphabet = sorted(set([alphabet[i] for i in range(len(alphabet))] + [joined_seqs[i] for i in range(len(joined_seqs))]))
	alphabetsize = len(alphabet)
	#print(seqs)
	inverted_seqs = {i:"".join([seqs[name][i] for name in seqs]) for i in range(len(seqs[genename]))}
	
	location = 1
	distribution = {}
	for i in range(len(seqs[genename])):
		if seqs[genename][i] != space:
			distribution[location] = {alpha:0 for alpha in alphabet}
			for j in range(len(inverted_seqs[i])):
				if inverted_seqs[i][j] != space:
					distribution[location][inverted_seqs[i][j]] += 1
			total = sum(distribution[location].values())
			if freq:
				distribution[location] = {k:v/float(total) for k,v in distribution[location].items()}
			location += 1
	return(distribution)

#generates nucleotides in an MSA at the position of variants
def get_positions(seqs, genename, variants, space="-", index=0):
	varset = sorted(set([(pos, (aas[0], aas[1])) for v in variants.values() for pos, aas in v.items()]))
	joined_seqs = "".join(seqs.values()).replace(space, "")
	alphabet = sorted(set([joined_seqs[i] for i in range(len(joined_seqs))]))
	alphabetsize = len(alphabet)

	inverted_seqs = {i:"".join([seqs[name][i] for name in seqs]) for i in range(len(seqs[genename]))}
	
	distribution = {}
	for (pos, aas) in varset:
		pos2 = compute_features.convert_position2(seqs[genename], pos-index)
		if seqs[genename][pos2] != aas[0]:
			print("Warning: Positions don't match for " + aas[0] + str(pos) + aas[1] +": " + str(pos2))
		if isinstance(pos, (int,)):
			distribution[pos] = inverted_seqs[pos2]
		else:
			distribution[pos] = None
	return(distribution, varset)

#finds distribution of nucleotides at each position of alignment
def variant_distribution(seqs, genename, variants, space="-", index=0):
	distribution, varset = get_positions(seqs, genename, variants, space=space, index=index)
	scores = {}
	for var, muts in variants.items():
		valid = 0.0
		correct = 0.0
		for i in range(len(seqs)):
			validflag = True
			correctflag = True
			for pos, aas in muts.items():
				if distribution[pos][i] != aas[1]:
					correctflag = False
				if distribution[pos][i] == space:
					validflag = False
			if validflag:
				valid += 1
			if correctflag:
				correct += 1
		try:
			scores[var] = correct/valid
		except ZeroDivisionError as e:
			scores[var] = np.nan
	return(scores)

def cond_pair_variants(seqs, genename, variants, space="-", index=0):
	distribution, varset = get_positions(seqs, genename, variants, space=space, index=index)
	scores = {}
	for (pos1, aas1) in varset:
		for (pos2, aas2) in varset:
			name = (aas1[0] + str(pos1) + aas1[1], aas2[0] + str(pos2) + aas2[1])
			valid = 0.0
			correct = 0.0
			for i in range(len(seqs)):
				if distribution[pos1][i] == aas1[1] and distribution[pos2][i] == aas2[1]:
					correct += 1
				if distribution[pos1][i] == aas1[1] and distribution[pos2][i] != space:
					valid += 1
			if correct > valid:
				print("Warning: An error occured while processing pair " + str(name))
			try:
				score = correct/valid
			except ZeroDivisionError as e:
				score = np.nan
			if name[0] not in scores:
				scores[name[0]] = {name[1]:score}
			else:
				scores[name[0]][name[1]] = score
	scores = sorted(scores.items())
	scores = {k:sorted(v.items()) for (k,v) in scores}
	scores = {k:{k1:v1 for (k1,v1) in v} for k,v in scores.items()}
	return(scores)

def varset_vs_ind(seqs, genename, variants, space="-", index=0):
	distribution, varset = get_positions(seqs, genename, variants, space=space, index=index)
	scores = {}
	for var, muts in variants.items():
		valid = 0.0
		correct = 0.0
		ind = 1.0
		for i in range(len(seqs)):
			validflag = True
			correctflag = True
			for pos, aas in muts.items():
				if distribution[pos][i] != aas[1]:
					correctflag = False
				if distribution[pos][i] == space:
					validflag = False
			if validflag:
				valid += 1
			if correctflag:
				correct += 1
		for pos, aas in muts.items():
			ind *= sum([1 for aa in distribution[pos] if aa == aas[1]])/sum([1 for aa in distribution[pos] if aa != space])
		try:
			scores[var] = correct/valid/ind
		except ZeroDivisionError as e:
			scores[var] = np.nan
	return(scores)

def pairset_vs_ind(seqs, genename, variants, space="-", index=0):
	distribution, varset = get_positions(seqs, genename, variants, space=space, index=index)
	scores = {}
	for i, (pos1, aas1) in enumerate(varset):
		for j in range(i):
			(pos2, aas2) = varset[j]
			try:
				valid = len([1 for seq in seqs.values() if seq[pos1-index] != space and seq[pos2-index] != space])
				correct = len([1 for seq in seqs.values() if seq[pos1-index] == aas1[1] and seq[pos2-index] == aas2[1]])
				ind = len([1 for seq in seqs.values() if seq[pos1-index] == aas1[1]])*len([1 for seq in seqs.values() if seq[pos2-index] == aas2[1]])/len([1 for seq in seqs.values() if seq[pos1-index] != space])/len([1 for seq in seqs.values() if seq[pos2-index] != space])
				scores[(aas2[0] + str(pos2) + aas2[1], aas1[0] + str(pos1) + aas1[1])] = correct/valid - ind
			except ZeroDivisionError as e:
				scores[(aas2[0] + str(pos2) + aas2[1], aas1[0] + str(pos1) + aas1[1])] = 0
	return(scores)

def pairset_vs_min(seqs, genename, variants, space="-", index=0):
	distribution, varset = get_positions(seqs, genename, variants, space=space, index=index)
	scores = {}
	for i, (pos1, aas1) in enumerate(varset):
		for j in range(i):
			(pos2, aas2) = varset[j]
			try:
				table = []
				for k in range(2):
					table.append([])
					for l in range(2):
						table[-1].append(len([1 for seq in seqs.values() if seq[pos1-index] == aas1[k] and seq[pos2-index] == aas2[l]]))
				chisq = scipy.stats.chisquare(table, axis=None)[0]
				sign = np.sign(table[0][0]*table[1][1] - table[0][1]*table[1][0])
				pair = sign * math.sqrt(chisq/(table[0][0] + table[0][1] + table[1][0] + table[1][1]))
				if np.isnan(pair):
					raise ZeroDivisionError("Chi-square test produced an error")
				scores[(aas2[0] + str(pos2) + aas2[1], aas1[0] + str(pos1) + aas1[1])] = pair
				print(str((pos1, aas1)) + "\t" + str((pos2, aas2)))
				print(scores[(aas2[0] + str(pos2) + aas2[1], aas1[0] + str(pos1) + aas1[1])])
			except ZeroDivisionError as e:
				scores[(aas2[0] + str(pos2) + aas2[1], aas1[0] + str(pos1) + aas1[1])] = 0
	return(scores)

def pairing_variants(seqs, genename, variants, space="-", index=0):
	distribution, varset = get_positions(seqs, genename, variants, space=space, index=index)
	varset = sorted(varset)
	scores = {}
	for (pos1, aas1) in varset:
		for (pos2, aas2) in varset:
			name = (aas1[0] + str(pos1) + aas1[1], aas2[0] + str(pos2) + aas2[1])
			scores[name] = pair_positions(distribution[pos1], distribution[pos2], space=space)
	newscores = {}
	for k, v in scores.items():
		if k[0] not in newscores:
			newscores[k[0]] = {k[1]:v}
		else:
			newscores[k[0]][k[1]] = v
	scores = {}
	varlist = list(newscores.keys())
	for i, k in enumerate(varlist):
		for j in range(i):
			scores[(varlist[i], varlist[j])] = math.sqrt(abs(newscores[varlist[i]][varlist[j]] * newscores[varlist[j]][varlist[i]]))
	return(scores)
	
def pair_positions(dist1, dist2, space="-"):
	score = 1
	alphabet = sorted(set(dist1).difference(set([space])))
	used = alphabet
	for alpha in sorted(alphabet, key=lambda alpha: dist1.count(alpha), reverse=True):
		subset = [dist2[i] for i in range(len(dist1)) if dist1[i] == alpha and dist2[i] != space and dist2[i] not in used]
		try:
			corr = max(set(subset), key=subset.count)
			score *= subset.count(corr)/len([1 for i in range(len(dist1)) if dist1[i] == alpha])
			used.append(corr)
		except ValueError as e:
			pass
	return(score)

def pairing(seqs, genename, space="-", index=0):
	seq = seqs[genename].replace(space, "")
	variants = {"All":{i+1:[seq[i], seq[i]] for i in range(len(seq))}}
	scores = pairing_variants(seqs, genename, variants, space=space, index=index)
	scores = {(k[1][1:-1], k[0][1:-1]):v for k,v in scores.items()}
	newscores = {}
	for (k1, k2), v in scores.items():
		try:
			newscores[(int(k1), int(k2))] = v
		except Exception as e:
			pass
	return(newscores)

def seq_pos_chisq(seqs, genename, space="-"):
	length = len(seqs[genename].replace(space, ""))
	couplings = {}
	dist = {}
	m = 0
	init_time = time.time()
	for i in range(len(seqs[genename])):
		if seqs[genename][i] != space:
			m += 1
			dist[m] = {}
			n = 0
			for j in range(i):
				if seqs[genename][j] != space:
					n += 1
					observed = {}
					for k in seqs:
						if seqs[k][i] not in dist[m]:
							dist[m][seqs[k][i]] = 1
						else:
							dist[m][seqs[k][i]] += 1
						if seqs[k][j] not in dist[n]:
							dist[n][seqs[k][j]] = 1
						else:
							dist[n][seqs[k][j]] += 1
						if seqs[k][i] != space and seqs[k][j] != space:
							if (seqs[k][i], seqs[k][j]) not in observed:
								observed[(seqs[k][i], seqs[k][j])] = 1
							else:
								observed[(seqs[k][i], seqs[k][j])] += 1
					observed = {k:v for k,v in observed.items() if k[0] != space and k[1] != space}
					expected = {k:dist[m][k[0]] * dist[n][k[1]] for k in observed.keys()}
					scale = sum(observed.values())/sum(expected.values())
					expected = {k:v*scale for k,v in expected.items()}
					couplings[(m,n)] = 1-scipy.stats.chisquare([x for x in observed.values()], [expected[x] for x in observed.keys()]).pvalue
			tools.update_time(sum([z for z in range(m)]), length * (length-1)/2, init_time)
	return(couplings)

#reads an EC file (evolutionary coupling output from PLMC)
def read_ec(df, ORF="", space="-"):
	data = {}
	seq = {}
	with open(df, "r") as inf:
		for line in inf:
			pieces = line.split()
			if pieces[1] != space and pieces[3] != space:
				seq[int(pieces[0])] = pieces[1]
				seq[int(pieces[2])] = pieces[3]
	sorted_seq = sorted(seq.items(), key = lambda kv: kv[0])
	if "".join([s[1] for s in sorted_seq]) != ORF:
		print("Given ssequence does not match input")
		print(ORF)
		print("".join([s[1] for s in sorted_seq]))
	sorted_seq = {k:(i+1,v) for i, (k,v) in enumerate(sorted_seq)}

	with open(df, "r") as inf:
		with open(os.path.splitext(df)[0] + "_aligned.ec", "w") as outf:
			for line in inf:
				pieces = line.split()
				pos1 = int(pieces[0])
				pos2 = int(pieces[2])
				if pos1 in sorted_seq and pos2 in sorted_seq:
					#if pieces[1] != ORF[sorted_seq[pos1][0]-1]:
						#print("Sequences don't match at position " + str(pos1))
					#if pieces[3] != ORF[sorted_seq[pos2][0]-1]:
						#print("Sequences don't match at position " + str(pos2))
					data[(sorted_seq[pos1][0], sorted_seq[pos2][0])] = float(pieces[-1])
					outf.write(str(sorted_seq[pos1][0]) + "\t" + str(sorted_seq[pos2][0]) + "\t" + pieces[-1] + "\n")
	return(sorted(data.items(), key=lambda kv: abs(kv[1])))

#given precomputed (domain) EVmutation scores from database, compiles for entire gene
#finds and requires Uniprot accession
def get_coupling(genename, ntseq, accession=None, directory="./couplings/", cdf=True, absval=True):
	if not os.path.isdir(directory):
		os.mkdir(directory)
	if accession == None or len(accession) < 1:
		uniprot = compute_features.parse_uniprot(genename, ntseq, taxid=9606)
		accession = uniprot["Accession"]
	files = [x for x in os.listdir(directory) if x.lower().startswith(accession.lower())]
	aaseq = compute_features.translate_seq(ntseq)
	imputed_seq = {}
	couplings = {}
	for f in files:
		df = pd.read_csv(os.path.join(directory, f), sep=",", index_col=None).T.to_dict()
		tmpcouplings = {}
		for k,v in df.items():
			imputed_seq[int(v["i"])] = v["A_i"]
			imputed_seq[int(v["j"])] = v["A_j"]
			if absval:
				tmpcouplings[(int(v["i"]), int(v["j"]))] = abs(float(v["cn"]))
			else:
				tmpcouplings[(int(v["i"]), int(v["j"]))] = float(v["cn"])
		if cdf:
			couplings.update(copy.deepcopy(stat_tools.cdf(tmpcouplings)))
		else:
			couplings.update(copy.deepcopy(tmpcouplings))
	seq = ""
	for i in range(len(aaseq)):
		if i+1 in imputed_seq.keys():
			seq += imputed_seq[i+1]
		else:
			seq += "X"
	aligned = tools.align_sequences({"seq":aaseq, "imputed":seq})
	new_couplings = {}
	for ((p1, p2),v) in couplings.items():
		q1, error1 = compute_features.convert_position(aligned["imputed"], aligned["seq"], p1)
		q2, error2 = compute_features.convert_position(aligned["imputed"], aligned["seq"], p2)
		if error1 == None and error2 == None:
			minq = min(q1, q2)
			maxq = max(q1, q2)
			new_couplings[(minq, maxq)] = copy.deepcopy(v)	
	return(new_couplings)

def pair_scores(variants, conservation, pairing, pair_only=False):
	if max(list(conservation)) <= 1.0 and min(list(conservation)) >= 0.0:
		conservation = {k:math.log(1-conservation[k]+0.001) for k in conservation.index}
	scores = {}
	for var, muts in variants.items():
		scores[var] = 0.0
		for pos1, nts1 in muts.items():
			if not pair_only:
				scores[var] += conservation[pos1]
			for pos2, nts2 in muts.items():
					pair = list((pairing[pairing[0] == pos1])[pairing[1] == pos2][2])
					if len(pair) >= 1:
						pair = pair[0]
						scores[var] += 0.5*pair*(conservation[pos1] + conservation[pos2])
	return(sorted(scores.items(), key=lambda kv: kv[0]))

#generates a clustering of sequences in an NT MSA, then finds examplar sequences from each cluster and uses examplars to compute conservation
def provean_nts(ntvariants, seqs, genename=None, outdir="./", space="-", frac=True, maxnum=None):
	muts = list(set([(pos, tuple(nts)) for muts in ntvariants.values() for (pos, nts) in muts.items()]))
	muts = {mut:[] for mut in muts}
	if not isinstance(seqs, (dict,)) and isinstance(seqs, (str,)):
		seqs = tools.read_fasta(seqs)
	if genename == None or len(genename) < 1:
		genename = list(seqs.keys())[-1]
	if maxnum != None and isinstance(maxnum, (int, float)):
		frac = float(maxnum)/len(seqs)
	write = True
	if outdir != None and len(outdir) > 0 and os.path.exists(os.path.join(outdir, genename + "_clust.pkl")):
		try:
			clusters = tools.pickleload(os.path.join(outdir, genename + "_clust.pkl"))
			write = False
		except Exception as e:
			print(e)
			write = True
	if write:
		newseqs, clusters = tools.cluster_seqs2(seqs, genename=genename, space=space, frac=frac)
		tools.pickledump(clusters, os.path.join(outdir, genename + "_clust.pkl"))
	clusters = {k:{k2:seqs[k2] for k2 in v} for k,v in clusters.items()}
	for k,v in clusters.items():
		clusters[k][genename] = seqs[genename]
	for k,clust in clusters.items():
		clust = tools.align_sequences({k:v.replace(space, "") for k,v in clust.items()})
		distribution = compute_distribution(clust, genename, space=space)
		for mut in muts:
			muts[mut].append(distribution[mut[0]][mut[1][1]])
	muts = {k:statistics.mean(v) for k,v in muts.items()}
	return(muts)

#queries VEP
def get_vep(ntvariants, genename, get_pairs=True, epsilon=0.001, outdir="./"):
	
	mutscores = {}
	pairings = {}
	muts = list(set([(pos, tuple(nts)) for muts in ntvariants.values() for (pos, nts) in muts.items()]))
	mutlist = [genename + ":c." + str(pos) + nts[0] + ">" + nts[1] for pos, nts in muts]

	if outdir != None and len(outdir) > 0 and os.path.exists(os.path.join(outdir, genename + "_VEP.pkl")) and os.path.exists(os.path.join(outdir, genename + "_LD.pkl")):
		transposed_mutscores = tools.pickleload(os.path.join(outdir, genename + "_VEP.pkl"))
		pairings = tools.pickleload(os.path.join(outdir, genename + "_LD.pkl"))
		stored_vars = transposed_mutscores["Polyphen"].keys()
		if set(muts).issubset(set(stored_vars)):
			return(transposed_mutscores, pairings)
		else:
			prev_transposed_mutscores = copy.deepcopy(transposed_mutscores)
			prev_pairings = copy.deepcopy(pairings)
			write = True
	else:
		write = True

	server = "https://rest.ensembl.org"
	ext = "/vep/human/hgvs"
	headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
	r = tools.retry_request(requests.post, positional_arguments=[server+ext], keyword_arguments={"headers":headers, "data":'{ "hgvs_notations" : ' + str(mutlist).replace("\'", "\"") + ' }'})


	if r == None or not r.ok:
		print(r.text)
		r.raise_for_status()
	else:
		decoded = r.json()
		init_time = time.time()
		for i in range(len(decoded)):
			mutname = decoded[i]["input"].split(":c.")[1]
			pos = int(max(re.findall("\d+", mutname), key=len))
			nts = re.findall("[A-Z]", mutname)
			mutname = (pos, tuple(nts[:2]))
			if mutname not in mutscores:
				mutscores[mutname] = {"Polyphen":[], "SIFT":[], "GnomAD":[], "dbSNP":None, "rsids":""}
				for j in range(len(decoded[i]["transcript_consequences"])):
					if "polyphen_score" in decoded[i]["transcript_consequences"][j]:
						mutscores[mutname]["Polyphen"].append(1-decoded[i]["transcript_consequences"][j]["polyphen_score"])
					if "sift_score" in decoded[i]["transcript_consequences"][j]:
						mutscores[mutname]["SIFT"].append(decoded[i]["transcript_consequences"][j]["sift_score"])
				if "colocated_variants" in decoded[i]:
					for j in range(len(decoded[i]["colocated_variants"])):
						if "frequencies" in decoded[i]["colocated_variants"][j]:
							try:
								mutscores[mutname]["GnomAD"].append(decoded[i]["colocated_variants"][j]["frequencies"][nts[1]]["gnomad"])
							except KeyError as e:
								pass
						try:
							if re.match("^rs\d+", decoded[i]["colocated_variants"][j]["id"]):
								mutscores[mutname]["rsids"] = decoded[i]["colocated_variants"][j]["id"]
						except:
							pass
			tools.update_time(i, len(decoded), init_time)
	if get_pairs:
		init_time = time.time()
		for i, mut1 in enumerate(mutscores):
			rsid1 = mutscores[mut1]["rsids"]
			for j in range(i):
				mut2 = list(mutscores.keys())[j]
				rsid2 = mutscores[mut2]["rsids"]
				
				togetherflag = False
				for muts in ntvariants.values():
					muts = [(k,tuple(v)) for k,v in muts.items()]
					if mut1 in muts and mut2 in muts:
						togetherflag = True
						break				
				if togetherflag:
					mutmin = min(mut1, mut2, key=lambda kv: kv[0])
					mutmax = max(mut1, mut2, key=lambda kv: kv[0])
					mutmin = mutmin[1][0] + str(mutmin[0]) + mutmin[1][1]
					mutmax = mutmax[1][0] + str(mutmax[0]) + mutmax[1][1]
					if (mutmin, mutmax) not in pairings and rsid1 != "" and rsid2 != "":
						ext = "/ld/human/pairwise/" + rsid1 + "/" + rsid2
						try:
							r2 = tools.retry_request(requests.get, positional_arguments=[server+ext], keyword_arguments={"headers":{ "Content-Type" : "application/json"}})
							time.sleep(1)
							if not r2.ok:
								raise Exception(r2.text)
						except Exception as e:
							print(e)
							continue
			
						decoded = r2.json()
						pairings[(mutmin, mutmax)] = 0.0
						for i in range(len(decoded)):
							try:
								pairings[(mutmin, mutmax)] = max(float(decoded[0]["r2"]), pairings[(mutmin, mutmax)])
							except KeyError as e:
								continue
						pairings[(mutmin, mutmax)] = math.sqrt(pairings[(mutmin, mutmax)])
					elif rsid1 != "" or rsid2 != "":
						pairings[(mutmin, mutmax)] = 0
			tools.update_time(i, len(mutscores), init_time)
	for mut in mutscores:
		if mutscores[mut]["dbSNP"] == None:
			try:
				rsid = max(re.findall("\d+", rsid))
				r = tools.retry_request(requests.get, ["https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/" + rsid])
				mutscores[mut]["dbSNP"] = parse_dict(r.json(), "clinical_significances")
				mutscores[mut]["dbSNP"] = statistics.mean([sigdict(val) for val in mutscores[mut]["dbSNP"]])
			except:
				mutscores[mut]["dbSNP"] = 0.5
		time.sleep(1)

	transposed_mutscores = {}
	for t in ["Polyphen","SIFT","GnomAD", "dbSNP"]:
		transposed_mutscores[t] = {}
		for mut in mutscores:
			transposed_mutscores[t][mut] = func_or_default(mutscores[mut][t], func=statistics.mean, default=0.5, verbose=False)
	if write:
		if 'prev_transposed_mutscores' in locals():
			for k in transposed_mutscores.keys():
				if k in prev_transposed_mutscores.keys():
					transposed_mutscores[k].update(copy.deepcopy(prev_transposed_mutscores[k]))
			del prev_transposed_mutscores
		if 'prev_pairings' in locals():
			for k in prev_pairings.keys():
				if k not in pairings:
					pairings.update(copy.deepcopy(prev_pairings[k]))
			del prev_pairings
		tools.pickledump(transposed_mutscores, os.path.join(outdir, genename + "_VEP.pkl"))
		tools.pickledump(pairings, os.path.join(outdir, genename + "_LD.pkl"))
	return(transposed_mutscores, pairings)	

#queries PROVEAN
def get_provean(aavariants, genename, aaseq, sleeptime=10, timeout=60, outdir="./"):
	aaseq = aaseq.replace("*", "").replace("-", "")
	muts = list(set([(pos, tuple(aas)) for muts in aavariants.values() for (pos, aas) in muts.items()]))
	mutlist = "\n".join([aas[0] + str(pos) + aas[1] for pos, aas in muts])
	
	if outdir != None and len(outdir) > 0 and os.path.exists(os.path.join(outdir, genename + "_PROVEAN.pkl")):
		mutscores = tools.pickleload(os.path.join(outdir, genename + "_PROVEAN.pkl"))
		if set(muts).issubset(set(mutscores.keys())):
			return(mutscores)
		else:
			prev_mutscores = copy.deepcopy(mutscores)
			write = True
	else:
		write = True
	r = tools.retry_request(requests.post, ["http://provean.jcvi.org/provean_seq_prg.php"], {'data':{"query":(">" + genename + "\n" + aaseq), "variant":mutlist, "database":"nr_sep_2012"}})
	try:	
		jobid = max(re.findall("job ID\: (\d+)", r.text), key=len)
	except Exception as e:
		with open(genename + "_test.html", "w") as outf:
			outf.write(r.text)
		print(aavariants)
	mutlist = mutlist.split("\n")
	success = False
	counter = 0
	mutscores = {}
	while counter * sleeptime < timeout*60 and not success:
		r = tools.retry_request(requests.get, ["http://provean.jcvi.org/provean_seq_report.php?jobid=" + jobid])
		with open("test.html", "w") as outf:
			outf.write(r.text)
		if '<div id="middle_content_template">' in r.text and 'PROVEAN Result' in r.text:
			table = r.text.split('<div id="middle_content_template">')[1]
			table = table.split("</div>")[0]
			rows = tools.get_element(r.text, k="tr")
			tds = [tools.get_element(row, k="td") for row in rows]
			tds = [td for td in tds if any([mutlist[i] in td[j] for i in range(len(mutlist)) for j in range(len(td))])]
			tds = [[tools.get_element(x, "p")[0].replace("<p>", "").replace("</p>", "") for x in td] for td in tds]
			for td in tds:
				if len(td) >= 2:
					pos, nts = tools.get_mutation_data(td[0])
					mutscores[(pos, tuple(nts))] = float(td[1])
			success = True
		else:
			counter = counter + 1
			time.sleep(sleeptime)
	if not success:
		return(None)
	if write:
		if 'prev_mutscores' in locals():
			mutscores.update(prev_mutscores)
		tools.pickledump(mutscores, os.path.join(outdir, genename + "_PROVEAN.pkl"))
	return(mutscores)

#queries dbSNP
def get_dbSNP(geneset, accids={}):
	if accids == None or len(accids) < 1:
		accids = {k:bio_tools.get_accids(k)["mRNA"][0][0] for k,pholder,v in geneset}
	data = {}
	init_time = time.time()
	for i, (gene, directory, f) in enumerate(geneset):
		df = pd.read_csv(f, sep="\t", index_col=0)
		variants = {}
		for var in [x for x in df.columns if str(x).startswith("c.")]:
			try:
				ids = tools.get_mutalyzer(accids[gene], str(var).split("c.")[-1], "c.")
				variants[str(var)] = bio_tools.query_dbSNP([x for x in ids if x.startswith("NC")][0])
			except Exception as e:
				variants[str(var)] = [[],[],[],[]]
				print(e)
		print(variants)
		df = df.T.to_dict()
		for k1,v1 in df.items():
			try:
				dataset = [vi for vi in v1 if vi.startswith("c.") and v1[vi] > 0]
				dataset = {vi:variants[vi] for vi in dataset}
				newdataset = []
				for ki,vi in dataset.items():
					if len(vi) > 1 and "pathogenic" in [str(x).lower() for x in vi[1]]:
						newdataset.append(ki)
				data[(f, k1)] = newdataset
			except Exception as e:
				data[(f, k1)] = []
				print(str(e))
				print(dataset)
		tools.update_time(i, len(geneset), init_time)
	return(data)

#finds best threshold (by BAS) given x scores and y categorization
def find_opt_thresh(x, y, invert=False):
	tmpx = list(x)
	tmpy = list(y)
	#tmpx = [x[i] for i in range(len(tmpx)) if tools.is_numeric(tmpx[i]) and not np.isnan(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpy[i])]
	#tmpy = [y[i] for i in range(len(tmpy)) if tools.is_numeric(tmpx[i]) and not np.isnan(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpy[i])]
	allvalues = []
	for xi in sorted(list(set(tmpx))):
		bas = balanced_accuracy_score(tmpx, tmpy, thresh=xi, invert=invert)
		allvalues.append((xi, bas))
	return(allvalues)

#converts scores x to categorization and computes BAS
def balanced_accuracy_score(x,y,thresh, invert=False):
	tmpx = list(x)
	tmpy = list(y)
	#tmpx = [x[i] for i in range(len(tmpx)) if tools.is_numeric(tmpx[i]) and not np.isnan(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpy[i])]
	#tmpy = [y[i] for i in range(len(tmpy)) if tools.is_numeric(tmpx[i]) and not np.isnan(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpy[i])]
	if not invert:
		pred = [1 if tmpx[i] >= thresh else 0 for i in range(len(tmpx))]
	else:
		pred = [1 if tmpx[i] <= thresh else 0 for i in range(len(tmpx))]

	return(sklearn.metrics.balanced_accuracy_score(tmpy, pred))

#randomly samples dataset num times and computes AUC and BAS for each subsample
def compute_subsets(f="merged_outputs.xlsx", response="Disease", num=1000, multonly=True, cols = {"Polyphen": True, "SIFT": True, "PROVEAN":True, "GnomAD": True, "Conservation (with PDB)": True, "AA conservation": True, "AA entropy": False, "AA variance": False, "BLOSUM conservation": False, "Rare codon enrichment": True, "Consurf (wo 3D info)": True, "Accessible surface area (PDBPisa)": False, "Solvation energy (PDBPisa)": False, "NT conservation": True, "NT entropy": False, "NT variance": False, "Relative surface area (NetsurfP2)": False, "MSA loglikelihood (AA)": True, "MSA loglikelihood (NT)": True, "EVmutation (AA)": True, "EVmutation (NT)": True, "Number of mutants": False}):
	wb = openpyxl.load_workbook(f, read_only=True, keep_links=False)
	sheetnames = list(wb.sheetnames)
	auc = {}
	bas = {}
	pre_rec = {}
	summary = {}
	init_time = time.time()
	sheets = {}
	for i, sheetname in enumerate(sheetnames):
		sheets[sheetname] = pd.read_excel(f, sheet_name=sheetname, index_col=0)
	for i, sheetname in enumerate(sheets.keys()):
		df = sheets[sheetname]
		if multonly:
			df = df.loc[df["Number of mutants"] >= 2]
		auc[sheetname] = {k:[] for k in cols}
		bas[sheetname] = {k:[] for k in cols}
		pre_rec[sheetname] = {k:[] for k in cols}
		summary[sheetname] = {k:{} for k in cols}
		for j in range(num):
			frac = np.random.uniform(low=0.3, high=0.7)
			rand = df.sample(frac = frac, replace=False)
			complement = df[~df.index.isin(rand.index)]
			tmpy = [float(x) for x in list(rand[response])]
			complementy = [float(x) for x in list(complement[response])]
			for col in cols:
				try:
					tmpx = [float(x) for x in list(rand[col])]
					auc[sheetname][col].append(combine_documents.plot_roc("pholder", tmpx, tmpy, invert=cols[col], plot=False))
					thresh = max(find_opt_thresh(tmpx, tmpy, invert=cols[col]), key=lambda kv: kv[1])[0]
					complementx = [float(x) for x in list(complement[col])]
					bas[sheetname][col].append(balanced_accuracy_score(complementx, complementy, thresh, invert=cols[col]))
				except ValueError as e:
					pass
				except Exception as e:
					traceback.print_exc()
					print(e)
			tools.update_time(i*num + j, len(sheetnames) * num, init_time)
	summary = {"AUC":{k:summarize_rand_data(auc[k]) for k in auc}, "BAS":{k:summarize_rand_data(bas[k]) for k in bas}}
	writer = pd.ExcelWriter("./randomized_analysis.xlsx", engine="openpyxl")
	for score, v in summary.items():
		for sheet, v2 in v.items():
			pd.DataFrame(v2).T.to_excel(writer, score + "_" + sheet)
	writer.save()
	writer.close()
	return(auc, bas, summary)

#computes summary statistics for data
def summarize_rand_data(data):
	summary = {k:{} for k in data.keys()}
	for k,v in data.items():
		try:
			summary[k] ={"Min":min(v), "5th %ile":np.percentile(v, 5), "Mean":statistics.mean(v), "Median":statistics.median(v), "95th %ile":np.percentile(v, 95),"Max":max(v), "StDev":statistics.stdev(v)}
		except ValueError as e:
			summary[k] = {"Min":"", "5th %ile":"", "Mean":"", "Median":"", "95th %ile":"", "Max":"", "StDev":""}
	return(summary)

def split_long(s, lim=10):
	news = ""
	for si in s.split():
		if len(news) < 1 or len(news.split("\n")[-1]) + 1 + len(si) < lim:
			news += " " + si
		else:
			news += "\n" + si
	return(news)

#generates box-and-whiskers plot for randomized trials
def box_whiskers_subsets(data, lim=15, x=None, y=20, bottom=0.35):
	plt.cla()
	plt.clf()
	if x == None:
		x = len(data)*1.25
	print(str(x) + "\t" + str(y))
	data = {k:v for k,v in data.items() if len(v) > 0}
	df = pd.DataFrame({split_long(k, lim=lim):{i:v[i] for i in range(len(v))} for k,v in data.items()})
	cols = sorted([x for x in df.columns], key = lambda kv: statistics.median(df[kv]), reverse=True)
	df = df[cols]
	df.boxplot(figsize=(x, y), rot=90)
	plt.gcf().subplots_adjust(bottom=bottom)
	plt.tick_params(labelsize=8, pad=6);
	plt.show()

#combines seperate input files (which contains genetic variants) from different genes/diseases
def merge_inputs(genesets, outfile="./merged_inputs.tsv"):
	data = {}
	for dataset in genesets:
		df = pd.read_csv(dataset[2], sep="\t", index_col=0).T.to_dict()
		if "/" in dataset[1]:
			dataset = (dataset[0], dataset[1].split("/")[-1], dataset[2])
		for k,v in df.items():
			if "Neutral" in k:
				k = dataset[1] + "_" + str(k)
			if k in data:
				if len(k.split("_")) < 2:
					k = k + "_2"
				else:
					k = "_".join(k.split("_")[:-1]) + "_" + str(int(k.split("_")[-1])+1)
			data[k] = copy.deepcopy({k2:v2 for k2, v2 in v.items()})
			data[k]["Gene"] = dataset[0]
			data[k]["Source"] = dataset[2]
			variants = []
			for k2, v2 in list(data[k].items()):
				if k2.startswith("c.") and tools.is_numeric(v2) and float(v2) > 0:
					variants.append(dataset[0] + ":" + k2)
			data[k]["Variants"] = frozenset(variants)
			if len(variants) < 2:
				del data[k]

	pd.DataFrame(data).T.to_csv(outfile, sep="\t")
	return(data)

#removes rows from output (score) files if they have the same input (combination of variants)
def remove_dups(outf, inf, outfile="./merged_output_nodups"):
	if outf.endswith(".xlsx") or outf.endswith(".xls"):
		wb = openpyxl.load_workbook(outf, read_only=True, keep_links=False)
		sheetnames = list(wb.sheetnames)
		outdata = pd.read_excel(outf, sheet_name=sheetnames[0], index_col=0)
	else:
		outdata = pd.read_csv(outf, sep="\t", index_col=0)
	if inf.endswith(".xlsx") or outf.endswith(".xls"):
		wb = openpyxl.load_workbook(inf, read_only=True, keep_links=False)
		sheetnames = list(wb.sheetnames)
		indata = pd.read_excel(inf, sheet_name=sheetnames[0], index_col=0)
	else:
		indata = pd.read_csv(inf, sep="\t", index_col=0)
	indata = indata.T.to_dict()
	dups = {varset:[k for k in indata if indata[k]["Variants"] == varset] for varset in [indata[k]["Variants"] for k in indata]}
	toremove = {varset:dups[varset][1:] for varset in dups if len(dups[varset]) > 1}
	toremove = [x for v in toremove.values() for x in v]
	toremove = [x for x in toremove if x in outdata.index]
	innotout = [x for x in indata.keys() if x not in outdata.index]
	outnotin = [x for x in outdata.index if x not in indata.keys()]
	
	if outf.endswith(".xlsx") or outf.endswith(".xls"):
		writer = pd.ExcelWriter(outfile + ".xlsx", engine="openpyxl")
		for sheet in sheetnames:
			outdata = pd.read_excel(outf, sheet_name=sheet, index_col=0)
			outdata = outdata.drop(toremove, axis=0)
			outdata.to_excel(writer, sheet)
		writer.save()
		writer.close()
	else:
		outdata = outdata.drop(toremove, axis=0)
		outdata = outdata.to_csv(outfile + ".tsv", sep="\t")
	return(toremove, outnotin, innotout)

#combines outputs (scores) across datasets
def merge_computed(genesets, infile="", outfile="./merged_outputs.xlsx"):
	data = {}
	for dataset in genesets:
		#print(dataset)
		if infile != None and len(infile) > 0:
			if isinstance(dataset, (list, tuple,)) and len(dataset) > 1:
				dataset = os.path.join(os.path.dirname(dataset[2]), dataset[1], infile)
		name = os.path.basename(os.path.dirname(dataset))
		wb = openpyxl.load_workbook(dataset, read_only=True, keep_links=False)
		sheetnames = list(wb.sheetnames)
		for sheet in sheetnames:
			if sheet not in data:
				data[sheet] = {}
			df = pd.read_excel(dataset, sheet_name=sheet, index_col=0)
			df = df.T.to_dict()
			for k,v in df.items():
				if "Neutral" in k:
					k = name + "_" + str(k)
				if k in data[sheet]:
					if len(k.split("_")) < 2:
						k = k + "_2"
					else:
						k = "_".join(k.split("_")[:-1]) + "_" + str(int(k.split("_")[-1])+1)
				data[sheet][k] = {k2:v2 for k2, v2 in v.items() if not re.match("[A-Z]\d+[A-Z]", k2)}
				data[sheet][k]["Name"] = name
				data[sheet][k]["Source"] = dataset

	writer = pd.ExcelWriter(outfile, engine="openpyxl")
	for sheet, v in data.items():
		pd.DataFrame(v).T.to_excel(writer, sheet)
	writer.save()
	writer.close()
	return(data)

def write_pairs(pairfile):
	with open(pairfile, "rb") as pklfile:
		pairs = pickle.load(pklfile)
	muts = list(set([k[0] for k in pairs] + [k[1] for k in pairs]))
	muts = sorted(muts, key=lambda k: int(max(re.findall("\d+", k))))
	mutnames = [mut[1:-1]+mut[0] + ">" + mut[-1] for mut in muts]

	arr = [[""] + mutnames]
	for i in range(len(muts)):
		row = [mutnames[i]]
		for j in range(len(muts)):
			if muts[i] == muts[j]:
				row += [1]
			else:
				try:
					row += [pairs[(muts[i], muts[j])]]
				except:
					row += [""]
		arr.append(row)
	with open(os.path.splitext(pairfile)[0] + ".tsv", "w") as outf:
		for i in range(len(arr)):
			for j in range(len(arr[i])):
				outf.write(str(arr[i][j]) + "\t")
			outf.write("\n")
'''
def get_imutant(aaseq, aavarset, path=ram_disk):
	mutscores = {}
	with open(path + "tmp.seq", "w") as outf:
		outf.write(aaseq)
	for (pos, aas) in aavarset:
		try:
			output = subprocess.getoutput("python2 -O " + i_mutant_cmd + "  -seqv " + path + "tmp.seq " + str(pos) + " " + aas[1])
			value = max(re.findall("Position   WT  NEW     DDG    pH    T\s+" + str(pos) + "\s+" + aas[0] + "\s+" + aas[1] + "\s+" + "(\S+)", output), key=len)
			mutscores[(pos, aas)] = float(value)
		except Exception as e:
			mutscores[(pos, aas)] = np.nan
	return(mutscores)

def get_maestro(pdbfile, chain, aavarset):
	mutscores = {"MAESTRO score":{}, "MAESTRO dG":{}}
	for (pos, aas) in aavarset:
		mutscores["MAESTRO score"][(pos, tuple(aas))] = {}
		mutscores["MAESTRO dG"][(pos, tuple(aas))] = {}
		output = subprocess.getoutput("maestro /home/dholcomb/software/MAESTRO_linux_x64/config.xml " + pdbfile + " --evalmut=\'" + aas[0] + str(pos) + "." + str(chain) + "{" + aas[1] + "}\' --bu")
		try:
			mutdata = output.split("\n")[2].split("\t")
			mutscores["MAESTRO score"][(pos, tuple(aas))] = float(mutdata[4])
			mutscores["MAESTRO dG"][(pos, tuple(aas))] = float(mutdata[5])
		except Exception as e:
			mutscores["MAESTRO score"][(pos, tuple(aas))] = np.nan
			mutscores["MAESTRO dG"][(pos, tuple(aas))] = np.nan
	return(mutscores)

def get_foldx(pdbfile, chain, aavariants, path=ram_disk):
	pdbname = os.path.splitext(os.path.basename(pdbfile))[0]
	excluded_muts = []
	aavarset = sorted(set([(pos, (aas[0], aas[1])) for v in aavariants.values() for pos, aas in v.items()]))
	if len(aavarset) > 1:
		for (pos, aas) in aavarset:
			score = get_foldx(pdbfile, chain, {"test":{pos:[aas[0], aas[1]]}}, path=path)
			if not isinstance(score, (dict)) or "test" not in score:
				excluded_muts.append((pos, tuple(aas)))
	try:
		mutscores = {}
		with open(path + "individual_list.txt", "w") as outf:
			for i, (var, muts) in enumerate(aavariants.items()):
				mutset = set([(pos, tuple(aas)) for pos, aas in muts.items()]) - set(excluded_muts)
				for j, (pos, aas) in enumerate(mutset):
					outf.write(aas[0] + str(chain) + str(pos) + aas[1])
					if not j == len(mutset) - 1:
						outf.write(",")
				if len(mutset) > 0:
					outf.write(";\n")
		output = subprocess.getoutput("foldx --command=BuildModel --pdb=" + pdbfile + " --mutant-file=" + path + "individual_list.txt")
		if "Specified residue not found." in output:
			raise ValueError("Specified residue not found.")
		total_energy = []
		with open("Dif_"+pdbname+".fxout", "r") as inf:
			for i, line in enumerate(inf.readlines()):
				if i > 8 and len(line.split("\t")) >= 23:
					total_energy.append(float(line.split("\t")[1]))
		rowind = 0
		for i, (var, muts) in enumerate(aavariants.items()):
			mutset = set([(pos, tuple(aas)) for pos, aas in muts.items()]) - set(excluded_muts)
			if len(muts) < 1:
				mutscores[var] = 0
			elif len(mutset) < 1:
				mutscores[var] = np.nan
			else:
				mutscores[var] = total_energy[rowind]
				rowind += 1
	except Exception as e:
		mutscores = {}
	finally:
		for f in os.listdir("./"):
			if pdbname in f and os.path.basename(pdbfile) != os.path.basename(f):
				os.remove(f)
		return(mutscores)
'''
#returns either a function on a value, or some default value if error
def func_or_default(l, func=statistics.mean, default=0.5, verbose=False):
	try:
		return(func(l))
	except Exception as e:
		if verbose:
			print(e)
		return(default)

#parses the clinical_significance from dbSNP
def parse_dict(d, label="clinical_significances", path=[]):
	if isinstance(d, (dict,)):
		l = list(d.keys())
	elif isinstance(d, (list, tuple,)):
		l = range(len(d))
	else:
		return([])
	returnval = []
	for k in l:
		v = d[k]
		try:
			if isinstance(v, dict) and label in v:
				if isinstance(v[label], (list,)):
					returnval += v[label]
				elif isinstance(v[label], (tuple,)):
					returnval += list(v[label])
				else:
					returnval += [v[label]]
			if isinstance(v, (dict, list, tuple,)):
				val = parse_dict(v, label=label, path = path + [k])
				if val != None:
					returnval += [val]
		except Exception as e:
			print(e)
			continue
	returnval = recur_merge(returnval)
	return(returnval)

def recur_merge(l):
	return_val = []
	if isinstance(l, (list, tuple,)):
		for i in range(len(l)):
			return_val += recur_merge(l[i])
	else:
		return([l])
	return(return_val)

#general function to score each individual variant
#mutscores represents the variable containing scores
def score_variant(muts, mutscores, pairing, mode="add", invert=False, pctflag=False, hetero=False, epsilon=0.001, genename="", parameter="", write=True):
	scores = {}
	for pos, aas in muts.items():
		pos = int(pos)
		try:
			scores[(pos, tuple(aas))] = mutscores[(pos, tuple(aas))] #try (position, (wt NT, mut NT)) as input
		except (KeyError,IndexError,TypeError) as e:
			try:
				scores[(pos, tuple(aas))] = mutscores["EVmutation(NT) " + aas[1]][pos] #try EVmutation (with variant NT)
			except (KeyError,IndexError,TypeError) as e:
				try:
					scores[(pos, tuple(aas))] = mutscores["EVmutation " + aas[1]][pos] #try EVmutation (with variant AA)
				except (KeyError,IndexError,TypeError) as e:
					try:
						scores[(pos, tuple(aas))] = mutscores[pos][aas[1]] #try [position][mut NT] as input
					except (KeyError,IndexError,TypeError) as e:
						try:
							scores[(pos, tuple(aas))] = mutscores[pos] #try only position
						except (KeyError,IndexError,TypeError) as e:
							try:
								scores[(pos, tuple(aas))] = mutscores[(aas[0], aas[1])] #try (wt NT, mut NT)
							except (KeyError,IndexError,TypeError) as e: #other error
								pass
	scores = {mut:float(score) for mut,score in scores.items() if tools.is_numeric(score) and not np.isnan(float(score))}
	if write:
		with open("./variant_scores.tsv", "a") as outf:
			for k,v in scores.items():
				outf.write(str(parameter) + "\t" + str(genename) + ":" + str(k[0]) + str(k[1][0]) + ">" + str(k[1][1]) + "\t" + str(v) + "\n")
	if scores == None or len(scores) < 1:
		return("")
	if pctflag: #score represents a likelihood, convert to loglikelihood
		if invert:
			scores = {mut:math.log(1-score + epsilon) for mut, score in scores.items()}
		else:
			scores = {mut:math.log(score + epsilon) for mut, score in scores.items()}
	scores = {mut:score * (1-pairing[mut]/2) for mut, score in scores.items()} #account for position pairing (pairing[mut]=0 is independent model)
	if mode.lower() in ["additive", "add"]: #additive model
		return(sum(scores.values()))
	elif mode.lower() in ["max", "maximum"]: #maximum model
		if invert or pctflag:
			return(min(scores.values()))
		else:
			return(max(scores.values()))
	elif mode.lower() in ["maxminmedian", "minmaxmedian", "mmm"]: #other experimental model
		return(max(scores.values()) + min(scores.values()) - statistics.mean(scores.values()))
	
#analyze outputs (scores) for AUC, Mann-Whitney or Chi-square
def analyze_data(datafile, response={"Activity":100, "Expression":100}, pairmodel="ind", combine_model="add", now=datetime.datetime.now(), multonly=True, sep="\t"):
	invert = {"Polyphen": True, "SIFT": True, "PROVEAN": True, "PROVEAN (NT)":True, "GnomAD": True, "Conservation (with PDB)": True, "AA percent identity": True, "AA entropy": False, "AA variance": False, "BLOSUM conservation": False, "Rare codon enrichment": True, "Consurf (wo 3D info)": True, "Accessible surface area (PDBPisa)": False, "Solvation energy (PDBPisa)": False, "NT percent identity": True, "NT entropy": False, "NT variance": False, "Relative surface area (NetsurfP2)": False, "MSA loglikelihood (AA)": True, "MSA loglikelihood (NT)": True, "EVmutation (AA)": True, "EVmutation (NT)": True, "Number of mutants": False}
	analysis = {}
	if datafile.endswith(".xls") or datafile.endswith(".xlsx"):
		try:
			scores = pd.read_excel(datafile, sheet_name=combine_model + pairmodel, index_col=0)
		except Exception as e:
			print(e)
	else:
		try:
			scores = pd.read_csv(os.path.splitext(datafile)[0] + "_" + combine_model + pairmodel + "_" + str(now.year) + str(now.month).rjust(2,'0') + str(now.day).rjust(2,'0') + ".tsv", sep=sep, index_col=0)
		except Exception as e:
			scores = pd.read_csv(datafile, sep=sep, index_col=0)
	if multonly:
		scores = scores.loc[scores["Number of mutants"] >= 2]
	scores = scores.to_dict()
	for col in list(scores):
		numeric = [k for k in scores[col] if tools.is_numeric(scores[col][k]) and not np.isnan(float(scores[col][k]))]
		if len(numeric) == 0:
			valset = list(set([scores[col][k] for k in scores[col] if isinstance(scores[col][k], str)]))
		else:
			valset = []
		for val in valset:
			scores[val] = {}
			for k in scores[col]:
				if scores[col][k] == val:
					scores[val][k] = 1
				else:
					scores[val][k] = 0
	for col in scores:
		if col not in response.keys() and not re.match("[A-Za-z]\d+[A-Za-z]", str(col)):
			numeric = [k for k in scores[col] if tools.is_numeric(scores[col][k]) and not np.isnan(scores[col][k])]
			analysis[col] = {"Number distinct":len(set([scores[col][k] for k in numeric])), "Number missing":len(scores[col].keys()) - len(numeric)}
			try:
				with warnings.catch_warnings():
					warnings.simplefilter("ignore")
					analysis[col]["Shapiro-Wilk p"] = scipy.stats.shapiro([scores[col][k] for k in numeric])[1]
			except:
				analysis[col]["Shapiro-Wilk p"] = ""
			for respcol in response.keys():
				rows = [k for k in numeric if k in scores[respcol] and tools.is_numeric(scores[respcol][k]) and not np.isnan(scores[respcol][k])]
				x = [float(scores[col][k]) for k in rows]
				y = [float(scores[respcol][k]) for k in rows]
				if response[respcol] != None and tools.is_numeric(response[respcol]):
					y = [y[k] - response[respcol] for k in range(len(y))]
				if len(set(y)) == 2:
					if len(set(x)) == 2:
						f = []
						for xi in [min(x), max(x)]:
							f.append([])
							for yi in [min(y), max(y)]:
								f[-1].append(len([k for k in scores[col] if scores[col][k] == xi and scores[respcol][k] == yi]))
						try:
							chisq = scipy.stats.fisher_exact(f)
							analysis[col][respcol + " Chisquare"] = chisq[0]
							analysis[col][respcol + " Chisquare p"] = chisq[1]
							analysis[col][respcol + " AUC"] = combine_documents.plot_roc(col + "_" + combine_model + pairmodel, x, y, outdir=os.path.dirname(datafile), plot=False, invert=(col in invert and invert[col]))
						except Exception as e:
							print(e)
							traceback.print_exc()
							analysis[col][respcol + " Chisquare"] = ""
							analysis[col][respcol + " Chisquare p"] = ""
							analysis[col][respcol + " AUC"] = ""
					else:
						x1 = [x[i] for i in range(len(x)) if y[i] == min(set(y))]
						x2 = [x[i] for i in range(len(x)) if y[i] == max(set(y))]
						try:
							mw = scipy.stats.mannwhitneyu(x1, x2)
							analysis[col][respcol + " Mann-Whitney"] = mw.statistic
							analysis[col][respcol + " Mann-Whitney p"] = mw.pvalue
							analysis[col][respcol + " stdev diff"] = (statistics.median(x2) - statistics.median(x1))/((statistics.stdev(x1) + statistics.stdev(x2))/2)
							analysis[col][respcol + " AUC"] = combine_documents.plot_roc(col + "_" + combine_model + pairmodel, x, y, outdir=os.path.dirname(datafile), plot=False, invert=(col in invert and invert[col]))
						except ValueError:
							analysis[col][respcol + " Mann-Whitney"] = ""
							analysis[col][respcol + " Mann-Whitney p"] = ""
							analysis[col][respcol + " stdev diff"] = ""
							analysis[col][respcol + " AUC"] = ""
						except Exception as e:
							print(e)
							traceback.print_exc()
							analysis[col][respcol + " Mann-Whitney"] = ""
							analysis[col][respcol + " Mann-Whitney p"] = ""
							analysis[col][respcol + " stdev diff"] = ""
							analysis[col][respcol + " AUC"] = ""
						try:
							rows = [k for k in scores[respcol] if tools.is_numeric(scores[respcol][k]) and not np.isnan(scores[respcol][k])]
							notnarows = [k for k in rows if tools.is_numeric(scores[col][k]) and not np.isnan(scores[col][k])]
							narows = [k for k in rows if not k in notnarows]
							mat = [[len([k for k in notnarows if scores[respcol][k] == min(set(y))]), len([k for k in notnarows if scores[respcol][k] == max(set(y))])], [len([k for k in narows if scores[respcol][k] == min(set(y))]), len([k for k in narows if scores[respcol][k] == max(set(y))])]]
							mw = scipy.stats.chisquare(mat, axis=None)					
							analysis[col][respcol + " NA vs not chi-square"] = mw.statistic
							analysis[col][respcol + " NA vs not chi-square p-value"] = mw.pvalue
						except Exception as e:
							analysis[col][respcol + " NA vs not chi-square"] = ""
							analysis[col][respcol + " NA vs not chi-square p-value"] = ""
				else:
					if len(set(x)) == 2:
						y1 = [y[i] for i in range(len(y)) if x[i] == min(set(x))]
						y2 = [y[i] for i in range(len(y)) if x[i] == max(set(x))]
					try:
						sm = scipy.stats.spearmanr(x, y)						
						analysis[col][respcol + " Spearman"] = sm.correlation
						analysis[col][respcol + " Spearman p"] = sm.pvalue
					except Exception as e:
						print(e)
						traceback.print_exc()
						analysis[col][respcol + " Spearman"] = ""
						analysis[col][respcol + " Spearman p"] = ""

	for respcol in response:
		if len(set(scores[respcol].values())) == 2:
			cols = [col for col in analysis if respcol + " AUC" in analysis[col] and tools.is_numeric(analysis[col][respcol + " AUC"])]
			subset = [analysis[col][respcol + " AUC"] for col in cols]
		else:
			cols = [col for col in analysis if respcol + " Spearman" in analysis[col] and tools.is_numeric(analysis[col][respcol + " Spearman"])]
			subset = [analysis[col][respcol + " Spearman"] for col in cols]
		try:
			mad = float(pd.DataFrame(subset).mad())
			maxabs = max(subset, key=abs)
			maxcol = cols[subset.index(maxabs)]
			#print(respcol + "\t" + str(statistics.median(subset)) + "\t" + str(mad) + "\t" + str(maxabs) + "\t" + maxcol)
		except:
			pass
	pd.DataFrame(analysis).T.to_csv(os.path.splitext(datafile)[0] + "_analyzed_" + combine_model + pairmodel + "_" + str(now.year) + str(now.month).rjust(2,'0') + str(now.day).rjust(2,'0') + ".tsv", sep="\t")

def regress_variants(mutdf, response, normal=None, genename="", write=False):
	if genename == None or len(genename) < 1:
		genename = os.path.splitext(os.path.basename(mutdf))[0]
	df = pd.read_csv(mutdf, sep="\t", index_col=0)
	x = df[[str(col) for col in list(df.columns) if str(col).startswith("c.")]]
	y = df[response]
	if normal != None and tools.is_numeric(normal):
		y = y - normal
	responseset = set(list(y))
	colnames = list(x.columns)
	imp = SimpleImputer(strategy="median")
	x = np.array(imp.fit_transform(x))
	y = np.array(imp.fit_transform([[k] for k in list(y)]))

	if len(responseset) > 2:
		model = sm.OLS(y,x).fit()
	elif len(responseset) == 2:
		model = sm.Logit(y,x).fit()
	
	results_as_html = model.summary(xname=colnames).tables[1].as_html()
	results = pd.read_html(results_as_html, header=0, index_col=0)[0]
	
	if genename != None and len(genename) > 0:
		muts = [str(col) for col in list(df.columns) if str(col).startswith("c.")]
		ntvars = {"All":{int(max(re.findall("\d+", mut), key=len)):re.findall("[ACGT]", mut) for mut in muts}}
		transposed_mutscores, pairings = get_vep(ntvars, genename, get_pairs=False, outdir=os.path.join(os.path.dirname(mutdf), genename))

		results = results.to_dict()
		transposed_mutscores = {p:{"c."+str(pos)+nts[0]+">"+nts[1]: val for ((pos, nts), val) in d.items()} for p,d in transposed_mutscores.items()}
		results.update(transposed_mutscores)
		results = pd.DataFrame(results)
		correlations = {}
		for k in transposed_mutscores:
			x = [list(results[k])[i] for i in range(len(list(results[k]))) if not tools.is_nan(list(results["coef"])[i]) and not tools.is_nan(list(results[k])[i])]
			y = [list(results["coef"])[i] for i in range(len(list(results["coef"]))) if not tools.is_nan(list(results["coef"])[i]) and not tools.is_nan(list(results[k])[i])]
			smr = scipy.stats.spearmanr(x,y)
			correlations[k] = {"r":smr.correlation, "p":smr.pvalue}
	
	newexcel = os.path.join(os.path.dirname(mutdf), genename, genename + "_" + response + "_regression.xlsx")
	writer = pd.ExcelWriter(newexcel, engine="openpyxl")
	results.to_excel(writer, "regression")
	pd.DataFrame(correlations).T.to_excel(writer, "correlations")	
	writer.save()
	writer.close()
	
def regress_scores(mutdf, response, normal=None, genename="", write=False):
	columns = ["Polyphen","SIFT","PROVEAN","PROVEAN (NT)", "GnomAD","Conservation (with PDB)","AA percent identity","AA entropy", "AA variance", "BLOSUM conservation", "Rare codon enrichment","Consurf (wo 3D info)","Accessible surface area (PDBPisa)", "Solvation energy (PDBPisa)", "NT percent identity","NT entropy","NT variance","Relative surface area (NetsurfP2)", "MSA loglikelihood (AA)","MSA loglikelihood (NT)","EVmutation (AA)","EVmutation (NT)","Number of mutants"]
	if genename == None or len(genename) < 1:
		genename = os.path.splitext(os.path.basename(mutdf))[0]
	wb = openpyxl.load_workbook(mutdf, read_only=True, keep_links=False)
	sheetnames = list(wb.sheetnames)
	sheets = {}
	for count in range(len(sheetnames)):
		try:
			sheets[count] = pd.read_excel(mutdf, sheet_name=count, header=0, index_col=0)
		except Exception as e:
			print(e)
	newexcel = os.path.join(os.path.dirname(mutdf), genename + "_" + response + "_regression.xlsx")
	writer = pd.ExcelWriter(newexcel, engine="openpyxl")
	for k,sheet in sheets.items():
		try:
			sheetname = sheetnames[k]
			x = sheet[[str(col) for col in list(sheet.columns) if col in columns]]
			y = sheet[response]
			if normal != None and tools.is_numeric(normal):
				y = y - normal
			responseset = set(list(y))
			colnames = [i for i in list(x.columns) if any([True for j in list(sheet[i]) if tools.is_numeric(j) and not np.isnan(j)])]
			imp = SimpleImputer(strategy="median")
			x = np.array(imp.fit_transform(x))
			y = np.array(imp.fit_transform([[k] for k in list(y)]))

			if len(responseset) > 2:
				model = sm.OLS(y,x).fit()
			elif len(responseset) == 2:
				model = sm.Logit(y,x).fit()
			results_as_html = model.summary(xname=colnames).tables[1].as_html()
			results = pd.read_html(results_as_html, header=0, index_col=0)[0]
			results.to_excel(writer, sheetname)
		except Exception as e:
			print(e)
	writer.save()
	writer.close()

def variant_scores(f="./variant_scores.tsv", seqdir="./gene_data/"):
	data = {}
	with open(f) as inf:
		lines = inf.readlines()
	for l in lines:
		pieces = l.strip().split("\t")
		if len(pieces) >= 3:
			if pieces[1] not in data:
				data[pieces[1]] = {pieces[0]:pieces[2]}
			else:
				data[pieces[1]][pieces[0]] = pieces[2]
	for k,v in list(data.items()):
		gene = k.split(":")[0]
		aavar = k.split(":")[-1]
		if "GnomAD" not in v or not tools.is_numeric(v["GnomAD"]) or np.isnan(float(v["GnomAD"])):
			seq = tools.read_fasta(os.path.join(seqdir, gene + ".fasta"))["ORF"].replace("-", "")
			for k2 in data.keys():
				ntvar = k2.split(":")[-1]
				if k2.split(":")[0] == gene and re.match("\d+[ACGT]\>[ACGT]", ntvar):
					ntvar = tools.get_mutation_data(ntvar)
					impact = tools.get_mutant_aa(ntvar, seq, aaseq=compute_features.translate_seq(seq), index=1)
					aamut = tools.get_mutation_data(aavar)
					if aamut[0] == impact[0] and aamut[1][0] == impact[1][0] and aamut[1][1] == impact[1][1]:
						data[k2].update({k3:v3 for k3, v3 in data[k].items()})
						del data[k]
						break
			
	pd.DataFrame(data).T.to_csv("./variant_scores_(processed).tsv", sep="\t")
	return(data)
	
#scores all patients within a dataset
def score_variants(mutdf, datadf, aafasta, ntfasta, pdb="6qig.pdb", chain="A", pairmodel="ind", combine_model="add", genename="", aapairing="", ntpairing="", output="./ADAMTS13_mut_scores.tsv", hetero=False, multonly=True, dups=False, respcols={"Activity":100, "Expression":100}):

	transition = {("C","T"):1, ("A","G"):1, ("T","C"):1, ("G","A"):1}
	transversion = {("C","T"):0, ("A","G"):0, ("T","C"):0, ("G","A"):0}
	for nt1 in ["A","C","G","T"]:
		for nt2 in ["A","C","G","T"]:
			if (nt1, nt2) not in transition:
				transition[(nt1, nt2)] = 0
			if (nt1, nt2) not in transversion:
				transversion[(nt1, nt2)] = 1
	aa_data = {
		'A' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 1, "Hydrophobicity": 1.8},
		'C' : {"Charge": 0, "Polar Neutral": 1, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": 2.5},
		'D' : {"Charge": -1, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": -3.5},
		'E' : {"Charge": -1, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": -3.5},
		'F' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 1, "Large": 1, "Small": 0, "Hydrophobicity": 2.8},
		'G' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": -0.4},
		'H' : {"Charge": 1, "Polar Neutral": 0, "Aromatic": 0, "Large": 1, "Small": 0, "Hydrophobicity": -3.2},
		'I' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": 4.5},
		'K' : {"Charge": 1, "Polar Neutral": 0, "Aromatic": 0, "Large": 1, "Small": 0, "Hydrophobicity": -3.9},
		'L' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": 3.8},
		'M' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 1, "Small": 0, "Hydrophobicity": 1.9},
		'N' : {"Charge": 0, "Polar Neutral": 1, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": -3.5},
		'P' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": -1.6},
		'Q' : {"Charge": 0, "Polar Neutral": 1, "Aromatic": 0, "Large": 1, "Small": 0, "Hydrophobicity": -3.5},
		'R' : {"Charge": 1, "Polar Neutral": 0, "Aromatic": 0, "Large": 1, "Small": 0, "Hydrophobicity": -4.5},
		'S' : {"Charge": 0, "Polar Neutral": 1, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": -0.8},
		'T' : {"Charge": 0, "Polar Neutral": 1, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": -0.7},
		'V' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 0, "Large": 0, "Small": 0, "Hydrophobicity": 4.2},
		'W' : {"Charge": 0, "Polar Neutral": 0, "Aromatic": 1, "Large": 1, "Small": 0, "Hydrophobicity": -0.9},
		'Y' : {"Charge": 0, "Polar Neutral": 1, "Aromatic": 1, "Large": 1, "Small": 0, "Hydrophobicity": -1.3}	
	}
	aa_delta = {("Delta " + k):{} for k in aa_data[list(aa_data.keys())[0]]}
	for k1 in aa_data:
		for k2 in aa_data:
			for k3 in aa_delta:
				aa_delta[k3][(k1,k2)] = aa_data[k2][k3.replace("Delta ", "")] - aa_data[k1][k3.replace("Delta ", "")]

	ntseq = tools.read_fasta(ntfasta)[genename].replace("-", "")
	aaseq = tools.read_fasta(aafasta)[genename].replace("-", "")
	now = datetime.datetime.now()

	mutdf = pd.read_csv(mutdf, sep="\t", index_col=0)
	datadf = pd.read_csv(datadf, sep="\t", index_col=0)
	response = mutdf[list(respcols.keys())].to_dict()
	other_cols = [col for col in list(mutdf.columns) if col not in response and "c." not in col and "p." not in col]
	other_cols = mutdf[other_cols].to_dict()
	thrombosis = ("Thrombosis" in response.keys())

	if not dups:
		delset = []
		varsets = {}
		for i, row in mutdf.iterrows():
			varsets[i] = set()
			for j in list(mutdf.columns):
				if str(j).startswith("c.") and tools.is_numeric(row[j]) and float(row[j]) > 0:
					varsets[i] = varsets[i].union(set([j]))
		vlist = list(varsets.keys())
		for i, p1 in enumerate(vlist):
			for j in range(i):
				try:
					p2 = vlist[j]
					if varsets[p1] == varsets[p2] and mutdf[list(respcols.keys())[0]][p1] == mutdf[list(respcols.keys())[0]][p2]:
						delset.append(p1)
				except Exception as e:
					pass
		mutdf = mutdf.drop(list(set(delset)), axis=0)	

	ntvariants = {}
	aavariants = {}
	heterozygous = {}
	for i, row in mutdf.iterrows():
		nonsense = False
		heterozygous[i] = {}
		ntvariants[i] = {}
		aavariants[i] = {}
		for j in list(mutdf.columns):
			if j not in response and "c." in str(j) and mutdf[j][i] > 0:
				mutcode = str(j).replace("c.", "")
				posnt, nts = tools.get_mutation_data(mutcode)
				posaa, aas = tools.get_mutant_aa((posnt, nts), ntseq, aaseq=aaseq, index=1)
				if aas[1] == "*" or aas[1].lower() == "stop": #ignore nonsense mutations
					nonsense = True
				if aas[0] != aas[1]:
					aavariants[i][int(posaa)] = [aas[0], aas[1]]
				ntvariants[i][int(posnt)] = [nts[0], nts[1]]
				if hetero and mutdf[j][i] < 1:
					heterozygous[i].append(int(posnt))
		if nonsense: #remove patient with nonsense mutation
			del heterozygous[i]
			del ntvariants[i]
			del aavariants[i]
	ntvarset = sorted(set([(pos, (nts[0], nts[1])) for v in ntvariants.values() for pos, nts in v.items()]))
	aavarset = sorted(set([(pos, (aas[0], aas[1])) for v in aavariants.values() for pos, aas in v.items()]))

	aacols = {"AA conservation":{"invert": True, "pct": True}, "AA entropy":{"invert": False, "pct": False}, "AA variance":{"invert": False, "pct": False}, "BLOSUM conservation":{"invert": False, "pct": False}, "Rare codon enrichment":{"invert": False, "pct": False}, "Consurf (wo 3D info)":{"invert": True, "pct": False}, "Accessible surface area (PDBPisa)":{"invert": False, "pct": False}, "Solvation energy (PDBPisa)":{"invert": False, "pct": False}, "Relative surface area (NetsurfP2)":{"invert": False, "pct": False}}
	ntcols = {"NT conservation":{"invert": True, "pct": True}, "NT entropy":{"invert": False, "pct": False}, "NT variance":{"invert": False, "pct": False}}
	vepcols = {"Polyphen":{"invert": False, "pct": True},"SIFT":{"invert": False, "pct": True},"GnomAD":{"invert": False, "pct": True}}

	consurf3d = []
	for i, row in datadf.iterrows():
		if tools.is_numeric(row["Conservation (with PDB)"]) and not np.isnan(float(row["Conservation (with PDB)"])):
			consurf3d.append(float(row["Conservation (with PDB)"]))
		else:
			consurf3d.append(float(row["Consurf (wo 3D info)"]))

	aadist = compute_distribution(tools.read_fasta(aafasta), genename, alphabet="ACDEFGHIKLMNPQRSTVWXY")
	ntdist = compute_distribution(tools.read_fasta(ntfasta), genename, alphabet="TACG")

	scores = {}

	aapairings = {}

	vep_scores, ntpairings = get_vep(ntvariants, genename, get_pairs=True, epsilon=0.001, outdir=os.path.dirname(output))
	provean = get_provean(aavariants, genename, aaseq, sleeptime=10, timeout=60, outdir=os.path.dirname(output))
	#provean_nt = provean_nts(ntvariants, tools.read_fasta(ntfasta), genename, outdir=os.path.dirname(output), frac=True)

	if pairmodel.upper() == "LD":
		ntseq = tools.read_fasta(ntfasta)[genename].replace("-", "")
		for (mut1, mut2), val in ntpairings.items():
			pos1, nts1 = tools.get_mutation_data(mut1)
			pos2, nts2 = tools.get_mutation_data(mut2)
			(pos1, aas1) = tools.get_mutant_aa((pos1, nts1), ntseq, index=1)
			(pos2, aas2) = tools.get_mutant_aa((pos2, nts2), ntseq, index=1)
			if (pos1, aas1) in aavarset and (pos2, aas2) in aavarset:
				aapairings[(aas1[0]+str(pos1)+aas1[1],aas2[0]+str(pos2)+aas2[1])] = val

	if pairmodel.upper() in ["EV", "PLMC", "EC"]:
		if aapairing != "":
			with open(aapairing, "r") as inf:
				for i, line in enumerate(inf):
					if i > 0:
						try:
							pos1 = int(line.split()[0])
							pos2 = int(line.split()[1])
							aas1 = {k:v for (k,v) in aavarset}[pos1]
							aas2 = {k:v for (k,v) in aavarset}[pos2]
							aapairings[(aas1[0] + str(pos1) + aas1[1], aas2[0] + str(pos2) + aas2[1])] = float(line.split()[-1])
						except:
							pass
			u = statistics.mean([x for x in aapairings.values() if tools.is_numeric(x) and not np.isnan(x)])
			o = statistics.stdev([x for x in aapairings.values() if tools.is_numeric(x) and not np.isnan(x)])
			aapairings = {k:scipy.stats.norm((v-u)/o) for k,v in aapairings.items()}
		if ntpairing != "":
			with open(ntpairing, "r") as inf:
				for i, line in enumerate(inf):
					if i > 0:
						try:
							pos1 = int(line.split()[0])
							pos2 = int(line.split()[1])
							nts1 = {k:v for (k,v) in ntvarset}[pos1]
							nts2 = {k:v for (k,v) in ntvarset}[pos2]
							ntpairings[(nts1[0] + str(pos1) + nts1[1], nts2[0] + str(pos2) + nts2[1])] = float(line.split()[-1])
						except:
							pass
			u = statistics.mean([x for x in aapairings.values() if tools.is_numeric(x) and not np.isnan(x)])
			o = statistics.stdev([x for x in aapairings.values() if tools.is_numeric(x) and not np.isnan(x)])
			aapairings = {k:scipy.stats.norm((v-u)/o) for k,v in aapairings.items()}
	elif genename != "" and pairmodel.upper() == "MSA":
		aapairings = pairset_vs_min(tools.read_fasta(aafasta), genename, aavariants, space="-", index=1)
		ntpairings = pairset_vs_min(tools.read_fasta(ntfasta), genename, ntvariants, space="-", index=1)
	elif pairmodel.lower() in ["ind", "independent"]:
		aapairings = {}
		ntpairings = {}

	ntpairscores = {}
	aapairscores = {}
	for var, muts in ntvariants.items():
		ntpairscores[var] = {}
		for pos, nts in muts.items():
			ntpairscores[var][(pos, tuple(nts))] = 0
			for pos2, nt2 in muts.items():
				try:
					pair_ij = (min(nts[0]+str(pos)+nts[1], nts2[0]+str(pos2)+nts2[1], key=lambda k: int(k[1:-1])), max(nts[0]+str(pos)+nts[1], nts2[0]+str(pos2)+nts2[1], key=lambda k: int(k[1:-1])))
					ntpairscores[var][(pos, tuple(nts))] += ntpairings[pair_ij]
				except Exception as e:
					pass
	for var, muts in aavariants.items():
		aapairscores[var] = {}
		for pos, aas in muts.items():
			aapairscores[var][(pos, tuple(aas))] = 0
			for pos2, aas2 in muts.items():
				try:
					pair_ij = (min(aas[0]+str(pos)+aas[1], aas2[0]+str(pos2)+aas2[1], key=lambda k: int(k[1:-1])), max(aas[0]+str(pos)+aas[1], aas2[0]+str(pos2)+aas2[1], key=lambda k: int(k[1:-1])))
					aapairscores[var][(pos, tuple(aas))] += aapairings[pair_ij]
				except:
					pass

	if len(vep_scores) > 0:
		for var, muts in ntvariants.items():
			for col in vepcols:
				if col not in scores:
					scores[col] = {}
				scores[col][var] = score_variant(muts, vep_scores[col], ntpairscores[var], invert=vepcols[col]["invert"], pctflag=vepcols[col]["pct"], mode=combine_model, hetero=False, genename=genename, parameter=col)
	if len(provean) > 0:
		scores["PROVEAN"] = {}
		for var, muts in aavariants.items():
			scores["PROVEAN"][var] = score_variant(muts, provean, aapairscores[var], invert=True, pctflag=False, mode=combine_model, hetero=False, genename=genename, parameter="PROVEAN")
	'''
	if len(provean_nt) > 0:
		scores["PROVEAN (NT)"] = {}
		for var, muts in ntvariants.items():
			scores["PROVEAN (NT)"][var] = score_variant(muts, provean_nt, ntpairscores[var], invert=True, pctflag=True, mode=combine_model, hetero=False, genename=genename, parameter="PROVEAN (NT)")
	'''
	for col in aacols:
		scores[col] = {}
		for var, muts in aavariants.items():
			scores[col][var] = score_variant(muts, list(datadf[col]), aapairscores[var], invert=aacols[col]["invert"], pctflag=aacols[col]["pct"], mode=combine_model, hetero=False, genename=genename, parameter=col)
	for col in ntcols:
		scores[col] = {}
		for var, muts in ntvariants.items():
			scores[col][var] = score_variant(muts, list(datadf[col]), ntpairscores[var], invert=ntcols[col]["invert"], pctflag=ntcols[col]["pct"], mode=combine_model, hetero=False, genename=genename, parameter=col)

	scores["Conservation (with PDB)"] = {}
	for var, muts in aavariants.items():
		scores["Conservation (with PDB)"][var] = score_variant(muts, consurf3d, aapairscores[var], invert=True, pctflag=False, mode=combine_model, hetero=False, genename=genename, parameter="Conservation (with PDB)")

	for (name, d) in [("Transition", transition), ("Transversion", transversion)]:
		scores[name] = {}
		for var, muts in ntvariants.items():
			scores[name][var] = score_variant(muts, d, ntpairscores[var], mode=combine_model, pctflag=False, hetero=False, invert=False, genename=genename, parameter=name)
	for name in aa_delta.keys():
		scores[name] = {}
		for var, muts in aavariants.items():
			scores[name][var] = score_variant(muts, aa_delta[name], aapairscores[var], mode=combine_model, pctflag=False, hetero=False, invert=False, genename=genename, parameter=name)
	scores["MSA loglikelihood (AA)"] = {}
	for var, muts in aavariants.items():
		scores["MSA loglikelihood (AA)"][var] = score_variant(muts, aadist, aapairscores[var], mode=combine_model, pctflag=True, hetero=False, invert=False, genename=genename, parameter="MSA loglikelihood (AA)")
	scores["MSA loglikelihood (NT)"] = {}
	for var, muts in ntvariants.items():
		scores["MSA loglikelihood (NT)"][var] = score_variant(muts, ntdist, ntpairscores[var], mode=combine_model, pctflag=True, hetero=False, invert=False, genename=genename, parameter="MSA loglikelihood (NT)")

	evdf = datadf[[col for col in datadf.columns if "EVmutation " in col]]
	scores["EVmutation (AA)"] = {}
	for var, muts in aavariants.items():
		scores["EVmutation (AA)"][var] = score_variant(muts, evdf, aapairscores[var], mode=combine_model, pctflag=False, hetero=False, invert=False, genename=genename, parameter="EVmutation (AA)")

	evdf = datadf[[col for col in datadf.columns if "EVmutation(nt)" in col]]
	scores["EVmutation (NT)"] = {}
	for var, muts in ntvariants.items():
		scores["EVmutation (NT)"][var] = score_variant(muts, evdf, ntpairscores[var], mode=combine_model, pctflag=False, hetero=False, invert=False, genename=genename, parameter="EVmutation (NT)")
	'''
	scores["Pairing AA"] = {}
	for var, muts in aavariants.items():
		tmpscores = {}
		for pos, aas in muts.items():
			tmpscores[(pos, tuple(aas))] = aapairscores[var][(pos, tuple(aas))]
		tmpscores = {mut:score for mut,score in tmpscores.items() if not np.isnan(score)}
		if combine_model.lower() in ["additive", "add"]:
			scores["Pairing AA"][var] = (sum(tmpscores.values()))
		elif combine_model.lower() in ["max", "maximum"]:
			scores["Pairing AA"][var] = (max(tmpscores.values()))
		elif combine_model.lower() in ["maxminmedian", "minmaxmedian", "mmm"]:
			scores["Pairing AA"][var] = (max(tmpscores.values()) + min(tmpscores.values()) - statistics.mean(tmpscores.values()))

	scores["Pairing NT"] = {}
	for var, muts in ntvariants.items():
		tmpscores = {}
		for pos, nts in muts.items():
			tmpscores[(pos, tuple(nts))] = ntpairscores[var][(pos, tuple(nts))]
		tmpscores = {mut:score for mut,score in tmpscores.items() if not np.isnan(score)}
		if combine_model.lower() in ["additive", "add"]:
			scores["Pairing NT"][var] = (sum(tmpscores.values()))
		elif combine_model.lower() in ["max", "maximum"]:
			scores["Pairing NT"][var] = (max(tmpscores.values()))
		elif combine_model.lower() in ["maxminmedian", "minmaxmedian", "mmm"]:
			scores["Pairing NT"][var] = (max(tmpscores.values()) + min(tmpscores.values()) - statistics.mean(tmpscores.values()))
	'''
	scores["Number of mutants"] = {}
	for var, muts in ntvariants.items():
		scores["Number of mutants"][var] = 0
		for pos, nts in muts.items():
			if hetero and pos in heterozygous[var]:
				mult=0.5
			else:
				mult=1
			scores["Number of mutants"][var] += mult

	for (pos, nts) in ntvarset:
		scores[nts[0] + str(pos) + nts[1]] = {}
		for var, muts in ntvariants.items():
			if pos in muts.keys() and muts[pos][0] == nts[0] and muts[pos][1] == nts[1]:
				if pos in heterozygous[var] and hetero:
					mult = 0.5
				else:
					mult = 1
				scores[nts[0] + str(pos) + nts[1]][var] = mult
			else:
				scores[nts[0] + str(pos) + nts[1]][var] = 0

	scores.update(other_cols)
	scores.update(response)

	scores_df = pd.DataFrame(scores)
	if multonly:
		scores_df = scores_df.loc[scores_df["Number of mutants"] >= 2]
	scores_df.to_csv(os.path.splitext(output)[0] + "_" + combine_model + pairmodel + "_" + str(now.year) + str(now.month).rjust(2,'0') + str(now.day).rjust(2,'0') + ".tsv", sep="\t")

	analyze_data(datafile=output, response=respcols, pairmodel=pairmodel, combine_model=combine_model, now=now, multonly=multonly, sep="\t")

	plot_all(os.path.splitext(output)[0] + "_" + combine_model + pairmodel + "_" + str(now.year) + str(now.month).rjust(2,'0') + str(now.day).rjust(2,'0') + ".tsv", response=response)
	return(scores_df)

def plot_all(df, response=["Activity", "Expression"], excluded_rows=["WT"], prefix="", sheetname=None):
	if prefix == None or len(prefix) < 1:
		try:
			prefix = os.path.basename(df).split("_")[1]
		except:
			prefix = ""
	directory = os.path.dirname(df)
	if sheetname == None or len(sheetname) < 1:
		df = pd.read_csv(df, sep="\t", index_col=0)
	else:
		df = pd.read_excel(df, sheet_name=sheetname, index_col=0)
	try:
		df.drop(excluded_rows, axis=0)
	except KeyError as e:
		pass
	
	figcount = 1
	for col1 in df.columns:
		if not re.search("[ACGT]\d+[ACGT]", str(col1)):
			if col1 not in response:
				for col2 in response:
					try:
						fig = plt.figure(figcount)
						ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
						for i, (row, pholder) in enumerate(df.iterrows()):
							if row != "WT":
								r,g,b = tools.hsv_to_rgb(float(i)/len(df.index), 1, 1)
								ax.scatter(float(df[col1][row]), float(df[col2][row]), marker=list(markers.keys())[i%len(markers)], c=[[r,g,b]], label=row)
						ax.set_xlabel(col1)
						ax.set_ylabel(col2)
						ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
						fig.savefig(os.path.join(directory, col1 + " " + col2 + "_" + prefix + ".png"), dpi=300)
						plt.cla()
						plt.clf()
						figcount += 1
					except ValueError as e:
						pass

#summarize the variants and patients of different datasets
def summary_stats(geneset, multonly=True, seqdir="./gene_data/"):
	data = {}
	for gene, directory, f in geneset:
		df = pd.read_csv(f, sep="\t", index_col=0).T.to_dict()
		if multonly:
			df = {k:v for k,v in df.items() if sum([v[x] for x in v.keys() if x.startswith("c.") and tools.is_numeric(v[x]) and not np.isnan(v[x])]) >= 2}
		data[f] = {"Number of patients":len([k for k in df if df[k]["Disease"] > 0]), "Number of controls":len([k for k in df if df[k]["Disease"] <= 0]), "Number of individual variants":len([x for x in pd.DataFrame(df).T.columns if x.startswith("c.") and sum([df[k][x] for k in df.keys() if tools.is_numeric(df[k][x]) and not np.isnan(df[k][x])]) > 0])}
		haplos = combine_documents.unique_haplo(f, binary="Disease", multonly=multonly)
		data[f].update({"Number of unique haplotypes/combinations (disease)":len(haplos[1]), "Number of unique haplotypes/combinations (neutral)":len(haplos[0])})
	for k,v in data.items():
		df = pd.read_csv(k, sep="\t", index_col=0).T.to_dict()
		if multonly:
			df = {k:v for k,v in df.items() if sum([v[x] for x in v.keys() if x.startswith("c.") and tools.is_numeric(v[x]) and not np.isnan(v[x])]) >= 2}
		try:
			if "LOVD" in k:
				ntseq = tools.read_fasta(seqdir + os.path.basename(k).split("_")[1].split(".")[0] + ".fasta")["ORF"].replace("-", "")
			else:
				ntseq = tools.read_fasta(seqdir + os.path.basename(k).split("_")[0] + ".fasta")["ORF"].replace("-", "")
		except IOError as e:
			if "Plaimauer" in k:
				ntseq = tools.read_fasta(seqdir + "ADAMTS13" + ".fasta")["ORF"].replace("-", "")
			else:
				ntseq = ""
		neutral = []
		disease = []
		for k2, v2 in df.items():
			for k3, v3 in v2.items():
				if k3.startswith("c.") and tools.is_numeric(v3) and not np.isnan(v3) and float(v3) > 0:
					if "Gene" in v2.keys():
						k3 = v2["Gene"] + ":" + k3
					if v2["Disease"] <= 0:
						neutral.append(k3)
					else:
						disease.append(k3)
		disease = set(disease)
		neutral = set(neutral)
		data[k].update({"Disease-only variants": len(set([x for x in disease if x not in neutral])), "Neutral-only variants":len(set([x for x in neutral if x not in disease])), "Number of variants (disease)":len(set(disease)), "Number of variants (neutral)": len(set(neutral)), "Number of individual variants seen in both disease & control":len([x for x in neutral if x in disease])})
		try:
			missense = {}
			for k2,l in [("Neutral", neutral), ("Disease", disease)]:
				print(k2 + "\t" + str(l))
				missense[k2] = 0
				for variant in set(l):
					if not variant.startswith("c.") and ":c." in variant:
						g = variant.split(":c.")[0]
						ntseq = tools.read_fasta(seqdir + g + ".fasta")["ORF"].replace("-", "")
					aavar = tools.get_mutant_aa(tools.get_mutation_data(variant.split("c.")[-1]), ntseq=ntseq, index=1)
					if aavar[1][0] != aavar[1][1]:
						missense[k2] += 1
			data[k].update({"Missense variants (disease)":missense["Disease"], "Missense variants (neutral)":missense["Neutral"]})
		except Exception as e:
			traceback.print_exc()
			data[k].update({"Missense variants (disease)":"", "Missense variants (neutral)":""})
	newdata = {}
	for k,v in data.items():
		if "LOVD" in k:
			name = ("LOVD", k.split("_")[-1].split(".")[0], os.path.basename(k).split("_")[0])
		else:
			name = ("Internal", str(os.path.basename(k).split("_")[0])) 
		if name == ("Internal", "VWF"):
			name = ("Goodeve et. al., 2007", "VWF")
		elif name == ("Internal", "ADAMTS13"):
			name = ("Hing et. al. (Internal)", "ADAMTS13")
		elif name == ("Internal", "Plaimauer"):
			name = ("Plaimauer et. al.", "ADAMTS13")
		newdata[name] = data[k]
	df = pd.DataFrame(newdata).T
	df = df[["Number of individual variants", "Number of individual variants seen in both disease & control", "Number of patients", "Number of unique haplotypes/combinations (disease)", "Number of variants (disease)", "Missense variants (disease)", "Disease-only variants", "Number of controls", "Number of unique haplotypes/combinations (neutral)", "Number of variants (neutral)", "Missense variants (neutral)", "Neutral-only variants"]]
	return(df)

#summarize for top 4 predictors, which are adjustable
def top_4_sum(top4=["Polyphen","SIFT","PROVEAN", "MSA loglikelihood (AA)"]):
	data = {}
	for directory in datasets:
		try:
			df = pd.read_excel(os.path.join(directory, "Disease_features_summarized.xlsx"), sheet_name="features by AUC", index_col=0).T.to_dict()
			for feat in top4:
				if feat not in data:
					data[feat] = []
				try:
					data[feat] += list(df[feat].values())
				except:
					pass
		except Exception as e:
			print(e)
	for i, key in enumerate(data.keys()):
		plt.scatter([i for j in range(len(data[key]))], data[key], c=[list(tools.hsv_to_rgb(random.random(), 1, 1))], alpha=0.05, s=0)
	plt.ylabel("AUC")

	plt.ylim([0,1.05])
	plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
	plt.title("Distribution of AUC values across all datasets and models")
	plt.xticks(np.arange(len(data.keys())), list(data.keys()))
	plt.savefig("top4_dist.tif", dpi=500)
	return(data)

def compute_ROC_bound(x, y, bound):
	tp = len([i for i in range(len(x)) if x[i] >= bound and y[i] == 1])
	fp = len([i for i in range(len(x)) if x[i] >= bound and y[i] != 1])	
	return(tp, fp)

def compare_predictions_bound(x0, x1, y, bound, names, patients=[]):
	comparison = []
	if patients == None or len(patients) < 1:
		patients = [i for i in range(len(y))]
	bound = [min(x0)+bound*(max(x0)-min(x0)), min(x1)+bound*(max(x1)-min(x1))]
	for i in range(len(y)):
		if (x0[i] < bound[0] and x1[i] > bound[1] and y[i] == 0) or (x0[i] > bound[0] and x1[i] < bound[1] and y[i] == 1):
			comparison.append([patients[i], [names[0], x0[i], bound[0]], [names[1], x1[i], bound[1]], y[i]])
		elif (x0[i] > bound[0] and x1[i] < bound[1] and y[i] == 0) or (x1[i] < bound[0] and x1[i] > bound[1] and y[i] == 1):
			comparison.append([patients[i], [names[1], x1[i], bound[1]], [names[0], x0[i], bound[0]], y[i]])
	return(comparison)

def compare_models_AUC(inf, ycol="Disease", step=0.0001):
	wb = openpyxl.load_workbook(inf, read_only=True, keep_links=False)
	sheetnames = list(wb.sheetnames)
	columns = list(pd.read_excel(inf, sheet_name=sheetnames[0], index_col=0).columns)
	predict_comparison = {}
	init_time = time.time()
	limits = {}
	data = {}
	for count, xcol in enumerate(columns):
		limits[xcol] = {}
		predict_comparison[xcol] = []
		try:
			data[xcol] = {}
			patients = set([])
			for sheetname in sheetnames:
				df = pd.read_excel(inf, sheet_name=sheetname, index_col=0).T.to_dict()
				data[xcol][sheetname] = {k:{xcol:v[xcol], ycol:v[ycol]} for k,v in df.items()}
				data[xcol][sheetname] = {k:v for k,v in data[xcol][sheetname].items() if tools.is_numeric(v[xcol]) and not np.isnan(v[xcol]) and tools.is_numeric(v[ycol]) and not np.isnan(v[ycol])}
				if len(patients) < 1:
					patients = set(data[xcol][sheetname].keys())
				else:
					patients = patients.intersection(data[xcol][sheetname].keys())
			patients = list(patients)
			flip = False
			if all([True for sheetname in data[xcol] if sklearn.metrics.roc_auc_score([data[xcol][sheetname][xi][ycol] for xi in patients],[data[xcol][sheetname][xi][xcol] for xi in patients]) < sklearn.metrics.roc_auc_score([data[xcol][sheetname][xi][ycol] for xi in patients],[-data[xcol][sheetname][xi][xcol] for xi in patients])]):
				flip = True
			#print(xcol + "\t" + str(len(patients)) + "\t" + str(flip))
			truey = [data[xcol][sheetname][xi][ycol] for xi in patients]
			optbas = [0.0, -np.inf]
			for i in np.arange(0, 1, step):
				for m, sheet1 in enumerate(sheetnames):
					if flip:
						x1 = [-data[xcol][sheet1][xi][xcol] for xi in patients]
					else:
						x1 = [data[xcol][sheet1][xi][xcol] for xi in patients]
					bas = sklearn.metrics.balanced_accuracy_score([data[xcol][sheetname][xi][ycol] for xi in patients], [1 if x1[j] > i else 0 for j in range(len(x1))])
					
					limits[xcol][sheet1] = {"min":min(x1), "max":max(x1), "opt_bound":optbas[1]}
					if bas > optbas[0]:
						optbas = [bas, i]
					for n in range(m):
						sheet2 = sheetnames[n]
						if flip:
							x2 = [-data[xcol][sheet2][xi][xcol] for xi in patients]
						else:
							x2 = [data[xcol][sheet2][xi][xcol] for xi in patients]
						predict_comparison[xcol] += compare_predictions_bound(x1, x2, [data[xcol][sheetname][xi][ycol] for xi in patients], i, names=[sheet1, sheet2], patients=patients)
		except Exception as e:
			print(xcol + "\t" + str(e))
		tools.update_time(count, len(columns), init_time)
	new_predict_comparison = {}
	init_time = time.time()
	for i, xcol in enumerate(predict_comparison):
		new_predict_comparison[xcol] = {}
		for xi in predict_comparison[xcol]:
			if xi[0] not in new_predict_comparison[xcol]:
				new_predict_comparison[xcol][xi[0]] = {"Add better than max": 0, "Max better than add": 0, "ld better than ind": 0, "ind better than ld": 0, "Disease":xi[-1]}
			if (("addind" == xi[1][0] and "maxind" == xi[2][0]) or ("addld" == xi[1][0] and "maxld" == xi[2][0])):
				new_predict_comparison[xcol][xi[0]]["Add better than max"] += 1
			elif (("maxind" == xi[1][0] and "addind" == xi[2][0]) or ("maxld" == xi[1][0] and "addld" == xi[2][0])):
				new_predict_comparison[xcol][xi[0]]["Max better than add"] += 1
			elif (("maxld" == xi[1][0] and "maxind" == xi[2][0]) or ("addld" == xi[1][0] and "addind" == xi[2][0])):
				new_predict_comparison[xcol][xi[0]]["ld better than ind"] += 1
			elif (("maxind" == xi[1][0] and "maxld" == xi[2][0]) or ("addind" == xi[1][0] and "addld" == xi[2][0])):
				new_predict_comparison[xcol][xi[0]]["ind better than ld"] += 1
		tools.update_time(i, len(predict_comparison), init_time)
	writer = pd.ExcelWriter("./model_variant_comparison_" + str(1/i) + ".xlsx", engine="openpyxl")
	for xcol in new_predict_comparison.keys():
		pd.DataFrame(new_predict_comparison[xcol]).T.to_excel(writer, xcol)	
	writer.save()
	writer.close()
	
	fig, ax = plt.subplots()
	ax2 = ax.twinx()
	for i, model in enumerate(sheetnames):
		y0 = [-data["SIFT"][model][k]["SIFT"] for k in data["SIFT"][model] if data["SIFT"][model][k][ycol] == 0]
		y1 = [-data["SIFT"][model][k]["SIFT"] for k in data["SIFT"][model] if data["SIFT"][model][k][ycol] == 1]
		if "max" in str(model):
			ax.scatter([i for j in range(len(y0))], y0, c="blue", alpha=0.1)
			ax.scatter([i for j in range(len(y1))], y1, c="red", alpha=0.1)
			ax.scatter([i], [limits["SIFT"][model]["opt_bound"]], c="black")
		else:
			ax2.scatter([i for j in range(len(y0))], y0, c="blue", alpha=0.1)
			ax2.scatter([i for j in range(len(y1))], y1, c="red", alpha=0.1)
			ax2.scatter([i], [limits["SIFT"][model]["opt_bound"]], c="black")
	plt.xticks([i for i in range(len(sheetnames))], [sheetnames[i] for i in range(len(sheetnames))])
	plt.show()
	
	return(predict_comparison, new_predict_comparison, limits)

#runs all scoring and analyses
#requires as input the gene, the directory to output to, and the spreadsheet containing the combinations
def pipeline(genes, dataloc = "./gene_data/", response="Disease", combine_row=None, dups=False, normal=0):
	try:
		if isinstance(genes, dict) or len(genes[0]) < 3:
			raise Exception("Input is not a list of triples")
	except:
		genes = genes.items()
	now = datetime.datetime.now()
	if combine_row == None or len(combine_row) < 1:
		combine_row = response + " Mann-Whitney p"
	#generate MSAs
	p_aa = {}
	p_nt = {}
	for genename, pholder, mutdf in genes:
		directory = os.path.join(os.path.dirname(mutdf))
		if not os.path.exists(os.path.join(directory, genename)):
			os.makedirs(os.path.join(directory, genename))
		fasta = dataloc + genename + ".fasta"
		ntseq = tools.read_fasta(fasta)["ORF"].replace("-", "")
		aaseq = compute_features.translate_seq(ntseq)
		if not os.path.exists(os.path.join(directory, genename, "aa_msa.fasta")) or len(tools.read_fasta(os.path.join(directory, genename, "aa_msa.fasta"))) < 2: #AA MSA not found
			p_aa[genename] = Process(target=compute_features.run_blast, args=(ntseq, genename, "aa", "nr", True, None, os.path.join(directory, genename), "-"))
			p_aa[genename].start()
		if not os.path.exists(os.path.join(directory, genename, "nt_msa.fasta")) or len(tools.read_fasta(os.path.join(directory, genename, "nt_msa.fasta"))) < 2:	#NT MSA not found 						
			p_nt[genename] = Process(target=compute_features.run_blast, args=(ntseq, genename, "nt", "nt", True, None, os.path.join(directory, genename), "-"))
			p_nt[genename].start()
	for genename, pholder, mutdf in genes:
		try:
			directory = os.path.join(os.path.dirname(mutdf))
			datadf = os.path.join(dataloc, genename + ".tsv")
			fasta = os.path.join(dataloc, genename + ".fasta")
			ntseq = tools.read_fasta(fasta)["ORF"].replace("-", "")
			aaseq = compute_features.translate_seq(ntseq)
			if genename in p_aa:
				p_aa[genename].join()
			if genename in p_nt:
				p_nt[genename].join()
			ntseqs = tools.read_fasta(os.path.join(directory, genename, "nt_msa.fasta"))
			aaseqs = tools.read_fasta(os.path.join(directory, genename, "aa_msa.fasta"))
			if any([True for seq in ntseqs.values() if len(seq) != len(ntseqs[genename])]):
				ntseqs = tools.align_sequences(ntseqs)
			if any([True for seq in aaseqs.values() if len(seq) != len(aaseqs[genename])]):
				aaseqs = tools.align_sequences(aaseqs)
			tools.write_fasta(aaseqs, os.path.join(directory, genename, "aa_msa.fasta"))
			tools.write_fasta(ntseqs, os.path.join(directory, genename, "nt_msa.fasta"))
			uniprot = compute_features.parse_uniprot(genename, ntseq)
			pdbdat = compute_features.parse_pdb(uniprot["Accession"], genename, ORF=ntseq)
			pdb = [pdbdat[0][0].upper(), pdbdat[0][1][0]]
			if not any(fname.endswith('.pdb') for fname in os.listdir(os.path.join(directory, genename))):
				r = requests.get("https://files.rcsb.org/download/"+pdb[0]+".pdb")
				with open(os.path.join(directory, genename, pdb[0] + ".pdb"), "w") as pdbfile:
					pdbfile.write(r.text)
			for pairmodel in ["ind", "ld"]:
				for combinemodel in ["add", "max"]:
					score_variants(mutdf, datadf, os.path.join(directory, genename, "aa_msa.fasta"), os.path.join(directory, genename, "nt_msa.fasta"), pairmodel=pairmodel, combine_model=combinemodel, genename=genename, pdb = os.path.join(directory, genename, pdb[0] + ".pdb"), chain = pdb[1], output=os.path.join(directory, genename, response+".tsv"), respcols = {response:normal}, dups=dups, multonly=True)
			combine_documents.merge_files(os.path.join(directory, genename), response+"\_analyzed\_([a-zA-Z]+)\_"+str(now.year) + str(now.month).rjust(2,'0') + str(now.day).rjust(2,'0')+"\.tsv")
			combine_documents.merge_files(os.path.join(directory, genename), response+"\_([a-zA-Z]+)\_"+str(now.year) + str(now.month).rjust(2,'0') + str(now.day).rjust(2,'0')+"\.tsv")
			
			try:
				regress_scores(os.path.join(directory, genename, response.lower()+"_a-za-z_"+str(now.year) + str(now.month).rjust(2,'0') + str(now.day).rjust(2,'0')+"tsv.xlsx"), response=response, normal=normal, genename=genename)
			except Exception as e:
				print(e)
			datamw, scoresmw, sheetscoresmw, bestsheetscoresmw = combine_documents.rank_rows(os.path.join(directory, genename, response.lower()+"_analyzed_a-za-z_"+str(now.year) + str(now.month).rjust(2,'0') + str(now.day).rjust(2,'0')+"tsv.xlsx"), [combine_row])
			datamw = datamw[combine_row]
			datamw = {k:{k1:v1 for (k1, v1) in v} for k,v in datamw.items()}
			scoresmw = {k:{"Median p-value":v} for k,v in scoresmw.items()}
			sheetscoresmw = {k:{"Median p-value":v} for k,v in sheetscoresmw.items()}
			bestsheetscoresmw = {k:{"Median p-value":v} for k,v in bestsheetscoresmw.items()}
			try:
				dataauc, scoresauc, sheetscoresauc, bestsheetscoresauc = combine_documents.rank_rows(os.path.join(directory, genename, response.lower()+"_analyzed_a-za-z_"+str(now.year) + str(now.month).rjust(2,'0') + str(now.day).rjust(2,'0')+"tsv.xlsx"), [response+" AUC"])
				dataauc = dataauc[response+" AUC"]
				dataauc = {k:{k1:v1 for (k1, v1) in v} for k,v in dataauc.items()}
				scoresauc = {k:{"Median AUC":v} for k,v in scoresauc.items()}
				sheetscoresauc = {k:{"Median AUC":v} for k,v in sheetscoresauc.items()}
				bestsheetscoresauc = {k:{"Median AUC":v} for k,v in bestsheetscoresauc.items()}
			except Exception as e:
				print(e)
				dataauc = {}
				scoresauc = {}
				sheetscoresauc = {}
				bestsheetscoresauc = {}
			newexcel = os.path.join(os.path.dirname(mutdf), genename, response + "_features_summarized.xlsx")
			writer = pd.ExcelWriter(newexcel, engine="openpyxl")
			pd.DataFrame(datamw).to_excel(writer, "features by p-value")
			pd.DataFrame(scoresmw).T.to_excel(writer, "overall p-values")
			pd.DataFrame(sheetscoresmw).to_excel(writer, "models by p-value")
			pd.DataFrame(bestsheetscoresmw).to_excel(writer, "models (best) by p-value")
			pd.DataFrame(dataauc).to_excel(writer, "features by AUC")
			pd.DataFrame(scoresauc).to_excel(writer, "overall AUCs")
			pd.DataFrame(sheetscoresauc).to_excel(writer, "models by AUC")	
			pd.DataFrame(bestsheetscoresauc).to_excel(writer, "models (best) by AUC")
			writer.save()
			writer.close()

			combine_documents.haplo_patients(mutdf, genename, response=response, binary=response)

		except FileNotFoundError as e:
			print("\033[93m" + genename + "\t" + str(e) + "\n" + str(traceback.format_exc()) + "\t" + "\033[0m")
		except Exception as e:
			print("\033[93m" + genename + "\t" + str(e) + "\n" + str(traceback.format_exc()) + "\t" + "\033[0m")
	
	genes_2 = []
	for (genename, pholder, f) in genes:
		genes_2.append((genename, pholder, f))
		wb = openpyxl.load_workbook(os.path.join(os.path.dirname(f), pholder, response + "_features_summarized.xlsx"), read_only=True, keep_links=False)
		sheetnames = list(wb.sheetnames)
		for sheetname in sheetnames:
			df = pd.read_excel(os.path.join(os.path.dirname(f), pholder, response + "_features_summarized.xlsx"), sheet_name=sheetname, index_col=0)
			if len(df.index) <= 0:
				genes_2.remove((genename, pholder, f))
				break
		
	combine_documents.merge_datasets([os.path.join(os.path.dirname(f), pholder) for (genename, pholder, f) in genes_2], response=response)
	merge_computed(genes, infile=(response.lower() + "_a-za-z_" + str(now.year) + str(now.month).rjust(2,'0') + str(now.day).rjust(2,'0') + "tsv.xlsx"))
	for pairmodel in ["ind", "ld"]:
		for combine_model in ["add", "max"]:
			analyze_data("merged_outputs.xlsx", pairmodel=pairmodel, combine_model=combine_model, response={response: normal}, multonly=True)
	combine_documents.merge_files("./", "merged\_outputs\_analyzed\_([a-zA-Z]+)\_" + str(now.year) + str(now.month).rjust(2,'0') + str(now.day).rjust(2,'0') + "\.tsv", sep="\t")
	merge_inputs(genes)
	remove_dups("./merged_outputs.xlsx", "./merged_inputs.tsv")
	remove_dups("./merged_inputs.tsv", "./merged_inputs.tsv", "./merged_inputs_nodups.tsv")
	summary = summary_stats(genes, multonly=True)
	summary.to_csv("./datasets_summary.tsv", sep="\t")
	variant_scores(seqdir=dataloc)

def analyze_indiv_vars(df1="./merged_inputs.tsv", df2="./variant_scores_(processed).tsv"):
	data = {}
	df1 = pd.read_csv(df1, sep="\t", index_col=0).T.to_dict()
	df2 = pd.read_csv(df2, sep="\t", index_col=0).T.to_dict()
	dfadd = pd.read_excel("merged_outputs.xlsx", sheet_name="addind", index_col=0).T.to_dict()
	dfmax = pd.read_excel("merged_outputs.xlsx", sheet_name="maxind", index_col=0).T.to_dict()
	for k,v in df1.items():
		data[k] = {}
		df1[k]["Variants"] = eval(df1[k]["Variants"])
		for xcol, inv in invert.items():
			try:
				if inv:
					maxvar = min(df1[k]["Variants"], key=lambda kv: df2[kv.replace("c.","")][xcol])
				else:
					maxvar = max(df1[k]["Variants"], key=lambda kv: df2[kv.replace("c.","")][xcol])
				dif = df2[maxvar.replace("c.", "")][xcol] - statistics.mean([df2[var.replace("c.", "")][xcol] for var in df1[k]["Variants"] if var != maxvar])
				data[k][xcol + " max variant"] = maxvar
				data[k][xcol + " max value"] = df2[maxvar.replace("c.", "")][xcol]
				data[k][xcol + " max distance to others"] = dif
				data[k][xcol + " add/max difference"] = dfadd[k][xcol] - dfmax[k][xcol]
			except KeyError as e:
				pass
			except Exception as e:
				print(e)
	df = pd.DataFrame(data).T
	corr = {}
	for xcol in invert.keys():
		try:
			l1, l2 = stat_tools.process_lists_numerical(list(df[xcol + " max distance to others"]), list(df[xcol + " add/max difference"]))
			corr[xcol] = scipy.stats.spearmanr(l1, l2)
		except KeyError as e:
			pass
	return(data, corr)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='This script reads gene data and creates a visualization of the NT and AA/codon data separately.', add_help=True)
	parser.add_argument('mutdf', type=str, metavar='PATH', help="path to mutation spreadsheet")
	parser.add_argument('datadf', type=str, metavar='PATH', help="path to spreadsheet of characteristics of gene")
	parser.add_argument('aafasta', type=str, metavar='PATH', help="path to amino acid MSA fasta file")
	parser.add_argument('ntfasta', type=str, metavar='PATH', default="", help="path to nucleotide MSA fasta file")
	parser.add_argument('--pairmodel', choices=["ind", "PLMC", "EV", "MSA", "LD"], default="ind", help="type of pairing model to combine effects of variants")
	parser.add_argument('--combinemodel', choices=["add", "additive", "minmaxmedian", "maxminmedian", "max", "maximum"], default="add", help="method to combine SNV scores")
	parser.add_argument('--genename', type=str, default="", help="name of gene in MSA file")
	parser.add_argument('--pdb', type=str, default="", help="location of pdb to use for crystal structure")
	parser.add_argument('--chain', type=str, default="A", help="corresponding chain in the pdb")
	parser.add_argument('--aapairing', type=str, metavar='FILE', default="", help="location of amino acid pairing file")
	parser.add_argument('--ntpairing', type=str, metavar='FILE', default="", help="location of nucleotide pairing file")
	parser.add_argument('--output', type=str, metavar='PATH', help="location of output file")
	parser.add_argument('--hetero', action="store_true", help="whether to consider heterozygosity")
	parser.add_argument('--response', type=str, default=["Activity", "Expression"], nargs="+", help="response variables")
	parser.add_argument('--multonly', action="store_true", help="whether to include only multiple mutants")
	clargs = parser.parse_args()

	score_variants(clargs.mutdf, clargs.datadf, clargs.aafasta, clargs.ntfasta, pairmodel=clargs.pairmodel, combine_model=clargs.combinemodel, genename=clargs.genename, pdb=clargs.pdb, chain=clargs.chain, aapairing=clargs.aapairing, ntpairing=clargs.ntpairing, output=clargs.output, hetero=clargs.hetero, multonly=clargs.multonly, response=clargs.response)

