import os, re, sys, time, os.path, argparse, math, statistics, pickle, subprocess, traceback, datetime, requests, warnings, openpyxl, bio_tools, copy, tools, stat_tools
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

#Input: string indicating clinical significance of variant
#Output: float variable in range 0-1 indicating pathogenicity of variant
def sigdict(k):
	if k in sigdict:
		return(sigdict[k])
	else:
		print("Warning: " + str(k) + " not recognized")
		return(0.5)

#For use in plotting
markers = {'o': 'circle', 'v': 'triangle_down', '^': 'triangle_up', 's': 'square', '*': 'star'}

#Names of (gene, directory, and file)
geneset = [('F8', 'F8', './LOVD/00495_F8.tsv'), ('F9', 'F9', './LOVD/02237_F9.tsv'), ('DMD', 'DMD', './LOVD/05324_DMD.tsv'), ('G6PD', 'G6PD', './LOVD/00308_G6PD.tsv'), ('HBB', 'HBB', './LOVD/03506_HBB.tsv'), ('VWF', 'VWF/VWD1', './LOVD/00498_VWF.tsv'), ('VWF', 'VWF/VWD2', './LOVD/03361_VWF.tsv'), ('VWF', 'VWF/VWD3', './LOVD/02126_VWF.tsv'), ('F9', 'F9', './F8_F9/F9_dataset_updated.tsv'), ('VWF', 'VWF', './VWD/VWF_dataset_updated.tsv'), ('ADAMTS13', 'ADAMTS13', './cTTP/ADAMTS13_dataset_updated.tsv'), ("ADAMTS13", "ADAMTS13", "./in_vitro_ADAMTS13/Plaimauer_variants.tsv")]

#Locations of datasets to use
datasets = ["./LOVD/F8/", "./LOVD/F9/", "./LOVD/DMD/", "./LOVD/HBB/", "./LOVD/VWF/VWD1/", "./LOVD/VWF/VWD2/", "./LOVD/VWF/VWD3/", "./LOVD/G6PD/", "./F8_F9/F8/", "./F8_F9/F9/", "./cTTP/ADAMTS13/", "./VWD/VWF/"]

#Indicates whether or not to "invert" a score for purposes of pathogenicity prediction
invert = {"Polyphen": True, "SIFT": True, "PROVEAN": True, "PROVEAN (NT)":True, "GnomAD": True, "Conservation (with PDB)": True, "AA percent identity": True, "AA entropy": False, "AA variance": False, "BLOSUM conservation": False, "Rare codon enrichment": True, "Consurf (wo 3D info)": True, "Accessible surface area (PDBPisa)": False, "Solvation energy (PDBPisa)": False, "NT percent identity": True, "NT entropy": False, "NT variance": False, "Relative surface area (NetsurfP2)": False, "MSA loglikelihood (AA)": True, "MSA loglikelihood (NT)": True, "EVmutation (AA)": True, "EVmutation (NT)": True, "Number of mutants": False}

#Input: multiple sequence alignment of sequences, name of focus sequence
#Output: distribution of characters in the MSA at each position in the focus sequence
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

#Input: Genetic variants in ORF coordinates, gene name
#Output: all data from VEP, including Polyphen and SIFT scores, MAFs, and linkage disequilibrium
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

#Input: Genetic variants in ORF coordinates, gene name, amino acid sequence
#Output: all data from PROVEAN
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
	#wait for Provean to finish, repeatedly query and check if success
	while counter * sleeptime < timeout*60 and not success:
		r = tools.retry_request(requests.get, ["http://provean.jcvi.org/provean_seq_report.php?jobid=" + jobid])
		with open("test.html", "w") as outf:
			outf.write(r.text)
		if '<div id="middle_content_template">' in r.text and 'PROVEAN Result' in r.text: #Success
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

#Inputs: list of continuous values x, list of binary values y, threshold that separates yes and no instances of x
#Outputs: balanced accuracy score of x values given threshold
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

#Input: dictionary of k: list of values
#Output: summary statistics of each list of values
def summarize_rand_data(data):
	summary = {k:{} for k in data.keys()}
	for k,v in data.items():
		try:
			summary[k] ={"Min":min(v), "5th %ile":np.percentile(v, 5), "Mean":statistics.mean(v), "Median":statistics.median(v), "95th %ile":np.percentile(v, 95),"Max":max(v), "StDev":statistics.stdev(v)}
		except ValueError as e:
			summary[k] = {"Min":"", "5th %ile":"", "Mean":"", "Median":"", "95th %ile":"", "Max":"", "StDev":""}
	return(summary)

#Input: long string, limiting number of characters per line
#Output: split up string that limits the number of characters per line
def split_long(s, lim=10):
	news = ""
	for si in s.split():
		if len(news) < 1 or len(news.split("\n")[-1]) + 1 + len(si) < lim:
			news += " " + si
		else:
			news += "\n" + si
	return(news)

#Input: dictionary of k: list of values
#Output: multiple box-and-whiskers plot for each key and list of values inputted
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

#Input: a value and a function
#Output: either the result of function on value, or some default if error
#useful when there is a good N/A value
def func_or_default(l, func=statistics.mean, default=0.5, verbose=False):
	try:
		return(func(l))
	except Exception as e:
		if verbose:
			print(e)
		return(default)

#Input: complex dictionary (usually from JSON), potential label of field of interest
#Output: clinical significance values from any fields labeled with the given label
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
