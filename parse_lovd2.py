import requests, re, os, sys, traceback, warnings, itertools, random, copy, bio_tools, stat_tools, statistics
import tools
from Bio import Entrez
import pandas as pd
import score_variants, time
import compute_features
from bs4 import BeautifulSoup

phenlabels = ["Diagnosis/Initial", "Diagnosis/Definite", "Inheritance", "Age/Examination", "Age/Diagnosis", "Phenotype details", "Age/Onset", "Phenotype/Onset"]
screenlabels = ["Template", "Technique", "Tissue", "Remarks", "Genes screened"]
#manual exclusion list
bad = [97732, 236883, 237268, 237282, 237283, 237789, 238701, 47105, 82622, 82675, 82876, 82877, 83042, 83660, 83795, 83801, 84150, 85187, 85538, 86837, 86842, 86843, 47307, 59243, 100302, 104076, 147153, 181005, 234464, 235372, 119100, 121892, 121894, 121932, 123231, 123248, 123308, 123309, 123314, 123321, 123322, 123494, 123497, 123505, 123506, 123510, 123512, 123514, 123516, 123521, 123526, 123536, 123542, 123544, 123547, 123548, 123550, 123551, 123552, 123560, 126239, 126240, 126242, 126243, 126246, 126248, 126249, 126250, 126251, 126253, 126978, 127010, 127798, 151782, 250289, 250292, 250293, 250294, 250400, 250566, 250567, 250666, 250687, 250749, 250754, 250761, 250767, 250772, 250808, 250809, 250810, 250848, 250869, 250876, 250902, 250918, 250945, 250949, 250954, 250975, 251013, 251027, 251050, 251057, 251066, 251074, 251084, 251106, 251121, 251148, 251149, 251180, 251181, 251189, 251199, 251206, 251227, 251230, 251235, 251257, 251266, 251273, 251300, 251325, 251339, 251344, 251346, 251352, 251360, 251362, 251380, 251383, 251397, 251441, 251444, 251464, 251479, 251505, 251522, 251544, 251581, 251584, 251590, 251598, 251601, 251624, 251631, 251660, 251682, 251692, 251696, 251722, 251737, 251739, 251755, 251829, 251835, 251837, 251838, 251841, 251875, 251879, 251884, 251887, 251890, 251918, 251935, 251947, 251949, 251968, 251970, 251971, 251984, 252029, 252086, 252092, 252112, 252166, 252217, 252349, 252355, 252357, 252364, 252378, 252387, 252dd417, 252418, 252430, 252443, 252460, 252467, 252479, 252484, 252508, 252551, 252592, 252616, 252638, 252643, 252657, 252687, 252696, 252719, 252751, 252756, 252765, 252768, 252801, 252828, 252851, 252864, 252870, 252890, 252916, 252938, 252962, 252970, 252974, 252977, 252988, 253011, 253045, 253051, 253106, 253110, 253135, 253137, 253138, 253142, 253181, 253231, 253267, 253327, 253355, 253359, 253381, 253408, 253417, 253423, 253426, 253445, 253446, 253499, 253527, 253530, 253534, 253549, 253599, 253661, 253664, 253810, 253834, 253835, 253873, 253876, 253907, 253974, 253979, 65317, 65346, 65349, 65362, 65363, 65365, 65372, 65393, 65399, 65428, 65438, 65440, 65465, 65480, 65554, 65565, 65588, 65656, 65671, 65692, 65745, 65774, 65794, 65819, 65862, 65868, 65876, 65903, 65917, 65924, 65935, 65964, 66016, 66065, 66114, 66131, 66176, 66196, 66210, 66275, 66279, 66296, 66310, 66311, 66323, 66324, 66340, 67074, 67155, 67325, 68201, 68567, 69057, 69741, 69750, 69805, 69835, 69857, 70605, 71003, 71548, 71583, 71614, 71680, 71700, 71717, 71740, 71782, 71809, 71813, 71814, 71859, 71860, 71892, 71898, 71946, 71954, 71967, 71974, 72023, 72033, 72048, 72055, 72059, 72067, 72127, 72128, 241659, 240096, 34389, 58828, 59231, 239002, 239028, 239034, 24255, 79024, 238993, 239040, 239041, 239042, 239091, 34050, 34392, 37583, 37587, 163818, 239048, 239064]

#association of genes with their transcript ID
runset = [("F8", "00495", "NM_000132.4"), ("F9", "02237", "NM_000133.4"), ("DMD", "05324", "NM_004006.3"), ("G6PD", "00308", "NM_000402.4"), ("HBB", "03506", "NM_000518.5"), ("VWF", "00498", "NM_000552.5"), ("VWF", "03361", "NM_000552.5"), ("VWF", "02126", "NM_000552.5")]
internal = [("F9", "./F8_F9/F9_dataset_updated.tsv", "NM_000133.4"), ("ADAMTS13", "./cTTP/Hing_w_neutral.tsv", "NM_139025.5"), ("VWF", "./VWD/VWF_dataset_updated.tsv", "NM_000552.5")]
manual_vars = {"03506":[["399A>T", "20A>T"], ["20A>T", "419A>G"], ["20A>T", "317T>C"], ["20A>T", "196A>G"], ["79G>A", "364G>C"]]}
for k,v in manual_vars.items():
	manual_vars[k] = {str(k) + "_" + str(i+1):{int(max(re.findall("\d+", mut), key=len)):re.findall("[ACGT]", mut) for mut in v[i]} for i in range(len(v))}

#manual inclusion list
good = {"00495":[82083, 82083, 82117, 82377, 82539, 82643, 82818, 84054, 84843, 84896], "02237":[235592, 235701, 235721, 235832, 235910, 235936, 235937, 235938, 235941, 235943, 235945, 236238, 236280, 236289, 236360, 236392, 236397, 236615, 236786, 236787, 237458, 237737, 237813, 237814], "05324":[122030, 122706, 123317, 123318, 123319, 123326, 126244, 126254, 144897, 174422], "00308":[4327, 4340, 4341, 250285, 250286, 250287, 250288, 250290, 250291, 250295, 250296, 250297], "00498":[46579, 47296, 47297, 47298, 47302, 59242, 59244, 59246, 79022, 79027, 79028, 180996, 181031, 181058, 181070, 234346, 234351, 234353, 234411, 234413, 234416, 234429, 234430, 234431, 234456, 234468, 235370], "03361":[151558, 240099, 240232, 240237, 240239, 240336, 240342, 244465, 34055, 47348, 74495, 79020, 100230, 100236, 239152, 239155, 239158, 239166, 239171, 239172, 240100, 240102, 240131, 240233, 240323, 240326, 240331, 240348, 240494, 241232, 241240, 241977, 244369, 244430, 244432, 244433, 244435, 244437, 244441, 244443, 244450, 244464, 244469, 244477, 244479, 244483, 244488, 244491, 244496, 244498], "02126":[239075, 239067, 24260, 239047]}

def parse_individual(individ, accid, ntseq=None, aaseq=None, index=1, multonly=True):
	accidversind = ".".join(accid.split(".")[:-1])
	data = {}
	translation = {}
	r2 = tools.retry_request(requests.get, ["https://databases.lovd.nl/shared/individuals/"+str(individ)])
	muts = re.findall(re.escape(accidversind) + "(\.\d+)?\:c\.(\d+[ACGT]\&gt\;[ACGT])", r2.text)
	muts = [mut[1].replace("&gt;", ">") for mut in muts]

	try:
		html_parsed = BeautifulSoup(str(r2.text))
		th = [x.split("</TH>")[0] for x in r2.text.split("<TH") if "DIV" in x]
		th = [th[i].replace("&nbsp;", " ") for i in range(len(th)) if "DIV" in th[i]]
		th = [re.findall("\>(.+)\<\/DIV\>", th[i]) for i in range(len(th))]
		th = [th[i] for i in range(len(th)) if len(th[i]) > 0]
		th = [th[i][j].strip() for i in range(len(th)) for j in range(len(th[i]))]
		phenorder = th[:th.index("Screening ID")-1]
		phenorder = {phenorder[j]:j for j in range(len(phenorder)) if phenorder[j] in phenlabels}
		screenorder = th[th.index("Screening ID")+1:th.index("Chr")]
		screenorder = {screenorder[j]:j for j in range(len(screenorder)) if screenorder[j] in screenlabels}
		varlabels = th[th.index("Chr"):]
		phenotypes = re.findall("\<td\>(.+)\<\/td", str(html_parsed.body.find_all("tr", attrs={"class":"data"})[0]))
		screenings = re.findall("\<td\>(.+)\<\/td", str(html_parsed.body.find_all("tr", attrs={"class":"data"})[1]))
		data = {k:phenotypes[phenorder[k]] for k in phenorder if not phenotypes[phenorder[k]].startswith("<span")}
		data.update({k:screenings[screenorder[k]] for k in screenorder if not screenings[screenorder[k]].startswith("<span")})
		
		for i in range(len(muts)):
			vardata = re.findall("\<td\>(.+)\<\/td", str(html_parsed.body.find_all("tr", attrs={"class":"data"})[i+2]))
			j = varlabels.index("Protein")
			if j >= 0 and j < len(vardata):
				translation[muts[i]] = vardata[j-2]
		
	except Exception as e:
		pass
		#tb = traceback.format_exc()
		#print("\033[93m" + str(tb) + "\033[0m")
		#print("\033[93m" + str(e) + "\033[0m")
	
	if ntseq != None and len(ntseq) > 0:
		for mut in muts:
			ntmut = tools.get_mutation_data(mut)
			if aaseq == None:
				aaseq = compute_features.translate_seq(ntseq)
			aamut = tools.get_mutant_aa(ntmut, ntseq, aaseq=aaseq, index=index)
			if aamut[1][0].lower() in ["stop", "*"] or aamut[1][1].lower() in ["stop", "*"]:#nonsense mutation
				return(None, None, None)
	muts = {int(max(re.findall("\d+", mut), key=len)):re.findall("[ACGT]", mut) for mut in muts}
	if multonly and len(muts) < 2:
		return(None, None, None)
	return(muts, data, translation)

def parse_disease(disid, accid, ntseq=None):
	page = 1
	pagegood = True
	ntvars = {}
	data = {}
	translation = {}
	accidversind = ".".join(accid.split(".")[:-1])
	while pagegood:
		try:
			r1 = tools.retry_request(requests.get, ["https://databases.lovd.nl/shared/ajax/viewlist.php?viewlistid=Individuals_for_D_VE&object=Individual&order=variants_%2CDESC&search_diseaseids_searched="+str(disid)+"&page_size=100&page="+str(page)])
			pagegood = False
			for line in r1.text.split("/TR"):
				cells = line.split("/TD")
				if len(cells) > 18:
					if int(re.findall("\>(\d+)\<", cells[18])[0]) > 1: #number of variants > 1
						pagegood = True
						individ = re.findall("id\=\"(\d+)\"", cells[0])[0] #individual id
						if int(individ) in bad:
							continue
						tmpntvars, tmpdata, tmptrans = parse_individual(individ, accidversind, ntseq=ntseq)
						if tmpntvars != None and tmpdata != None and tmptrans != None:
							ntvars[individ] = copy.deepcopy(tmpntvars)
							data[individ] = copy.deepcopy(tmpdata)
							translation.update(tmptrans)
						else:
							continue
			page += 1
		except Exception as e:
			print(str(page) + "\t" + str(e))
			page += 1
	ntvars = {k:v for k,v in ntvars.items() if len(v) > 1}
	
	for k,v in translation.items():
		posnt = re.findall("\d+", k)
		posaa = re.findall("\d+", v)
		if posnt != None and posaa != None and len(posnt) > 0 and len(posaa) > 0:
			posnt = int(max(posnt, key=len))
			posaa = int(max(posaa, key=len))
			if (posnt - 1)//3+1 != posaa:
				warnings.warn("Potential mismatch! " + str(k) + "\t" + str(v))

	return(ntvars, data, translation)

def parse_neutral(genename, accid, email="insert_here", conservative=True):
	Entrez.email=email
	if conservative:
		clinreq = "AND benign[CLIN] "
	else:
		clinreq = ""
	handle = tools.retry_func(Entrez.esearch, [], {'db':"snp", 'retmax':10000, 'term':genename+"[GENE] " + clinreq + "AND ( ( missense variant[Function_Class] OR synonymous variant[Function_Class] ) )"})
	record=Entrez.read(handle)
	ids = record["IdList"]
	prob = {}
	init_time = time.time()
	accidversind = ".".join(accid.split(".")[:-1])
	for i, snpid in enumerate(ids):
		try:
			freqd, clinsig, nm, rsid = bio_tools.query_dbSNP(snpid, refineorder=[])
			if any([True for i in range(len(clinsig)) if "pathogenic" in clinsig[i][1].lower()]) and not conservative:
				continue
			elif any([True for i in range(len(clinsig)) if "benign" not in clinsig[i][1].lower()]) and conservative:
				continue
			if ('Total',) not in freqd.keys() or freqd[('Total',)] == 0:
				freq = statistics.mean([v for k,v in freqd.items() if len(v) > 1 and v[1] == "Global"])
			else:
				freq = freqd[("Total",)]
			r = tools.retry_request(requests.get, ["https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/" + snpid])
			
			snpvar = re.findall(re.escape(accidversind) + "(\.\d+)?\:c\.(\d+[ACGT]>[ACGT])", str(r.text))
			snpvar = [snp[1] for snp in snpvar]
			if len(snpvar) > 1:
				if len(set([int(max(re.findall("\d+", snp),key=len)) for snp in snpvar])) > 1: #more than one position found
					print(snpvar)
				snpvar = snpvar[0]
			pos = int(max(re.findall("\d+", snpvar), key=len))
			nts = re.findall("[ACGT]", snpvar)
			prob[snpvar] = freq
		except (TypeError, statistics.StatisticsError,):
			pass
		except Exception as e:
			print(e)
		tools.update_time(i, len(ids), init_time)
	return(prob)


def generate_neutral(prob, disease, dups=True, ntseq="", aaseq=None, num=10000):
	count = {}
	vartype = {}
	if ntseq != None and len(ntseq) > 0:
		for ntvar in prob.keys():
			aavar = tools.get_mutant_aa(tools.get_mutation_data(ntvar), ntseq, aaseq=aaseq, index=1)
			if aavar[1][0] == aavar[1][1]:
				vartype[ntvar] = "Synonymous"
			else:
				vartype[ntvar] = "Missense"
	for patient, varset in disease.items():
		l = len(varset)
		if l not in count:
			count[l] = 1
		else:
			count[l] += 1
	t_n = stat_tools.range_tree(count)
	if dups:
		t_var = stat_tools.range_tree(prob)
		neutral = {}
		for i in range(len(disease)):
			num = stat_tools.sample_from_tree(t_n)
			neutral[i] = []
			while len(neutral[i]) < num:
				neutral[i].append(stat_tools.sample_from_tree(t_var))
				neutral[i] = list(set(neutral[i]))
		neutral = {"Neutral_"+str(i+1):{int(max(re.findall("\d+", mut), key=len)):re.findall("[ACGT]", mut) for mut in neutral[i]} for i in neutral}
		return(neutral)
	else:
		neutralset = list(itertools.chain.from_iterable(itertools.combinations(list(prob.keys()), r) for r in count.keys()))
		neutralset = {frozenset(k):np.prod([prob[ki] for ki in k]) for k in neutralset}
		neutralset = {(i,j):{k:v for k,v in neutralset.items() if len(k) == i+j and v > 0 and len([x for x in k if vartype[x] == "Synonymous"]) == i and len([x for x in k if vartype[x] == "Missense"]) == j} for j in range(max(count.keys())) for i in range(max(count.keys()))}
		t_var = {(i,j):stat_tools.range_tree(neutralset[(i,j)]) for (i,j) in neutralset}
		neutral = {}
		for i, (k,v) in enumerate(disease.items()):
			syn = 0
			mis = 0
			for pos, nts in v.items():
				aavar = get_mutant_aa((pos, nts), ntseq, aaseq=aaseq, index=1)
				if aavar[1][0] == aavar[1][1]:
					syn += 1
				else:
					mis += 1
			success = False
			if (syn, mis) in t_var.keys():
				try:
					neutral[i] = stat_tools.sample_from_tree(t_var[(syn, mis)])
					del neutralset[(syn, mis)][neutral[i]]
					t_var[(syn, mis)] = stat_tools.range_tree(neutralset[(syn, mis)])
					success = True
				except (TypeError, KeyError) as e:
					del t_var[(syn, mis)]
				if (syn, mis) in neutralset and len(neutralset[(syn, mis)]) == 0:
					del neutralset[(syn, mis)]
				if (syn, mis) in t_var and len(t_var[(syn, mis)]) == 0:
					del t_var[(syn, mis)]
			if not success:
				mlist = [x for x in (range(syn+mis+1))]
				random.shuffle(mlist)
				for m in mlist:
					n = syn+mis - m
					if (m,n) in t_var and (m,n) in neutralset:
						try:
							neutral[i] = stat_tools.sample_from_tree(t_var[(m, n)])
							del neutralset[(m, n)][neutral[i]]
							t_var[(m, n)] = stat_tools.range_tree(neutralset[(m, n)])
							success = True
						except (TypeError, KeyError) as e:
							del t_var[(m, n)]
						if (m, n) in neutralset and len(neutralset[(m, n)]) == 0:
							del neutralset[(m, n)]
						if (m, n) in t_var and len(t_var[(m, n)]) == 0:
							del t_var[(m, n)]
					if success:
						break
			
		neutral = {"Neutral_"+str(i+1):{tools.get_num(mut):re.findall("[ACGT]", mut) for mut in neutral[i]} for i in neutral}
		return(neutral)

def merge_disease_neutral(disease, neutral):
	variants = {}
	for patient, varset in (list(disease.items()) + list(neutral.items())):
		for pos, nts in varset.items():
			if (pos, tuple(nts)) not in variants:
				variants[(pos, tuple(nts))] = "c." + str(pos) + nts[0] + ">" + nts[1]
	data = {}
	for patient, varset in (list(disease.items()) + list(neutral.items())):
		if patient in disease:
			data[patient] = {"Disease":1}
		else:
			data[patient] = {"Disease":0}
		for variant, varname in variants.items():
			if int(variant[0]) in varset and list(variant[1]) == varset[variant[0]]:
				data[patient][varname] = 1
			else:
				data[patient][varname] = 0
	return(data)

def merge_neutrals(neutral1, neutral2, l=None):
	if l == None:
		l = len(neutral2)
	neutral1_set = {k:frozenset(["c." + str(k2) + v2[0] + ">" + v2[1] for k2, v2 in v.items()]) for k,v in neutral1.items()}
	neutral2_set = {k:frozenset(["c." + str(k2) + v2[0] + ">" + v2[1] for k2, v2 in v.items()]) for k,v in neutral2.items()}
	i = len(neutral1)
	for k,v in neutral2.items():
		if i >= l:
			break
		if neutral2_set[k] not in neutral1_set.values():
			neutral1["Neutral_" + str(i+1)] = neutral2[k]
		i = len(neutral1)
	return(neutral1)

def remove_dups(disease):
	test = list(disease[list(disease.keys())[0]].keys())[0]
	if isinstance(test, str):
		diseaseset = {k:frozenset([k2 for k2,v2 in v.items() if k2.startswith("c.") and tools.is_numeric(v2) and not np.isnan(v2) and v2 > 0]) for k,v in disease.items()}
	else:
		diseaseset = {k:frozenset(["c." + str(k2)  + v2[0] + ">" + v2[1] for k2, v2 in v.items()]) for k,v in disease.items()}
	diseasel = list(disease.keys())
	for i, k in enumerate(diseasel):
		for j in range(i):
			if diseaseset[diseasel[i]] == diseaseset[diseasel[j]]:
				try:
					del disease[diseasel[j]]
				except Exception as e:
					pass
	return(disease)

#parses combinations of variants from existing (local) table and generates corresponding neutral combinations
def pipeline_internal(internal=internal, email="insert_here"):
	data = {}
	for (gene, f, accid) in internal:
		ntseq = tools.read_fasta("./gene_data/" + gene + ".fasta")["ORF"].replace("-", "")
		disease = pd.read_csv(f, sep="\t", index_col=0)
		disease = disease.loc[disease["Disease"] == 1].T.to_dict()
		disease = remove_dups(disease)
		ntvars = {k:[k2 for k2, v2 in v.items() if k2.startswith("c.") and v2 > 0] for k,v in disease.items()}
		ntvars = {k:{tools.get_num(k2):re.findall("[ACGT]",k2) for k2 in v} for k,v in ntvars.items()}
		other = {k:{k2:v2 for k2, v2 in v.items() if not k2.startswith("c.")} for k,v in disease.items()}
		for k,v in list(ntvars.items()):
			for pos, nts in v.items():
				aamut = get_mutant_aa((pos, nts), ntseq, index=1)
				if aamut[1][0].lower() in ["stop", "*"] or aamut[1][1].lower() in ["stop", "*"]:#nonsense mutation
					del ntvars[k]
					del other[k]
		prob = parse_neutral(gene, accid, email=email, conservative=True)
		neutral = generate_neutral(prob, ntvars, dups=False, ntseq = ntseq)
		if len(neutral) < len(ntvars):
			prob = parse_neutral(gene, accid, email=email, conservative=False)
			tmpneutral = generate_neutral(prob, ntvars, dups=False, ntseq = ntseq)
			neutral = merge_neutrals(neutral, tmpneutral, l=len(ntvars))
		data[f] = merge_disease_neutral(ntvars, neutral)
		for k,v in data[f].items():
			if k in other:
				data[f][k].update(other[k])
		df = pd.DataFrame(data[f]).T
		df.to_csv(f, sep="\t")
	return(data)

#parses combinations of variants from LOVD, using manually curated combinations/patients where listed
def pipeline_manual(runset=runset, ntseq=None, email="insert_here"):
	for (gene, disid, accid) in runset:
		ntseq = tools.read_fasta("./gene_data/" + gene + ".fasta")["ORF"].replace("-", "")
		try:
			if disid in good:
				disease = {}
				data = {}
				translation = {}
				for individ in good[disid]:
					tmpntvars, tmpdata, tmptrans = parse_individual(individ, accid, ntseq=ntseq)
					if tmpntvars != None and tmpdata != None and tmptrans != None:
						disease[individ] = copy.deepcopy(tmpntvars)
						data[individ] = copy.deepcopy(tmpdata)
						translation.update(tmptrans)
					else:
						continue
			elif disid in manual_vars:
				disease = manual_vars[disid]
			disease = remove_dups(disease)
			prob = parse_neutral(gene, accid, email=email, conservative=True)
			neutral = generate_neutral(prob, disease, dups=False, ntseq = ntseq)
			if len(neutral) < len(disease):
				prob = parse_neutral(gene, accid, email=email, conservative=False)
				tmpneutral = generate_neutral(prob, disease, dups=False, ntseq = ntseq)
				neutral = merge_neutrals(neutral, tmpneutral, l=len(disease))
			#if len(prob) < max([len(v) for v in disease.values()]):
			#	raise ValueError("Not enough neutral variants found to generate control set")
			data = merge_disease_neutral(disease, neutral)
			df = pd.DataFrame(data).T
			#df.reindex_axis(sorted(df.columns, key=lambda x: int(max(re.findall("\d+", str(x)), key=len))), axis=1)
			if not os.path.isdir("./LOVD/"):
				os.mkdir("./LOVD/")
			df.to_csv("./LOVD/" + str(disid) + "_" + gene + ".tsv", sep="\t")
		except Exception as e:
			exc_type, exc_value, exc_traceback = sys.exc_info()
			traceback.print_tb(exc_traceback)
			print(str(gene) + "\t" + str(disid) + "\t" + str(e))

#parses combinations of variants from LOVD and generates corresponding set of neutral combinations
#needs disease ID from LOVD and Refseq accession ID	
def pipeline(genename, disid, accid, email="insert_here"):
	disease, data, translation = parse_disease(disid, accid)
	disease = remove_dups(disease)
	prob = parse_neutral(genename, accid, email=email, conservative=True)
	neutral = generate_neutral(prob, disease)
	if len(neutral) < len(disease):
		prob = parse_neutral(gene, accid, email=email, conservative=False)
		tmpneutral = generate_neutral(prob, disease, dups=False, ntseq = ntseq)
		neutral = merge_neutrals(neutral, tmpneutral, l=len(disease))
	#if len(prob) < max([len(v) for v in disease.values()]):
	#	raise ValueError("Not enough neutral variants found to generate control set")
	data = merge_disease_neutral(disease, neutral)
	df = pd.DataFrame(data).T
	snpcols = [x for x in df.columns if x.startswith("c.")]
	othercols = [x for x in df.columns if not x.startswith("c.")]
	#df.reindex_axis(sorted(df.columns, key=lambda x: int(max(re.findall("\d+", str(x)), key=len))), axis=1)
	if not os.path.isdir("./LOVD/"):
		os.mkdir("./LOVD/")
	df.to_csv("./LOVD/" + str(disid) + "_" + genename + ".tsv", sep="\t")
	print("Finished with " + genename + ": " + disid)

