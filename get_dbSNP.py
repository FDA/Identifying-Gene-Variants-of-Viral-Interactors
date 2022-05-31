import requests, os, sys, re, time, math, tools, bio_tools, compute_features, warnings, psutil, traceback
from multiprocessing import Process
from Bio import Entrez
curpath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, curpath + '/EVmutation')
import score_variants, compute_features, create_analysis

animals = ["Papio anubis", "canis lupus familiaris", "rhinolophus sinicus", "Mesocricetus auratus",  "mus musculus", "Rattus norvegicus", "Callithrix jacchus", "Oryctolagus cuniculus", "Cavia porcellus", "Mustela putorius furo", "Felis catus", "Manis javanica", "Macaca mulatta"]

def is_cds(seq):
	if len(seq) % 3 == 0 and compute_features.translate_seq(seq)[0] == "M" and compute_features.translate_seq(seq)[-1] == "*":
		return True
	else:
		return False

def get_animal_seq(genename, organism, taxid=None, wait_time=1):
	searchhandle = tools.retry_func(Entrez.esearch, [], {'db':"nuccore", 'term':genename + "[Gene Name] AND " + organism + "[Organism]"}, wait=5)
	if searchhandle != None:
		searchhandle = searchhandle.read()
	else:
		return(None)
	ids = re.findall("\<Id\>(\d+)\<\/Id\>", str(searchhandle))
	#print(str(organism) + ":\t" + str(ids))
	if len(ids) > 0:
		for i in range(len(ids)):
			time.sleep(wait_time)
			try:
				handle = tools.retry_func(Entrez.efetch, keyword_arguments={'db':"nuccore", 'id':ids[i], 'retmode':"xml"}, wait=5).read().decode('utf-8')
				seq = max(re.findall("\<GBSeq\_sequence\>([a-zA-Z]+)\<\/GBSeq\_sequence\>", handle), key=len).upper()
				if is_cds(seq):
					return(seq)
				else:
					locations = re.findall("\<GBFeature\_location\>(\d+\.+\d+)\<\/GBFeature\_location\>", handle)
					#print(ids[i] + ":\t" + str(locations))
					for location in locations:
						location = [int(x) for x in re.findall("\d+", location)]
						if is_cds(seq[location[0]-1:location[1]]):
							return(seq[location[0]-1:location[1]])
			except Exception as e:
				print(str(ids[i]) + ":\t" + str(e))
	return(None)

def get_frequencies(variants):
	server = "https://rest.ensembl.org"
	ext = "/vep/human/hgvs"
	headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
	data = {}
	init_time = time.time()
	for i, var in enumerate(variants):
		try:
			r = tools.retry_request(requests.post, positional_arguments=[server+ext], keyword_arguments={"headers":headers, "data":'{ "hgvs_notations" : ' + str([var]).replace("\'", "\"") + ' }'})
			decoded = r.json()
			for j in range(len(decoded)):
				for k in range(len(decoded[j]["colocated_variants"])):
					try:
						data[var] = decoded[j]["colocated_variants"][k]["frequencies"]
						break
					except KeyError as e:
						pass
					except Exception as e:
						print(e)
		except Exception as e:
			print(e)
		update_time(i, len(variants), init_time)
	return(data)

def get_GOs(pdbs):
	if isinstance(pdbs, dict) and not isinstance(pdbs, (list, tuple,)):
		pdbs = [[k.lower(), v["chain"].upper()] for k,v in pdbs.items()]
	data = {}
	d = {}
	init_time = time.time()
	for i, pdb in enumerate(pdbs):
		accids = []
		r = tools.retry_request(requests.get, ["https://www.rcsb.org/pdb/rest/describeMol?structureId=" + pdb[0].lower()])
		pieces = r.text.split("<polymer")
		for piece in pieces:
			try:
				chains = re.findall("\<chain id\=\"(.+)\"", piece)
				if len(chains) > 0 and pdb[1].upper() in [x.upper() for x in chains]:
					accid = max(re.findall("\<accession id\=\"([A-Z]\d+)\"", r.text), key=len)
					gos = tools.retry_request(requests.get, ["http://api.geneontology.org/api/bioentity/gene/" + accid.upper() + "/function"])
					data[accid] = re.findall("GO\:\s*(\d+)", str(gos.text))
					label = re.findall("\"id\"\: \"GO\:\d+\"\, \"label\"\: \".+?\"\,", str(gos.text))
					d.update({max(re.findall("GO\:(\d+)", k), key=len):max(re.findall("label\"\: \"(.+?)\"", k), key=len) for k in label})
			except Exception as e:
					print(str(pdb) + "\t" + str(e))
		update_time(i, len(pdbs), init_time)
	return(data, d)

def get_genename_from_pdb(pdblist, taxid="9606"):
	letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	names = {}
	for i in range(math.ceil(len(pdblist)/20)):
		pdbsubset = pdblist[20*i:min(len(pdblist), 20*(i+1))]
		r = requests.get("https://data.rcsb.org/graphql?query=%7B%0A%20%20polymer_entities(entity_ids%3A%5B%22" + "%22%2C%22".join([pdb[0].upper() + "_" + str(letters.index(pdb[1])+1) for pdb in pdbsubset])+ "%22%5D)%20%7B%0A%20%20%20%20rcsb_id%0A%20%20%20%20entity_src_gen%20%7B%0A%20%20%20%20%20%20pdbx_host_org_gene%0A%20%20%20%20%7D%0A%20%20%20%20rcsb_entity_source_organism%20%7B%0A%20%20%20%20%20%20ncbi_taxonomy_id%0A%20%20%20%20%20%20ncbi_scientific_name%0A%20%20%20%20%20%20rcsb_gene_name%20%7B%0A%20%20%20%20%20%20%20%20provenance_source%0A%20%20%20%20%20%20%20%20value%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%7D%0A%20%20%20%20rcsb_cluster_membership%20%7B%0A%20%20%20%20%20%20cluster_id%0A%20%20%20%20%20%20identity%0A%20%20%20%20%7D%0A%20%20%7D%0A%7D")
		d = r.json()["data"]["polymer_entities"]
		for i in range(len(d)):
			name = d[i]['rcsb_id'].split("_")[0].lower()
			name = name + "_" + [xi[1] for xi in pdblist if name.lower() == xi[0].lower()][0]
			names[name] = []
			d2 = d[i]["rcsb_entity_source_organism"]
			for j in range(len(d2)):
				if str(d2[j]["ncbi_taxonomy_id"]) == str(taxid):
					d3 = d2[j]["rcsb_gene_name"]
					for k in range(len(d3)):
						names[name].append(d3[k]["value"])
	return({k:list(set(v)) for k,v in names.items()})
				
		

def get_syn(genename, seq=None, refid=None, path="/media/temp/", outdir="/media/home/workspace/Coronavirus/interactome/", space = "-", index=1, quiet=False):
	if seq == None or refid == None:
		seqids = bio_tools.get_accids(genename, organism="homo sapiens")
		refids = sorted(seqids["mRNA"], key=lambda kv: int(max(re.findall("\d+", kv[0]), key=len)))
		refids = [x[0] for x in refids]
		seqs, cdses, mRNA, protein = bio_tools.pipeline(genename)
		seq = seqs["ORF"]
	else:
		refids = [refid]

	blastFlag = False
	if not os.path.exists(outdir + genename + "_nt_msa.fasta") or tools.read_fasta(outdir + genename + "_nt_msa.fasta")[genename].replace("-", "") != seq:
		p_nt = Process(target=compute_features.run_blast, args=(seq, genename, "nt", "nt", "blastn", True, None, path, space))
		p_nt.start()
		blastFlag = True

	r = tools.retry_request(requests.get, ['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=(((' + genename.lower() + '%5BGene%20Name%5D)%20AND%20"snv"%5BSNP%20Class%5D)%20AND%20"synonymous%20variant"%5BFunction%20Class%5D)&retmax=' + str(len(seq)*4)], {})
	ids = re.findall("\<Id\>(\d+)\<\/Id\>", r.text)
	
	init_time = time.time()
	variants = []
	for i in range(len(ids)//50):
		idsub = ids[50*i:min(50*(i+1), len(ids)+1)]
		r = tools.retry_request(requests.post, ["https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=" + ",".join(idsub)], {})
		possible_variants = re.findall("NM\_\d+\.?\d*\:c\.\d+[A-Z]\&gt\;[A-Z]", r.text) + re.findall("NM\_\d+\.?\d*\:c\.\d+[A-Z]\>[A-Z]", r.text)
		for i, refid in enumerate(refids):
			variants += [x.replace("&gt;", ">") for x in possible_variants if refid.split(".")[0] == x.split(".")[0] or refid.split(":")[0] == x.split(":")[0]]
			variants = [re.sub("NM\_.+\:c\.", refid + ":c.", x) for x in variants]
			variants = sorted(set(variants), key=lambda k: int(max(re.findall("\:c\.(\d+)", k), key=len)))
			if len(variants) > 0:
				break
		if len(variants) == 0:
			print(genename)
			print(possible_variants)
		tools.update_time(i, len(ids)//50, init_time)

	server = "https://rest.ensembl.org"
	ext = "/vep/human/hgvs"
	headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
	data = {}
	init_time = time.time()
	varquery = [genename + ":c." + var.split(":c.")[-1] for var in variants]
	try:
		r = tools.retry_request(requests.post, positional_arguments=[server+ext], keyword_arguments={"headers":headers, "data":'{ "hgvs_notations" : ' + str(varquery).replace("\'", "\"") + ' }'})
		decoded = r.json()
		for i in range(len(decoded)):
			try:
				var = tools.seek_in_struct(decoded[i], ['input'])["input"]
				data[var] = tools.seek_in_struct(decoded[i], ["frequencies"])["frequencies"]
				data[var].update(tools.seek_in_struct(decoded[i], ['motif_score_change', 'polyphen_score', 'sift_score']))
				break
			except KeyError as e:
				pass
			except Exception as e:
				print(e)
	except KeyError as e:
		pass
	except Exception as e:
		print(e)
	data = {k1.split(">")[0] + k2:v2 for k1 in data for k2, v2 in data[k1].items()}
	#pd.DataFrame(data).T.to_csv(genename + "_syn.tsv", sep="\t")
	conservation = {}
	if blastFlag:
		p_nt.join()
		nt_seqs = tools.read_fasta(path + "nt_msa.fasta", aformat="WHOLE")
	else:
		nt_seqs = tools.read_fasta(outdir + genename + "_nt_msa.fasta")
	for organism in animals:
		if organism not in nt_seqs.keys():
			animalseq = get_animal_seq(genename, organism, taxid=None, wait_time=1)
			if animalseq != None:
				nt_seqs[organism] = animalseq
	if any([True if len(nt_seqs[k]) != len(nt_seqs[genename]) else False for k in nt_seqs.keys()]):
		nt_seqs = {k:v.replace("-", "") for k,v in nt_seqs.items()}
	ntalignment, conservation["NT conservation"], conservation["NT entropy"], conservation["NT variance"] = compute_features.compute_conservation(nt_seqs, genename, mode="nt")
	dist = score_variants.compute_distribution(ntalignment, genename, space="-", alphabet="")
	conservation["Rare codon enrichment"] = compute_features.rc_enrichment(nt_seqs, genename, quiet=True)

	if blastFlag:
		write_fasta(ntalignment, outdir + genename + "_nt_msa.fasta")
	'''
	if not os.path.exists(outdir + genename + "_outfile.csv"):
		try:
			available_memory = (psutil.virtual_memory().available/1073741824)
			if available_memory <= (1/806.0)**2*(len(seq)**2):
				raise RuntimeError("Empirical memory limit reached for PLMC.")
			create_analysis.create_analysis(fasta=outdir + genename + "_nt_msa.fasta", focus=genename, skip_align_muts=True, model_params=outdir + genename + "_nt.model_params", mode="matrix", alphabet="-TAGC")
		except RuntimeError as e:
			print("\033[93m" + str(e) + "\033[0m")
	evscores = {}
	
	try:
		ev_data = pd.read_csv(outdir + genename + "_outfile.csv", sep=",")
		for i, col in enumerate(ev_data):
			if col.strip() in ["A", "G", "C", "T"] and i > 0:
				evscores[col] = list(ev_data[col])
	except FileNotFoundError as e:
		print("\033[93m" + str(e) + "\033[0m")
	'''
	codondata = pd.read_csv("/media/home/workspace/DB/sources/codondata.csv", sep="\t", header=0, index_col=0) #file containing the genomic codon data
	cpdata = pd.read_csv("/media/home/workspace/DB/sources/codonpairdata.csv", sep="\t", header=0, index_col=0) #file containing the genomic codon pair data
	minmax = compute_features.get_minmax(seq, genename)

	init_time = time.time()
	for j, var in enumerate(data.keys()):
		mut = get_mutation_data(str(var.split(":c.")[-1]))
		if len(seq) <= mut[0] - index:
			warnings.warn("Variant " + str(var) + " beyond end of sequence " + str(genename))
			continue
		elif seq[mut[0] - index] != mut[1][0]:
			warnings.warn("NTs don't match for " + str(var) + ", actual value is " + seq[mut[0] - index])
			
		for k in ["NT conservation", "NT entropy", "NT variance"]:
			data[var][k] = conservation[k][mut[0]-index]
		data[var]["Rare codon enrichment"] = conservation["Rare codon enrichment"][(mut[0]-index)//3]
		try:
			data[var]["MSA likelihood"] = dist[mut[0]][mut[1][1]]
		except KeyError as e:
			data[var]["MSA likelihood"] = 0
		try:
			data[var]["%MinMax"] = minmax["%minmax"][(mut[0]-index)//3]
			data[var]["%MinMax control"] = minmax["%minmax control"][(mut[0]-index)//3]
		except KeyError as e:
			data[var]["%MinMax"] = data[var]["%MinMax control"] = ""
		substring = subseq(seq, mut[0]-index, 75)
		mutstring = subseq(update_str(seq, mut[1][1], mut[0]-index), mut[0]-index, 75)
		try:
			data[var]["delta mRNA MFE (Kinefold)"] = compute_features.run_kinefold(mutstring) - compute_features.run_kinefold(substring)
		except RuntimeError as e:
			print(str(var) + "\t" + str(e))
		data[var]["delta mRNA MFE (RNAfold)"] = compute_features.get_RNAfold(mutstring) - compute_features.get_RNAfold(substring)
		data[var]["delta mRNA MFE (NUPACK)"] = compute_features.run_nupack(mutstring) - compute_features.run_nupack(substring)

		mutseq = update_str(seq, mut[1][1], mut[0]-index)

		#compute codon and codon pair for WT and mutant
		posincodon = (mut[0] -index) % 3
		codonstart = (mut[0] - index) - posincodon
		c_WT = seq[codonstart:codonstart+3]
		cp1_WT = seq[codonstart-3:codonstart+3]
		cp2_WT = seq[codonstart:codonstart+6]

		c_mut = mutseq[codonstart:codonstart+3]
		cp1_mut = mutseq[codonstart-3:codonstart+3]
		cp2_mut = mutseq[codonstart:codonstart+6]

		for column in codondata:
			try:
				data[var]["Δ " + column] = codondata.at[c_mut, column] - codondata.at[c_WT, column]
			except:
				data[var]["Δ " + column] = ""
		for column in cpdata:
			try: 
				data[var]["Δ " + column + " 1"] = cpdata.at[cp1_mut, column] - cpdata.at[cp1_WT, column]
			except:
				data[var]["Δ " + column + " 1"] = ""
			try:
				data[var]["Δ " + column + " 2"] = cpdata.at[cp2_mut, column] - cpdata.at[cp2_WT, column]
			except:
				data[var]["Δ " + column + " 2"] = ""
		tools.update_time(j, len(data), init_time)
	'''
		if len(evscores) > 0:
			data[var]["EVmutation(nt)"] = evscores[mut[1][1]][mut[0] - index]
	'''
	pd.DataFrame(data).T.to_csv(outdir + genename + "_syn.tsv", sep="\t")
	return(data)

def compute_codon(seq, data, genename, index=1):
	codondata = pd.read_csv("/media/home/workspace/DB/sources/codondata.csv", sep="\t", header=0, index_col=0) #file containing the genomic codon data
	cpdata = pd.read_csv("/media/home/workspace/DB/sources/codonpairdata.csv", sep="\t", header=0, index_col=0) #file containing the genomic codon pair data
	minmax = compute_features.get_minmax(seq, genename)
	for var in data.keys():
		mut = get_mutation_data(str(var.split(":c.")[-1]))
		if len(seq) <= mut[0] - index:
			warnings.warn("Variant " + str(var) + " beyond end of sequence " + str(genename))
			continue
		elif seq[mut[0] - index] != mut[1][0]:
			warnings.warn("NTs don't match for " + str(var) + ", actual value is " + seq[mut[0] - index])
		try:
			data[var]["%MinMax"] = minmax["%minmax"][(mut[0]-index)//3]
			data[var]["%MinMax control"] = minmax["%minmax control"][(mut[0]-index)//3]
		except KeyError as e:
			data[var]["%MinMax"] = data[var]["%MinMax control"] = ""
		substring = subseq(seq, mut[0]-index, 75)
		mutstring = subseq(update_str(seq, mut[1][1], mut[0]-index), mut[0]-index, 75)
		try:
			data[var]["delta mRNA MFE (Kinefold)"] = compute_features.run_kinefold(mutstring) - compute_features.run_kinefold(substring)
		except RuntimeError as e:
			print(str(var) + "\t" + str(e))
		data[var]["delta mRNA MFE (RNAfold)"] = compute_features.get_RNAfold(mutstring) - compute_features.get_RNAfold(substring)
		data[var]["delta mRNA MFE (NUPACK)"] = compute_features.run_nupack(mutstring) - compute_features.run_nupack(substring)

		mutseq = update_str(seq, mut[1][1], mut[0]-index)

		#compute codon and codon pair for WT and mutant
		posincodon = (mut[0] -index) % 3
		codonstart = (mut[0] - index) - posincodon
		c_WT = seq[codonstart:codonstart+3]
		cp1_WT = seq[codonstart-3:codonstart+3]
		cp2_WT = seq[codonstart:codonstart+6]

		c_mut = mutseq[codonstart:codonstart+3]
		cp1_mut = mutseq[codonstart-3:codonstart+3]
		cp2_mut = mutseq[codonstart:codonstart+6]

		for column in codondata:
			try:
				data[var]["Δ " + column] = codondata.at[c_mut, column] - codondata.at[c_WT, column]
			except:
				data[var]["Δ " + column] = ""
		for column in cpdata:
			try: 
				data[var]["Δ " + column + " 1"] = cpdata.at[cp1_mut, column] - cpdata.at[cp1_WT, column]
			except:
				data[var]["Δ " + column + " 1"] = ""
			try:
				data[var]["Δ " + column + " 2"] = cpdata.at[cp2_mut, column] - cpdata.at[cp2_WT, column]
			except:
				data[var]["Δ " + column + " 2"] = ""
	return(data)

def get_nonsyn(genename, seq=None, refid=None, path="/media/temp/", outdir="/media/home/workspace/Coronavirus/interactome/", space = "-", index=1, quiet=False):
	if seq == None or refid == None:
		seqids = bio_tools.get_accids(genename, organism="homo sapiens")
		refids = sorted(seqids["mRNA"], key=lambda kv: int(max(re.findall("\d+", kv[0]), key=len)))
		refids = [x[0] for x in refids]
		seqs, cdses, mRNA, protein = bio_tools.pipeline(genename)
		seq = seqs["ORF"]
	else:
		refids = [refid]
	aaseq = compute_features.translate_seq(seq)

	blastFlag = False
	if not os.path.exists(outdir + genename + "_aa_msa.fasta") or tools.read_fasta(outdir + genename + "_aa_msa.fasta")[genename].replace("-", "") != aaseq:
		p_aa = Process(target=compute_features.run_blast, args=(seq, genename, "aa", "nr", "blastp", True, None, path, space))
		p_aa.start()
		blastFlag = True
	else:
		aa_seqs = tools.read_fasta(outdir + genename + "_aa_msa.fasta")

	r = tools.retry_request(requests.get, ['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=snp&term=(((' + genename.lower() + '%5BGene%20Name%5D)%20AND%20"snv"%5BSNP%20Class%5D)%20AND%20"missense%20variant"%5BFunction%20Class%5D)&retmax=' + str(len(seq)*4)], {})
	ids = re.findall("\<Id\>(\d+)\<\/Id\>", r.text)

	init_time = time.time()
	variants = []
	for i in range(len(ids)//50):
		idsub = ids[50*i:min(50*(i+1), len(ids)+1)]
		r = tools.retry_request(requests.post, ["https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=" + ",".join(idsub)], {})
		possible_variants = re.findall("NM\_\d+\.?\d*\:c\.\d+[A-Z]\&gt\;[A-Z]", r.text) + re.findall("NM\_\d+\.?\d*\:c\.\d+[A-Z]\>[A-Z]", r.text)
		for i, refid in enumerate(refids):
			variants += [x.replace("&gt;", ">") for x in possible_variants if refid.split(".")[0] == x.split(".")[0] or refid.split(":")[0] == x.split(":")[0]]
			variants = [re.sub("NM\_.+\:c\.", refid + ":c.", x) for x in variants]
			variants = sorted(set(variants), key=lambda k: int(max(re.findall("\:c\.(\d+)", k), key=len)))
			if len(variants) > 0:
				break
		if len(variants) == 0:
			print(genename)
			print(possible_variants)
		tools.update_time(i, len(ids)//50, init_time)

	server = "https://rest.ensembl.org"
	ext = "/vep/human/hgvs"
	headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
	data = {}
	init_time = time.time()
	varquery = [genename + ":c." + var.split(":c.")[-1] for var in variants]
	try:
		r = tools.retry_request(requests.post, positional_arguments=[server+ext], keyword_arguments={"headers":headers, "data":'{ "hgvs_notations" : ' + str(varquery).replace("\'", "\"") + ' }'})
		decoded = r.json()
		for i in range(len(decoded)):
			try:
				var = tools.seek_in_struct(decoded[i], ['input'])["input"]
				data[var] = tools.seek_in_struct(decoded[i], ["frequencies"])["frequencies"]
				data[var].update(tools.seek_in_struct(decoded[i], ['motif_score_change', 'polyphen_score', 'sift_score']))
				break
			except KeyError as e:
				pass
			except Exception as e:
				print(e)
	except KeyError as e:
		pass
	except Exception as e:
		print(e)

	data = {k1.split(">")[0] + k2:v2 for k1 in data for k2, v2 in data[k1].items()}
	#pd.DataFrame(data).T.to_csv(genename + "_nonsyn.tsv", sep="\t")
	codondata = pd.read_csv("/media/home/workspace/DB/sources/codondata.csv", sep="\t", header=0, index_col=0) #file containing the genomic codon data
	cpdata = pd.read_csv("/media/home/workspace/DB/sources/codonpairdata.csv", sep="\t", header=0, index_col=0) #file containing the genomic codon pair data
	minmax = compute_features.get_minmax(seq, genename)
	conservation = {}

	try:
		conservation["O-linked Glycosylation potential"] = compute_features.run_netOglyc(seq)
		conservation["Phosphorylation potential"] = compute_features.run_netphos(seq, path=path)
		conservation["N-linked Glycosylation potential"] = compute_features.run_netNglyc(seq, path=path)
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")
	
	if not quiet:
		print("Computing surface area with NetSurfP2")
	try:
		conservation["Relative surface area (NetsurfP2)"], conservation["Relative surface area (NetsurfP2)"], conservation["SS3 (NetsurfP2)"], conservation["SS8 (NetsurfP2)"], conservation["Disorder"] = compute_features.run_netsurfp(seq, genename)
	except Exception as e:
		tb = traceback.format_exc()
		print("\033[93m" + str(tb) + "\033[0m")
		print("\033[93m" + str(e) + "\033[0m")

	if blastFlag:
		p_aa.join()
		aa_seqs = tools.read_fasta(path + "aa_msa.fasta", aformat="WHOLE")
	for organism in animals:
		if organism not in aa_seqs.keys():
			animalseq = get_animal_seq(genename, organism, taxid=None, wait_time=1)
			if animalseq != None:
				aa_seqs[organism] = compute_features.translate_seq(animalseq.replace("-", ""))
	if any([True if len(aa_seqs[k]) != len(aa_seqs[genename]) else False for k in aa_seqs.keys()]):
		aa_seqs = {k:v.replace("-", "") for k,v in aa_seqs.items()}
	aaalignment, conservation["AA conservation"], conservation["AA entropy"], conservation["AA variance"], conservation["BLOSUM"] = compute_features.compute_conservation(aa_seqs, genename, mode="aa")
	dist = score_variants.compute_distribution(aaalignment, genename, space="-", alphabet="")

	if blastFlag:
		write_fasta(aaalignment, outdir + genename + "_aa_msa.fasta")
	'''
	if not os.path.exists(outdir + genename + "_outfile.csv"):
		try:
			available_memory = (psutil.virtual_memory().available/1073741824)
			if available_memory <= (1/806.0)**2*(len(seq)**2):
				raise RuntimeError("Empirical memory limit reached for PLMC.")
			create_analysis.create_analysis(fasta=outdir + genename + "_aa_msa.fasta", focus=genename, skip_align_muts=True, model_params=outdir + genename + "_aa.model_params", mode="matrix", alphabet="-TAGC")
		except RuntimeError as e:
			print("\033[93m" + str(e) + "\033[0m")
	evscores = {}
	
	try:
		ev_data = pd.read_csv(outdir + genename + "_outfile.csv", sep=",")
		for i, col in enumerate(ev_data):
			if re.match("[A-Z]", str(col.strip())) and i > 0:
				evscores[col] = list(ev_data[col])
	except FileNotFoundError as e:
		print("\033[93m" + str(e) + "\033[0m")
	'''
	
	init_time = time.time()
	for j, var in enumerate(data.keys()):
		mut = get_mutation_data(str(var.split(":c.")[-1]))
		aamut = get_mutant_aa(mut, seq, aaseq=aaseq, index=1)
		if len(seq) <= mut[0] - index:
			warnings.warn("Variant " + str(var) + " beyond end of sequence " + str(genename))
			continue
		elif seq[mut[0] - index] != mut[1][0]:
			warnings.warn("NTs don't match for " + str(var) + ", actual value is " + seq[mut[0] - index])
			
		for k in ["AA conservation", "AA entropy", "AA variance"]:
			data[var][k] = conservation[k][aamut[0]-index]
		try:
			data[var]["MSA likelihood"] = dist[aamut[0]][aamut[1][1]]
		except KeyError as e:
			data[var]["MSA likelihood"] = 0
		try:
			data[var]["%MinMax"] = minmax["%minmax"][(mut[0]-index)//3]
			data[var]["%MinMax control"] = minmax["%minmax control"][(mut[0]-index)//3]
		except KeyError as e:
			data[var]["%MinMax"] = data[var]["%MinMax control"] = ""
		substring = subseq(seq, mut[0]-index, 75)
		mutstring = subseq(update_str(seq, mut[1][1], mut[0]-index), mut[0]-index, 75)
		try:
			data[var]["delta mRNA MFE (Kinefold)"] = compute_features.run_kinefold(mutstring) - compute_features.run_kinefold(substring)
		except RuntimeError as e:
			print(str(var) + "\t" + str(e))
		data[var]["delta mRNA MFE (RNAfold)"] = compute_features.get_RNAfold(mutstring) - compute_features.get_RNAfold(substring)
		data[var]["delta mRNA MFE (NUPACK)"] = compute_features.run_nupack(mutstring) - compute_features.run_nupack(substring)

		mutseq = update_str(seq, mut[1][1], mut[0]-index)

		#compute codon and codon pair for WT and mutant
		posincodon = (mut[0] -index) % 3
		codonstart = (mut[0] - index) - posincodon
		c_WT = seq[codonstart:codonstart+3]
		cp1_WT = seq[codonstart-3:codonstart+3]
		cp2_WT = seq[codonstart:codonstart+6]

		c_mut = mutseq[codonstart:codonstart+3]
		cp1_mut = mutseq[codonstart-3:codonstart+3]
		cp2_mut = mutseq[codonstart:codonstart+6]

		for column in codondata:
			try:
				data[var]["Δ " + column] = codondata.at[c_mut, column] - codondata.at[c_WT, column]
			except:
				data[var]["Δ " + column] = ""
		for column in cpdata:
			try: 
				data[var]["Δ " + column + " 1"] = cpdata.at[cp1_mut, column] - cpdata.at[cp1_WT, column]
			except:
				data[var]["Δ " + column + " 1"] = ""
			try:
				data[var]["Δ " + column + " 2"] = cpdata.at[cp2_mut, column] - cpdata.at[cp2_WT, column]
			except:
				data[var]["Δ " + column + " 2"] = ""
		vep = compute_features.get_vep(mut[0], mut[1], genename, epsilon=0.001)
		if vep != None:
			data[var].update(vep)
		tools.update_time(j, len(data), init_time)
	'''
		if len(evscores) > 0:
			data[var]["EVmutation(aa)"] = evscores[mut[1][1]][mut[0] - index]
	'''
	pd.DataFrame(data).T.to_csv(outdir + genename + "_nonsyn.tsv", sep="\t")
	return(data)

def gen_seqs(genomic, locations, space="-"):
	orf = "".join([genomic[locations[i][0]:locations[i][1]] for i in range(len(locations))])
	transcript = orf
	aligned_seqs = tools.align_sequences({"genomic_aligned":genomic[locations[0][0]:locations[-1][1]], "ORF_aligned":orf, "transcript_aligned": transcript})
	aligned_seqs["ORF"] = aligned_seqs["ORF_aligned"].replace(space, "")
	aligned_seqs["transcript"] = aligned_seqs["transcript_aligned"].replace(space, "")
	aligned_seqs["genomic"] = aligned_seqs["genomic_aligned"].replace(space, "")
	return(aligned_seqs)
