import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse, openpyxl, re, os, os.path, sys, time, datetime, scipy.stats, statistics, random, statsmodels, copy, openpyxl, tools
import statsmodels.api as sm
import sklearn.metrics
from sklearn.impute import SimpleImputer
from sklearn import preprocessing
from tools import slugify, is_numeric, pickleload, padstring
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['font.family'] = "arial"

#used to sort patients
nt_to_frac={"A":0.2, "C":0.4, "G":0.6, "T":0.8}

def merge_files(directory, regex, sep="\t"):
	now = datetime.datetime.now()	
	newexcel = os.path.join(directory, tools.slugify(regex)) + ".xlsx"
	writer = pd.ExcelWriter(newexcel, engine="openpyxl")
	for f in os.listdir(directory):
		f = os.path.basename(f)
		if re.search(regex, f) and not "f".endswith(".xlsx"):
			try:
				df = pd.read_csv(os.path.join(directory, f), sep=sep, index_col=0)
				model = max(re.findall(regex, f), key=len)	
				df.to_excel(writer, str(model))
			except Exception as e:
				print(e)
	print(newexcel)
	writer.save()
	writer.close()

def unique_haplo(df, binary=None, multonly=True):
	df = pd.read_csv(df, sep="\t", index_col=0).T.to_dict()
	if multonly:
		df = {k:v for k,v in df.items() if sum([v[x] for x in v.keys() if x.startswith("c.") and tools.is_numeric(v[x]) and not np.isnan(v[x])]) >= 2}
	for k,v in df.items():
		df[k]["Variants"] = frozenset([x for x in v.keys() if x.startswith("c.") and tools.is_numeric(v[x]) and float(v[x]) > 0])
	if binary != None and len(binary) > 0:
		data = {k:{} for k in set([df[k][binary] for k in df.keys()])}
	else:
		data = {}
	if binary != None and len(binary) > 0:
		haplotypes = {status:set([df[k]["Variants"] for k,v in df.items() if df[k][binary] == status]) for status in data}
	else:
		haplotypes = set([df[k]["Variants"] for k,v in df.items()])
	return(haplotypes)

def haplo_patients(dffile, genename, response="Disease", binary=None):
	df = pd.read_csv(dffile, sep="\t", index_col=0)
	if binary != None:
		data = {k:{} for k in set(list(df[binary]))}
	else:
		data = {}
	for i, row in df.iterrows():
		variants = []
		for col in df.columns:
			if str(col).startswith("c.") and row[col] > 0:
				variants.append(col)
		variants = tuple(sorted(variants, key=lambda kv: float(max(re.findall("\d+", str(kv)), key=len))+nt_to_frac[re.findall("[ACGT]", str(kv))[1]]))
		if binary != None:
			if variants not in data[row[binary]]:
				data[row[binary]][variants] = [i]
			else:
				data[row[binary]][variants].append(i)
		else:
			if variants not in data:
				data[variants] = [i]
			else:
				data[variants].append(i)
	if binary != None:
		for k in data:
			newdf = os.path.join(os.path.dirname(dffile), genename, response + "_" + str(k) + "_haplos.tsv")
			tmpdf = {", ".join(x):{"Patients":", ".join(y)} for x,y in data[k].items()}
			pd.DataFrame(tmpdf).T.to_csv(newdf, sep="\t")
	else:
		newdf = os.path.join(os.path.dirname(dffile), genename, response + "_haplos.tsv")
		tmpdf = {", ".join(x):{"Patients":", ".join(y)} for x,y in data.items()}
		pd.DataFrame(tmpdf).T.to_csv(newdf, sep="\t")
	return(data)

def unique_vars(df, binary=None):
	df = pd.read_csv(df, sep="\t", index_col=0)
	if binary != None and len(binary) > 0:
		data = {k:[] for k in set(list(df[binary]))}
	else:
		data = []
	for i, row in df.iterrows():
		if binary != None and len(binary) > 0:
			data[row[binary]] += [col for col in df.columns if str(col).startswith("c.") and row[col] > 0]
			data[row[binary]]  = list(set(data[row[binary]]))
		else:
			data += [col for col in df.columns if str(col).startswith("c.") and row[col] > 0]
			data = list(set(data))
	if binary != None and len(binary) > 0:
		for k in data:
			data[k] = list(set(data[k]))
			print(str(k) + "\t" + str(len(data[k])))
	else:
		print("All:" + "\t" + str(len(data)))
	return(data)

def compare_sheets(f, response="Disease"):
	wb = openpyxl.load_workbook(f, read_only=True, keep_links=False)
	sheetnames = list(wb.sheetnames)
	sheets = {}
	for count in range(len(sheetnames)):
		try:
			sheets[count] = pd.read_excel(f, sheet_name=count, header=0, index_col=0)
		except Exception as e:
			print(e)
	data = {}
	for i in range(len(sheetnames)):
		for j in range(i):
			for col in response:
				k = [ki for ki in range(len(list(sheets[i][col]))) if tools.is_numeric(list(sheets[i][col])[ki]) and not np.isnan(list(sheets[i][col])[ki]) and tools.is_numeric(list(sheets[j][col])[ki]) and not np.isnan(list(sheets[j][col])[ki])]
				x = [list(sheets[i][col])[ki] for ki in k]
				y = [list(sheets[j][col])[ki] for ki in k]
				xisnum = [list(sheets[i][col])[ki] for ki in range(len(list(sheets[i][col]))) if tools.is_numeric(list(sheets[i][col])[ki]) and not np.isnan(list(sheets[i][col])[ki])]
				yisnum = [list(sheets[j][col])[ki] for ki in range(len(list(sheets[j][col]))) if tools.is_numeric(list(sheets[j][col])[ki]) and not np.isnan(list(sheets[j][col])[ki])]
				data[(sheetnames[j], sheetnames[i])] = [scipy.stats.wilcoxon(x,y), statistics.median(yisnum), statistics.median(xisnum)]
	return(data)

def rank_rows(f, response="Disease"):
	bestcolumns = ["Polyphen","SIFT", "MSA loglikelihood (AA)","PROVEAN"]
	wb = openpyxl.load_workbook(f, read_only=True, keep_links=False)
	sheetnames = list(wb.sheetnames)
	sheets = {}
	for count in range(len(sheetnames)):
		try:
			sheets[count] = pd.read_excel(f, sheet_name=count, header=0, index_col=0)
		except Exception as e:
			print(e)
	data = {}
	scores = {}
	sheet_scores = {}
	best_sheet_scores = {}
	for col in response:
		data[col] = {}
		for i in range(len(sheetnames)):
			if sheetnames[i] not in sheet_scores:
				sheet_scores[sheetnames[i]] = []
			if sheetnames[i] not in best_sheet_scores:
				best_sheet_scores[sheetnames[i]] = []
			try:
				data[col][sheetnames[i]] = sorted({list(sheets[i].index)[j]:list(sheets[i][col])[j] for j in range(len(list(sheets[i][col]))) if tools.is_numeric(list(sheets[i][col])[j]) and not np.isnan(list(sheets[i][col])[j])}.items(), key=lambda kv:kv[1])
				for row in list(sheets[i].index):
					if not re.search("[ACGT]\d+[ACGT]", str(row)):
						row = str(row)
						if row not in scores:
							scores[row] = []
						scores[row] += [j[1] for j in data[col][sheetnames[i]] if j[0] == row]
						sheet_scores[sheetnames[i]] += [j[1] for j in data[col][sheetnames[i]] if j[0] == row]
						if str(row) in bestcolumns:
							best_sheet_scores[sheetnames[i]] += [j[1] for j in data[col][sheetnames[i]] if j[0] == row]
			except:
				pass
	scores = {row:statistics.mean(s) for row, s in scores.items() if len(s) > 0}
	scores = sorted(scores.items(), key=lambda kv: kv[1])
	scores = {k:v for (k,v) in scores}
	sheet_scores = {row:statistics.mean(s) for row, s in sheet_scores.items() if len(s) > 0}
	sheet_scores = sorted(sheet_scores.items(), key=lambda kv: kv[1])
	sheet_scores = {k:v for (k,v) in sheet_scores}
	best_sheet_scores = {row:statistics.mean(s) for row, s in best_sheet_scores.items() if len(s) > 0}
	best_sheet_scores = sorted(best_sheet_scores.items(), key=lambda kv: kv[1])
	best_sheet_scores = {k:v for (k,v) in best_sheet_scores}
	return(data, scores, sheet_scores, best_sheet_scores)

def sum_rows(f, response="Disease"):
	wb = openpyxl.load_workbook(f, read_only=True, keep_links=False)
	sheetnames = list(wb.sheetnames)
	sheets = {}
	for count in range(len(sheetnames)):
		try:
			sheets[count] = pd.read_excel(f, sheet_name=count, header=0, index_col=0)
		except Exception as e:
			print(e)
	data = {}
	scores = {}
	for col in response:
		data[col] = {}
		for i in range(len(sheetnames)):
			data[col][sheetnames[i]] = {list(sheets[i].index)[j]:list(sheets[i][col])[j] for j in range(len(list(sheets[i][col]))) if tools.is_numeric(list(sheets[i][col])[j]) and not np.isnan(list(sheets[i][col])[j])}
	for col in data:
		scores[col] = {}
		for sheet in data[col]:
			for row in data[col][sheet]:
				if row not in scores[col]:
					scores[col][row] = [data[col][sheet][row]]
				else:
					scores[col][row].append(data[col][sheet][row])
		for row in scores[col]:
			scores[col][row] = statistics.mean(scores[col][row])
		scores[col] = statistics.mean(scores[col].values())
	return(scores)

def print_rows(data):
	for k1,v1 in data.items():
		print(k1)
		for i in range(len(v1[list(v1.keys())[0]])):
			for k2,v2 in v1.items():
				try:
					print(v2[i][0], end="\t")
				except:
					print("", end="\t")
			print("")

def merge_datasets(dirs, f=None, response="Disease"):
	if f == None or len(f) < 1:
		f = response + "_features_summarized.xlsx"
	columns = ["Polyphen","SIFT","PROVEAN", "PROVEAN (NT)", "GnomAD","dbSNP","Conservation (with PDB)","AA percent identity","AA entropy", "AA variance", "BLOSUM conservation", "Rare codon enrichment","Consurf (wo 3D info)","Accessible surface area (PDBPisa)", "Solvation energy (PDBPisa)", "NT percent identity","NT entropy","NT variance","Relative surface area (NetsurfP2)", "MSA loglikelihood (AA)","MSA loglikelihood (NT)","EVmutation (AA)","EVmutation (NT)","Number of mutants", "Delta Small", "Delta Charge", "Delta Aromatic", "Delta Large", "Delta Hydrophobicity", "Transition", "Transversion", "Delta Polar Neutral"]
	bestcolumns = ["Polyphen","SIFT", "MSA loglikelihood (AA)","PROVEAN"]
	data = {"overall AUCs":{}, "models by AUC":{}, "overall p-values":{}, "models by p-value":{}, "models (best) by p-value":{}, "models (best) by AUC":{}}
	for d in dirs:
		if os.path.isfile(d) and not os.path.isdir(d):
			d = os.path.dirname(d)
		for k in ["overall AUCs", "models by AUC", "overall p-values", "models by p-value", "models (best) by p-value", "models (best) by AUC"]:
			df = pd.read_excel(os.path.join(d, f), sheet_name=k, index_col=0)
			data[k][d] = df.to_dict()
	for d in data["overall p-values"]:
		if "Median p-value" not in data["overall p-values"][d]:
			print(d)
	newdata = {"overall AUCs":{col:{d:data["overall AUCs"][d][col]["Median AUC"] for d in data["overall AUCs"] if col in data["overall AUCs"][d]} for col in columns}, "models by AUC":{}, "overall p-values":{col:{d:data["overall p-values"][d]["Median p-value"][col] for d in data["overall p-values"] if col in data["overall p-values"][d]["Median p-value"]} for col in columns}, "models by p-value":{}, "models (best) by p-value":{}, "models (best) by AUC":{}}
	for col in ["addind", "addld", "maxind", "maxld"]:
		try:
			newdata["models by AUC"][col] = {}
			newdata["models (best) by AUC"][col] = {}
			for d in data["models by AUC"].keys():
				newdata["models by AUC"][col][d] = data["models by AUC"][d][col]["Median AUC"]
			for d in data["models (best) by AUC"].keys():
				newdata["models (best) by AUC"][col][d] = data["models (best) by AUC"][d][col]["Median AUC"]
		except KeyError as e:
			print(e)
		try:
			newdata["models by p-value"][col] = {}
			newdata["models (best) by p-value"][col] = {}
			for d in data["models by p-value"].keys():
				newdata["models by p-value"][col][d] = data["models by p-value"][d][col]["Median p-value"]
			for d in data["models (best) by p-value"].keys():
				newdata["models (best) by p-value"][col][d] = data["models (best) by p-value"][d][col]["Median p-value"]
		except KeyError as e:
			print(e)
			
	writer = pd.ExcelWriter(os.path.join("./",response + "_merged_datasets.xlsx"), engine='xlsxwriter')
	for k,v in newdata.items():
		pd.DataFrame(v).to_excel(writer, sheet_name=k)
	writer.save()
	return(newdata)

def plot_roc(prefix, tmpx, tmpy, outdir=None, step=0.001, invert=False, plot=True):
	if outdir == None or len(outdir) < 1:
		outdir = "./"

	x = [tmpx[i] for i in range(len(tmpx)) if tools.is_numeric(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpx[i]) and not np.isnan(tmpy[i])]
	if invert:
		x = [-xi for xi in x]
	y = [tmpy[i] for i in range(len(tmpx)) if tools.is_numeric(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpx[i]) and not np.isnan(tmpy[i])]
	'''
	#fpr, tpr, thresholds = sklearn.metrics.roc_curve(y,x)
	fpr = [0.0]
	tpr = [0.0]
	for ix in np.arange(min(x), max(x), step=step*(max(x)-min(x))):
		predicted = []
		for i  in range(len(x)):
			if x[i] <= ix:
				predicted += [1]
			else:
				predicted += [0]
		tp = sum([1 for i in range(len(y)) if predicted[i] == 1 and y[i] == 1])
		tn = sum([1 for i in range(len(y)) if predicted[i] == 0 and y[i] == 0])
		fn = sum([1 for i in range(len(y)) if predicted[i] == 0 and y[i] == 1])
		fp = sum([1 for i in range(len(y)) if predicted[i] == 1 and y[i] == 0])
		try:
			fpri = fp/(fp+tn)
			tpri = tp/(tp+fn)
			if tools.is_numeric(fpri) and tools.is_numeric(tpri) and not np.isnan(fpri) and not np.isnan(tpri):
				fpr += [fpri]
				tpr += [tpri]
		except ZeroDivisionError as e:
			pass
		except Exception as e:
			print(e)
	fpr += [1.0]
	tpr += [1.0]

	total = sum([(fpr[i]-fpr[i-1])*(1/2)*(tpr[i]+tpr[i-1]) for i in range(1, len(tpr))])
	'''
	plt.rcParams.update({'font.size': 13})
	fpr, tpr, thresholds = sklearn.metrics.roc_curve(y, x)
	total = sklearn.metrics.roc_auc_score(y,x)
	if plot:
		plt.cla()
		plt.clf()
		plt.plot(fpr, tpr)
		plt.xlabel("False positive rate")
		plt.ylabel("True positive rate")	
		plt.text(0.7,0,"AUC="+str(round(total,3)))
		plt.savefig(os.path.join(outdir, prefix + "_ROC" + ".png"))
	return(total)

def fig3(spreadsheet="merged_output_nodups.xlsx", sheets=["addind", "maxind"], xcols={"SIFT": True, "Polyphen": True, "MSA loglikelihood (AA)": True, "PROVEAN": True}, dim=(2,2), ycol="Disease", outdir="./"):
	colors = ['b', 'g', 'r', 'c', 'm', 'y']
	sheetdict = {"maxind":"Maximum independent", "addind":"Additive independent", "maxld":"Maximum paired", "addld":"Additive paired"}
	plt.rcParams.update({'font.size': 11})
	plt.cla()
	plt.clf()
	if dim[0] * dim[1] != len(xcols):
		raise Exception("Dimension of figure does not match number of subplots")
	if outdir == None or len(outdir) < 1:
		outdir = "./"
	for plti, (xcol, invert) in enumerate(xcols.items()):
		ax = plt.subplot(dim[0], dim[1], plti+1)
		auc = []
		for pltj, sheet in enumerate(sheets):
			df = pd.read_excel(spreadsheet, sheet_name=sheet, index_col=0)
			tmpx = list(df[xcol])
			tmpy = list(df[ycol])
			x = [tmpx[i] for i in range(len(tmpx)) if tools.is_numeric(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpx[i]) and not np.isnan(tmpy[i])]
			if invert:
				x = [-xi for xi in x]
			y = [tmpy[i] for i in range(len(tmpx)) if tools.is_numeric(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpx[i]) and not np.isnan(tmpy[i])]
			fpr, tpr, thresholds = sklearn.metrics.roc_curve(y, x)
			total = sklearn.metrics.roc_auc_score(y,x)
			auc.append(total)
			#print(sheet + "\t" + xcol + "\t" + str(total))
			ax.plot(fpr, tpr, colors[pltj])
		#ax.xlabel("False positive rate")
		#ax.ylabel("True positive rate")	
		if xcol == "Polyphen":
			xcol = "Polyphen-2"
		ax.set_title(xcol)
		plt.legend(handles=[mpatches.Patch(color=colors[pltj], label=sheetdict[sheets[pltj]] + "\n" + "AUC: " + str(round(auc[pltj], 3))) for pltj in range(len(sheets))], fontsize=8)
	plt.subplots_adjust(wspace=0.37, hspace=0.37)
	plt.savefig(os.path.join(outdir, "Fig3" + ".tiff"), dpi=300)

def graphic_abstract(spreadsheet="merged_output_nodups.xlsx", sheets=["addind", "maxind"], xcols={"SIFT": True, "Polyphen": True, "MSA loglikelihood (AA)": True, "PROVEAN": True}, dim=(2,2), ycol="Disease", randauc="AUC.pkl", bas = "BAS.pkl", outdir="./"):
	randauc = tools.pickleload(randauc)
	randauc = {k2:[x for k in randauc for x in randauc[k][k2]] for k2 in randauc["maxind"]}
	bas = tools.pickleload(bas)
	bas = {k2:[x for k in bas for x in bas[k][k2]] for k2 in bas["maxind"]}
	colors = ['b', 'g', 'r', 'c', 'm', 'y']
	sheetdict = {"maxind":"Maximum independent", "addind":"Additive independent", "maxld":"Maximum paired", "addld":"Additive paired"}
	plt.rcParams.update({'font.size': 11})
	plt.rcParams['font.family'] = "Arial"
	plt.cla()
	plt.clf()
	if dim[0] * dim[1] != len(xcols):
		raise Exception("Dimension of figure does not match number of subplots")
	if outdir == None or len(outdir) < 1:
		outdir = "./"
	for plti, (xcol, invert) in enumerate(xcols.items()):
		ax = plt.subplot(dim[0], dim[1], plti+1)
		c = 'c'
		bp = ax.boxplot(randauc[xcol], positions=[0.05], patch_artist=True, boxprops={'facecolor':c, 'color':c}, flierprops = dict(marker='o', markerfacecolor='w', markersize=4, linestyle='none', markeredgecolor=c, alpha=0.5))
		for element in ['boxes', 'whiskers', 'means', 'caps']:
			plt.setp(bp[element], color=c)
		plt.setp(bp["medians"], color='k')
		c = 'm'
		bp = ax.boxplot(bas[xcol], positions=[0.95], patch_artist=True, boxprops={'facecolor':c, 'color':c}, flierprops = dict(marker='o', markerfacecolor='w', markersize=4, linestyle='none', markeredgecolor=c, alpha=0.5))
		for element in ['boxes', 'whiskers', 'means', 'caps']:
			plt.setp(bp[element], color=c)
		plt.setp(bp["medians"], color='k')
		
		auc = []
		for pltj, sheet in enumerate(sheets):
			df = pd.read_excel(spreadsheet, sheet_name=sheet, index_col=0)
			tmpx = list(df[xcol])
			tmpy = list(df[ycol])
			x = [tmpx[i] for i in range(len(tmpx)) if tools.is_numeric(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpx[i]) and not np.isnan(tmpy[i])]
			if invert:
				x = [-xi for xi in x]
			y = [tmpy[i] for i in range(len(tmpx)) if tools.is_numeric(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpx[i]) and not np.isnan(tmpy[i])]
			fpr, tpr, thresholds = sklearn.metrics.roc_curve(y, x)
			total = sklearn.metrics.roc_auc_score(y,x)
			auc.append(total)
			#print(sheet + "\t" + xcol + "\t" + str(total))
			ax.plot(fpr, tpr, colors[pltj])
		ax.set_xlim(-0.05, 1.05)
		tickl = [0.0, 0.25, 0.5, 0.75, 1.0]
		ax.set_xticks(tickl)
		ax.set_xticklabels([str(ticki) for ticki in tickl])
		if xcol == "Polyphen":
			xcol = "Polyphen-2"
		if dim == (1,1) or len(dim) < 2:
			title = "ROC of disease prediction of simultaneous multiple variants"
			ax.text(0.72, 0.33, "Predictor: " + xcol)
		else:
			title = xcol
		#ax.set_title(title)
		#ax.xlim(0.0, 1.0)

		plt.legend(handles=[mpatches.Patch(color=colors[pltj], label=sheetdict[sheets[pltj]] + "\n" + "AUC: " + tools.padstring(str(round(auc[pltj], 3)), n=5, space='0')) for pltj in range(len(sheets))]+[mpatches.Patch(color='c', label='AUC distribution'),  mpatches.Patch(color='m', label='BAS distribution')], fontsize=10)
	#plt.legend(handles=[mpatches.Patch(color='c', label="AUC"), mpatches.Patch(color='m', label="BAS")], fontsize=8)
	plt.subplots_adjust(wspace=0.37, hspace=0.37)
	plt.savefig(os.path.join(outdir, "graphic_abstract" + ".tiff"), figsize=(1920/2.0/300, 1440/2.0/300), dpi=300)

def compute_mse(f, sheet, xcol, ycol, error="p", verbose=False):
	df = pd.read_excel(f, sheet_name=sheet, index_col=0)
	tmpx = list(df[xcol])
	tmpy = list(df[ycol])
	x = [tmpx[i] for i in range(len(tmpx)) if tools.is_numeric(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpx[i]) and not np.isnan(tmpy[i])]
	y = [tmpy[i] for i in range(len(tmpx)) if tools.is_numeric(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpx[i]) and not np.isnan(tmpy[i])]
	model = sm.OLS(y,[[x[i]] for i in range(len(x))]).fit()
	predictions = model.predict([[x[i]] for i in range(len(x))])
	if verbose:
		print(str(model.summary()))
	return(model)

#identify the model name within a string
def find_substr_opt(s, options=["maxind", "maxLD", "addind", "addLD"]):
	for opt in options:
		if opt.lower() in s.lower():
			return(opt)
	return("")

def test_columns_mse(f, ycol, xcols = {"Polyphen":False, "SIFT":False, "PROVEAN": False, "PROVEAN (NT)":False, "GnomAD":False, "dbSNP":False, "Conservation (with PDB)":False, "AA percent identity":False, "AA entropy":True, "AA variance":True, "BLOSUM conservation":True, "Rare codon enrichment":True, "Consurf (wo 3D info)":False, "Accessible surface area (PDBPisa)":False, "Solvation energy (PDBPisa)":False, "NT percent identity":False, "NT entropy":True, "NT variance":True, "Relative surface area (NetsurfP2)":False, "MSA loglikelihood (AA)":False, "MSA loglikelihood (NT)":False, "EVmutation (AA)":True, "EVmutation (NT)":True}, error="p-value", outdir=None, verbose=False):
	if outdir == None or len(outdir) < 1:
		outdir = os.path.dirname(f)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	wb = openpyxl.load_workbook(f, read_only=True, keep_links=False)
	sheetnames = list(wb.sheetnames)

	writer = pd.ExcelWriter(os.path.join(outdir,ycol + "_" + error.upper() + "_results.xlsx"), engine='xlsxwriter')
	for sheet in sheetnames:
		data = {}
		df = pd.read_excel(f, sheet_name=sheet, index_col=0)
		for col, invert in xcols.items():
			try:
				if error.lower() in ["p", "pvalue", "p-value"]:
					data[col] = compute_mse(f, sheet, col, ycol, error=error, verbose=verbose).f_pvalue
				else:
					data[col] = compute_mse(f, sheet, col, ycol, error=error, verbose=verbose).mse_resid
			except Exception as e:
				print(e)
		data["Mean"] = statistics.mean([k for k in data.values() if tools.is_numeric(k) and not np.isnan(k)])
		pd.DataFrame({k:{error.upper():v} for k,v in data.items()}).T.to_excel(writer, sheet_name=find_substr_opt(sheet))
	writer.save()

#compute the AUC for columns
def test_columns_roc(f, ycol, xcols={"Polyphen":False, "SIFT":False, "PROVEAN": False, "PROVEAN (NT)":False, "GnomAD":False, "dbSNP":False, "Conservation (with PDB)":False, "AA percent identity":False, "AA entropy":True, "AA variance":True, "BLOSUM conservation":True, "Rare codon enrichment":True, "Consurf (wo 3D info)":False, "Accessible surface area (PDBPisa)":False, "Solvation energy (PDBPisa)":False, "NT percent identity":False, "NT entropy":True, "NT variance":True, "Relative surface area (NetsurfP2)":False, "MSA loglikelihood (AA)":False, "MSA loglikelihood (NT)":False, "EVmutation (AA)":True, "EVmutation (NT)":True}, plot=True, outdir=None, step=0.001):
	if outdir == None or len(outdir) < 1:
		outdir = os.path.dirname(f)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	wb = openpyxl.load_workbook(f, read_only=True, keep_links=False)
	sheetnames = list(wb.sheetnames)

	writer = pd.ExcelWriter(os.path.join(outdir,ycol + "_AUC_results.xlsx"), engine='xlsxwriter')
	for sheet in sheetnames:
		data = {}
		df = pd.read_excel(f, sheet_name=sheet, index_col=0)
		for col, invert in xcols.items():
			try:
				if not os.path.exists(os.path.join(outdir, find_substr_opt(sheet))):
					os.makedirs(os.path.join(outdir, find_substr_opt(sheet)))
				data[col] = plot_roc(f, sheet, col, ycol, outdir=os.path.join(outdir, find_substr_opt(sheet)), step=step, invert=invert, plot=plot)
			except Exception as e:
				print(e)
		data["Mean"] = statistics.mean(data.values())
		pd.DataFrame({k:{"AUC":v} for k,v in data.items()}).T.to_excel(writer, sheet_name=find_substr_opt(sheet))
	writer.save()

def plot_rocs_fig(f, ycol, xcols={"Polyphen":True, "SIFT":True, "PROVEAN": True, "GnomAD":True, "Conservation (with PDB)":True, "AA conservation":False, "AA entropy":True, "AA variance":True, "BLOSUM conservation":True, "Rare codon enrichment":True, "Consurf (wo 3D info)":True, "Accessible surface area (PDBPisa)":False, "Solvation energy (PDBPisa)":False, "NT conservation":False, "NT entropy":True, "NT variance":True, "Relative surface area (NetsurfP2)":False, "MSA loglikelihood (AA)":True, "MSA loglikelihood (NT)":True, "EVmutation (AA)":True}, plot=True, outdir=None, step=0.001):
	plt.rcParams.update({'font.size': 3.5})
	colormap = {0:"b", 1:"g", 2:"r", 3:"y", 4:"m", 5:"c"}
	if outdir == None or len(outdir) < 1:
		outdir = os.path.dirname(f)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	wb = openpyxl.load_workbook(f, read_only=True, keep_links=False)
	sheetnames = list(wb.sheetnames)

	writer = pd.ExcelWriter(os.path.join(outdir,ycol + "_AUC_results.xlsx"), engine='xlsxwriter')
	fig, axs = plt.subplots(len(xcols), 1, figsize=(1.5, 50*len(xcols)), dpi=300)
	for i, (col, invert) in enumerate(xcols.items()):
		data = {}
		ax = axs[i]
		ax.plot([i*0.01 for i in range(100)], [i*0.01 for i in range(100)], c="black", alpha=0.5)
		for j, sheet in enumerate(sheetnames):
			df = pd.read_excel(f, sheet_name=sheet, index_col=0)
			try:
				if not os.path.exists(os.path.join(outdir, find_substr_opt(sheet))):
					os.makedirs(os.path.join(outdir, find_substr_opt(sheet)))
				data[col] = plot_roc_2(find_substr_opt(sheet) + " " + col, list(df[col]), list(df[ycol]), ax, c=colormap[j], outdir=outdir, step=step, invert=invert, plot=plot)
			except Exception as e:
				print(e)
				traceback.print_exc()
		#data["Mean"] = statistics.mean(data.values())
		#pd.DataFrame({k:{"AUC":v} for k,v in data.items()}).T.to_excel(writer, sheet_name=find_substr_opt(sheet))
	plt.xlabel("False positive rate")
	plt.ylabel("True positive rate")
	plt.legend(handles=[mpatches.Patch(color=colormap[i], label=sheetnames[i]) for i in range(len(sheetnames))])
	plt.show()
	writer.save()

def plot_roc_2(prefix, tmpx, tmpy, ax, c="b", outdir=None, step=0.001, invert=False, plot=True):
	if outdir == None or len(outdir) < 1:
		outdir = "./"

	x = [tmpx[i] for i in range(len(tmpx)) if tools.is_numeric(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpx[i]) and not np.isnan(tmpy[i])]
	if invert:
		x = [-xi for xi in x]
	y = [tmpy[i] for i in range(len(tmpx)) if tools.is_numeric(tmpx[i]) and tools.is_numeric(tmpy[i]) and not np.isnan(tmpx[i]) and not np.isnan(tmpy[i])]
	fpr, tpr, thresholds = sklearn.metrics.roc_curve(y, x)
	total = sklearn.metrics.roc_auc_score(y,x)
	if plot:
		ax.plot(fpr, tpr, color=c)
		ax.set_title(prefix, y=0.02, x=0.5, alpha=0.8)
	return(total)

def gen_p_values(f, ycol="Disease", xcols={"Polyphen":True, "SIFT":True, "PROVEAN": True, "GnomAD":True, "Conservation (with PDB)":True, "AA conservation":False, "AA entropy":True, "AA variance":True, "BLOSUM conservation":True, "Rare codon enrichment":True, "Consurf (wo 3D info)":True, "Accessible surface area (PDBPisa)":False, "Solvation energy (PDBPisa)":False, "NT conservation":False, "NT entropy":True, "NT variance":True, "Relative surface area (NetsurfP2)":False, "MSA loglikelihood (AA)":True, "MSA loglikelihood (NT)":True, "EVmutation (AA)":True}, plot=True, outdir=None):
	if outdir == None or len(outdir) < 1:
		outdir = os.path.dirname(f)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	wb = openpyxl.load_workbook(f, read_only=True, keep_links=False)
	sheetnames = list(wb.sheetnames)
	data = {}
	for sheet in sheetnames:
		df = pd.read_excel(f, sheet_name=sheet, index_col=0)
		if xcols == None or len(xcols) < 1:
			xcols = [xi for xi in df.columns if xi != ycol]
		k1 = find_substr_opt(sheet)
		data[k1] = {}
		yset = set([yi for yi in list(df[ycol])])
		for xcol in xcols:
			if ycol in df.columns:
				x = {yi:[row[xcol] for i, row in df.iterrows() if tools.is_numeric(row[xcol]) and not np.isnan(row[xcol]) and row[ycol] == yi] for yi in yset}
			else:
				x = {}
				x[0] = [row[xcol] for i, row in df.iterrows() if tools.is_numeric(row[xcol]) and not np.isnan(row[xcol]) and i.startwith("Neutral")]
				x[1] = [row[xcol] for i, row in df.iterrows() if tools.is_numeric(row[xcol]) and not np.isnan(row[xcol]) and not i.startwith("Neutral")]
			if len(x.keys()) == 2:
				try:
					data[k1][xcol] = scipy.stats.mannwhitneyu(x[list(x.keys())[0]], x[list(x.keys())[1]])
				except Exception as e:
					pass
	data = {k2:{k1:v2 for k1 in data.keys()} for k2,v2 in data[list(data.keys())[0]].items()} 
	data = {k:max(v.values(), key=lambda kv: kv.pvalue) for k,v in data.items()}
	reject, pvals_corrected, alphaSidak, alphaBonf = statsmodels.stats.multitest.multipletests([data[k].pvalue for k in data.keys()], method="bonferroni")
	data = {list(data.keys())[i]:{"Reject":reject[i], "p-value (Bonferroni)":pvals_corrected[i]} for i in range(len(data.keys()))}
	return(data)

def gen_p_values_2(f, ycol="Disease", dirs=["./LOVD/F8/", "./LOVD/F9/", "./LOVD/DMD/", "./LOVD/HBB/", "./LOVD/VWF/VWD1/", "./LOVD/VWF/VWD2/", "./LOVD/VWF/VWD3/", "./LOVD/G6PD/", "./F8_F9/F8/", "./F8_F9/F9/", "./cTTP/ADAMTS13/", "./VWD/VWF/"], most_recent=True):
	data = {}
	for d in dirs:
		if most_recent:
			files = [fi for fi in os.listdir(d) if re.match(f, fi)]
			inf = max(files,key=lambda kv:int(max(re.findall("\d+", kv), key=len)))
		else:
			inf = f
		data[d] = gen_p_values(os.path.join(d,inf), ycol, xcols=None)
	return(data)

def combine_p_values(data):
	keys = [k for k in data[list(data.keys())[0]].keys()]
	d = {}
	for k2 in keys:
		vals = {k1:data[k1][k2]["p-value (Bonferroni)"] for k1 in data if k2 in data[k1] and "p-value (Bonferroni)" in data[k1][k2] and data[k1][k2]["p-value (Bonferroni)"] != 0.0}
		d[k2] = {"p-value (Fisher)":scipy.stats.combine_pvalues(list(vals.values()), method="fisher")[1], "p-value (Pearson)":scipy.stats.combine_pvalues(list(vals.values()), method="pearson")[1], "p-value (Stouffer)":scipy.stats.combine_pvalues(list(vals.values()), method="stouffer")[1]}
		try:
			d[k2]["p-value (Tippett)"] = scipy.stats.combine_pvalues(list(vals.values()), method="tippett")[1]
		except:
			pass
		try:
			d[k2]["p-value (Mudholkar/George)"] = scipy.stats.combine_pvalues(list(vals.values()), method="mudholkar_george")[1]
		except:
			pass
		for k,v in vals.items():
			d[k2][k] = v
	return(d)

def combine_files(directory):
	data = {}
	for f in os.listdir(directory):
		wb = openpyxl.load_workbook(os.path.join(directory, f), read_only=True, keep_links=False)
		sheetnames = list(wb.sheetnames)
		filename = f.replace("_AUC_results.xlsx", "")
		data[filename] = {}
		for sheet in sheetnames:
			model = [x for x in ["addind", "addLD", "maxLD", "maxind"] if x in str(sheet)][0]
			df = pd.read_excel(os.path.join(directory, f), sheet_name=sheet, index_col=0)
			df = df.to_dict()["AUC"]
			data[filename][model] = {}
			for key in df:
				data[filename][model][key] = df[key]
	return(data)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='This script reads all similarly named files spreadsheets in a directory and combines them into one excel file', add_help=True)
	parser.add_argument('directory', type=str, metavar='PATH', help="path to directory containing spreadsheets")
	parser.add_argument('regex', type=str, help="regex for filename")
	parser.add_argument('--sep', choices=["comma", "tab", "space"], default="tab")
	clargs = parser.parse_args()

	sep = {"comma":",", "tab":"\t", "space":" "}[clargs.sep]
	merge_files(clargs.directory, clargs.regex, sep=sep)
