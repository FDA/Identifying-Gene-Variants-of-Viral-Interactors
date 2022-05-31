import os, sys, re, subprocess, time, scipy.stats, copy, seaborn
import random
import numpy as np
import warnings
import matplotlib.pyplot as plt
from scipy.stats import *
from sklearn.ensemble import *
from sklearn.decomposition import PCA
from tools import dim, is_numeric, infer_numeric, update_time
from intervaltree import *
from statsmodels.stats.descriptivestats import sign_test
from statsmodels.stats.multitest import multipletests
import statistics
import math

def cdf(d):
	d = sorted(d.items(), key=lambda kv: kv[1])
	d = {d[i][0]:i/len(d) for i in range(len(d))}
	return(d)


def convert_to_m(d,sym=True):
	l1 = sorted(set([k[0] for k in d.keys()]))
	l2 = sorted(set([k[1] for k in d.keys()]))
	init_time = time.time()
	m = [[np.nan for j in range(min(l2), max(l2)+1)] for i in range(min(l1), max(l1)+1)]
	for i, k1 in enumerate(l1):
		for j, k2 in enumerate(l2):
			if (k1, k2) in d.keys():
				m[i][j] = d[(k1, k2)]
			elif sym and (k2, k1) in d.keys():
				m[i][j] = d[(k2, k1)]
		tools.update_time(i, len(l1), init_time)
	return(m)

def compare_dicts(d1, d2):
	alphabet = {}
	for k in d1.keys():
		for i in range(len(k)):
			if i not in alphabet:
				alphabet[i] = []
			if k[i] not in alphabet[i]:
				alphabet[i].append(k[i])

	print("", end="\t")
	for k in sorted(alphabet[0]):
		print(k, end="\t")
	print("")
	for i in sorted(alphabet.keys()):
		print(i+1, end="\t")
		for k in sorted(alphabet[i]):
			set1 = [d1[c] for c in d1 if c[i] == k]
			set2 = [d2[c] for c in d1 if c[i] == k]
			wc = wilcoxon(set1, set2)
			print(round(wc.pvalue, 5), end="\t")
		print("")
	wc = wilcoxon([d1[k] for k in d1], [d2[k] for k in d1])
	print("Wilcoxon p-value: " + str(round(wc.pvalue, 5)))

#computes entropy of an array
def shannon_entropy(arr):
	tot = sum(arr)
	new_arr = [arr[i]/tot for i in range(len(arr))]
	return(sum([-new_arr[i]*math.log(new_arr[i]) for i in range(len(new_arr)) if new_arr[i] > 0.0]))

def trans_mat(mat):
	return([[mat[i][j] for i in range(len(mat))] for j in range(len(mat[0]))])

def convert_to_freqs(d):
	tot = sum(d.values())
	d = copy.deepcopy({k:v/tot for k,v in d.items()})
	return(d)

#creates a sampling tree whose distribution is based on the difference target - actual
#designed to move toward the target distribution
def weighted_tree(d_target, d_actual, limit=0.001):
	for d in [d_target, d_actual]:
		tot = sum(d.values())
		if tot != 0:
			for c in d:
				d[c] = d[c]/tot
		else:
			for c in d:
				d[c] = 0

	d = {c:(d_target[c]-d_actual[c]) for c in d_target}
	min_d = min(d.values())
	if max(d.values()) - min_d > limit:
		d = {c:(d[c]-min_d) for c in d}
	else: #if range of differences is too small
		d = {c:d_target[c] for c in d_target}
	return(range_tree(d))

#uses only numeric values within x and y
#useful to avoid errors when computing statistics
def get_numeric(x,y=None):
	singleFlag = False
	if y == None:
		singleFlag = True
		y = copy.deepcopy(x)
	x = list(x)
	y = list(y)
	if not len(x) == len(y):
		raise Exception("Input lists have different lengths. Unable to resolve")
	newx = [x[i] for i in range(len(x)) if is_numeric(x[i]) and not np.isnan(x[i]) and is_numeric(y[i]) and not np.isnan(y[i])]
	newy = [y[i] for i in range(len(y)) if is_numeric(x[i]) and not np.isnan(x[i]) and is_numeric(y[i]) and not np.isnan(y[i])]
	if singleFlag:
		return(newx)
	else:
		return(newx, newy)

#generate a range tree from a dictionary of counts or probabilities to be used for sampling
def range_tree(d):
	tot = sum(d.values())
	d = {k:d[k]/tot for k in d if d[k] > 0} #convert all numbers to probabilities
	tree = IntervalTree()
	min_float = np.nextafter(0,1) #minimum positive number
	cum = 0
	for k,v in d.items():
		try:
			tree.addi(cum, cum+v-min_float, k) #add node to tree with appropriately long (unique) interval
			cum += v
		except ValueError as e:
			pass
	return(tree)

#sample from a range tree
def sample_from_tree(t):
	val = random.uniform(0,1)
	return(t[val].pop()[2]) #return the element having val in its interval

def run_mwu(df, xcol, ycol, rounding="UP"):
	x = list(df[xcol])
	y = list(df[ycol])
	if rounding.upper() == "UP":
		x = [math.ceil(xi) for xi in x]
	elif rounding.upper() =="DOWN":
		x = [math.floor(xi) for xi in x]
	if len(set(x)) != 2:
		raise Exception("x column is not a binary variable, even with rounding")
	z0 = list(set(x))[0]
	z1 = list(set(x))[1]
	y0 = [y[i] for i in range(len(x)) if x[i] == z0]
	y1 = [y[i] for i in range(len(x)) if x[i] == z1]
	return(scipy.stats.mannwhitneyu(y0, y1))
	
#calculate difference in explained variance with and without columns
def calc_explanatory_df(df, columns=[]):
	pca = PCA(n_components=1, random_state=0)
	df = df.fillna(df.median())
	pca.fit(df)
	fullmodel = pca.explained_variance_[0]

	df2 = df[df.columns.difference(columns)]
	partmodel = pca.fit(df).explained_variance_[0]
	return((fullmodel - partmodel)/fullmodel)
			
#calculates the reduction in the MSE due to each feature
def calc_explanatory(x, y, regr, n=10, direction="down"):
	F = []
	#for each feature/variable
	init_time = time.time()
	for i in range(n):

		if direction.lower() in ["down"]:
			train_x, train_y, train_weights = bootstrap(x, y, size=0.8)
			trans_x = trans_mat(train_x)
			regr.fit(train_x, train_y)
			mse_full = calc_mse(train_y, regr.predict(train_x))
			for j in range(dim(x)[1]):
				if len(F) < dim(x)[1]:
					F.append([])
				tmp_x = trans_mat(trans_x[:j]+trans_x[j+1:]) #remove the ith feature/variable
				regr.fit(tmp_x, train_y)
				mse_part = calc_mse(train_y, regr.predict(tmp_x))
				F[j].append((mse_part-mse_full)/mse_full)
		elif direction.lower() in ["up"]:
			train_x, train_y, train_weights = bootstrap(x, y, size=0.8)
			trans_x = trans_mat(train_x)
			mse_full = calc_mse(train_y, [0 for i in range(len(train_y))])
			for j in range(dim(x)[1]):
				if len(F) < dim(x)[1]:
					F.append([])
				tmp_x = trans_mat([trans_x[j]]) #include only the ith feature/column
				regr.fit(tmp_x, train_y)
				mse_part = calc_mse(train_y, regr.predict(tmp_x))
				F[j].append((mse_full-mse_part)/mse_full)
		elapsed = time.time() - init_time
		estimated_time = elapsed/(i+1)*(n-(i+1))
		print(str(round(100.0*(i+1)/n, 2)) + "% done. Estimated time remaining: " + str(int(estimated_time/60)) + " minutes, " + str(int(estimated_time % 60)) + " seconds.", end="\r", flush=True)
	for j in range(len(F)):
		F[j] = np.mean(F[j])
	return(F)

def calc_mse(arr1, arr2):
	if len(arr1) != len(arr2):
		raise Exception("Arrays have different lengths.")
	err = sum([(arr1[i]-arr2[i])**2 for i in range(len(arr1))])/len(arr1)
	return(err)

#computes n boostrap cross-validation
def cross_validate(regr, x, y, mode="regr", n=100):
	accuracies = []
	init_time = time.time()
	for i in range(n):
		train_x, train_y, train_weights = bootstrap(x, y, size=0.8)
		test_x, test_y, test_weights = bootstrap(x, y, size=0.2)
		regr.fit(train_x, train_y)
		predicted = regr.predict(test_x)
		if mode.lower() in ["reg", "regr", "regressor", "regression"]:
			accuracy = calc_mse(predicted, test_y)
		else:
			accuracy = sum([1 for j in range(len(test_y)) if test_y[j] == predicted[j]])/len(test_y)
		accuracies.append(accuracy)
		elapsed = time.time() - init_time
		estimated_time = elapsed/(i+1)*(n-(i+1))
		print("Accuracy: " + str(round(accuracy, 3)) + ", " + str(round(100.0*(i+1)/n, 2)) + "% done. Estimated time remaining: " + str(int(estimated_time/60)) + " minutes, " + str(int(estimated_time % 60)) + " seconds.", end="\r", flush=True)
	return(accuracies)

#generates a bootstrap sample
def bootstrap(x, y, weights=[], size=0.2):
	if weights == None or len(weights) < 1:
		weights = [1 for i in range(len(y))]

	indices = []
	for i in range(int(size*len(y))):
		indices.append(random.randint(0, len(y)-1))
	sample_x = [x[i] for i in indices]
	sample_y = [y[i] for i in indices]
	sample_weights = [weights[i] for i in indices]
	return(sample_x, sample_y, sample_weights)

def compare_datasets(datadict, dataref):
	for mode in datadict:
		pairs = {}
		for d in datadict[mode].keys():	
			pairs[d] = [[datadict[mode][d][c] for c in datadict[mode][d]], [dataref[d][c] for c in datadict[mode][d]]]
		print("\n" + mode)
		for name, pair in pairs.items():		   
			err = 100*sum([abs((pair[0][i]-pair[1][i]))/pair[1][i] for i in range(len(pair[0])) if pair[1][i] != 0])/len([pair[1][i] for i in range(len(pair[1])) if pair[1][i] != 0])
			cor = pearsonr(pair[0], pair[1])			 
			mw = mannwhitneyu(pair[0], pair[1])
			wc = wilcoxon(pair[0], pair[1])
			print(name, end=" \t")
			print(err, end=" \t")
			print(cor[0], end=" \t")
			print(cor[1], end=" \t")
			print(mw.pvalue, end=" \t")
			print(wc.pvalue)

#assumes angles are measured in degrees and set between -180 and 180
def convert_angle(l):
	if any([True for i in range(len(l)) if l[i] < -180 or l[i] < 180]): #not angles
		return(l)
	else:
		for i in range(l):
			if abs(l[i]) > abs(l[i] - -180): #closer to -180 than zero
				l[i] += 360
		return(l)

def compare_lists(l1, l2, forbidden=["", None, np.nan], numeric=True):
	if numeric and (infer_numeric(l1) and infer_numeric(l2)): #both lists are (mostly) numeric
		m1 = [float(k) for k in l1 if is_numeric(k) and not np.isnan(float(k)) and not k in forbidden]
		m2 = [float(k) for k in l2 if is_numeric(k) and not np.isnan(float(k)) and not k in forbidden]
		corrm1 = convert_angle(m1)
		corrm2 = convert_angle(m2)
		mw = scipy.stats.mannwhitneyu(corrm1, corrm2)
		medians = (statistics.median(corrm1), statistics.median(corrm2))
		return(mw, m1, m2, medians)
	else: #lists are not numeric
		categories = sorted(set([k for k in l1 + l2 if not k in forbidden]))
		m1 = {category:0 for category in categories}
		m2 = {category:0 for category in categories}
		for i in range(len(l1)):
			m1[l1[i]] += 1
		for i in range(len(l2)):
			m2[l2[i]] += 1
		m1 = [m1[category] for category in categories]
		m2 = [m2[category] for category in categories]
			
		scale = sum(m1)/sum(m2)
		chisq = scipy.stats.chisquare(m1, [scale*k for k in m2])
		maxnum = max([i for i in range(len(m1))], key=lambda j: abs(m1[j] - scale*m2[j]))
		return(chisq, {categories[i]:m1[i] for i in range(len(m1))}, {categories[i]:m2[i] for i in range(len(m2))}, categories[maxnum], m1[maxnum]/sum(m1) - m2[maxnum]/sum(m2))

def hist_lists(data):
	l = {k:[] for k in data.keys()}
	for k,c in data.items():
		total = sum(c.values())
		c = {str(k2):int(1000*v2/total) for k2,v2 in c.items()}
		for k2,v2 in c.items():
			l[k] += [k2 for i in range(v2)]
		plt.hist(sorted(l[k]), label=k, alpha=1.0/len(data))
	plt.ylabel("Per thousand")
	plt.legend()
	plt.show()

def process_lists_numerical(l1, l2):
	newl1 = [float(l1[i]) for i in range(len(l1)) if is_numeric(l1[i]) and is_numeric(l2[i]) and not np.isnan(float(l1[i])) and not np.isnan(float(l2[i]))]
	newl2 = [float(l2[i]) for i in range(len(l1)) if is_numeric(l1[i]) and is_numeric(l2[i]) and not np.isnan(float(l1[i])) and not np.isnan(float(l2[i]))]
	return(newl1, newl2)

#compare values of x and y columns
def compare_columns(x,y):
	x = list(x)
	y = list(y)
	new_x = [x[i] for i in range(len(x)) if is_numeric(x[i]) and is_numeric(y[i]) and not np.isnan(x[i]) and not np.isnan(y[i])]
	new_y = [y[i] for i in range(len(x)) if is_numeric(x[i]) and is_numeric(y[i]) and not np.isnan(x[i]) and not np.isnan(y[i])]
	st = sign_test([new_x[i]-new_y[i] for i in range(len(new_x))])
	wc = scipy.stats.wilcoxon(new_x, new_y)
	pct_dif = [(new_y[i] - new_x[i])/new_x[i] for i in range(len(new_x)) if new_x[i] != 0]
	print("Sign test " + str(st))
	print(str(wc))
	print(str(statistics.mean(pct_dif)) + "\t" + str(statistics.stdev(pct_dif)))

#splits x and y arrays into random training/testing split, assuming x and y have same length
def random_split(x, y, weights=[], testsize=0.2):
	if weights == None or len(weights) < 1:
		weights = [1 for i in range(len(y))]
	
#randomly choose a sample to be the test set
	test_i = random.sample(list(range(len(y))), int(len(y)*testsize))

	test_x = [x[i] for i in test_i]
	test_y = [y[i] for i in test_i]
	test_weights = [weights[i] for i in test_i]

	train_x = [x[i] for i in range(len(y)) if i not in test_i]
	train_y = [y[i] for i in range(len(y)) if i not in test_i]
	train_weights = [weights[i] for i in range(len(y)) if i not in test_i]
	
	return(train_x, train_y, test_x, test_y, train_weights, test_weights)

def correlation_lists(x, y, forbidden=["", None], threshold = 3, mode="sm"):
	x = list(x)
	y = list(y)
	x = [x[i] for i in range(len(x)) if x[i] not in forbidden and y[i] not in forbidden]
	y = [y[i] for i in range(len(x)) if x[i] not in forbidden and y[i] not in forbidden]
	if (infer_numeric(x) and len(set(y)) <= 2) or (infer_numeric(y) and len(set(x)) <= 2):
		u = [xi for xi in [x,y] if len(set(xi)) <= 2][0]
		v = [xi for xi in [x,y] if infer_numeric(xi) and xi != u][0]
		new_u = [float(u[i]) for i in range(len(u)) if is_numeric(u[i]) and is_numeric(v[i]) and not np.isnan(float(u[i])) and not np.isnan(float(v[i]))]
		new_v = [float(v[i]) for i in range(len(v)) if is_numeric(u[i]) and is_numeric(v[i]) and not np.isnan(float(u[i])) and not np.isnan(float(v[i]))]
		u = new_u
		v = new_v
		corr = scipy.stats.pointbiserialr(u,v)
		corr = {"Correlation":corr[0], "p-value":corr[1]}
		return(corr)
	elif infer_numeric(x) and infer_numeric(y) and len(set(x)) > threshold:
		new_x = [float(x[i]) for i in range(len(x)) if is_numeric(x[i]) and is_numeric(y[i]) and not np.isnan(float(x[i])) and not np.isnan(float(y[i]))]
		new_y = [float(y[i]) for i in range(len(x)) if is_numeric(x[i]) and is_numeric(y[i]) and not np.isnan(float(x[i])) and not np.isnan(float(y[i]))]
		x = new_x
		y = new_y
		if len(x) >= threshold and len(y) >= threshold:
			if mode == "infer":
				if scipy.stats.shapiro(x).pvalue < 0.05 or scipy.stats.shapiro(y).pvalue < 0.05:
					tmpmode="sm"
				else:
					tmpmode="pearson"
			else:
				tmpmode = mode
			if tmpmode.lower() in ["sm", "spearman"]:
				corr = scipy.stats.spearmanr(x, y)
				corr = {"Correlation":corr.correlation, "p-value":corr.pvalue}
			else:
				corr = scipy.stats.pearsonr(x,y)
				corr = {"Correlation":corr[0], "p-value":corr[1]}
			return(corr)
	elif not infer_numeric(x) and not infer_numeric(y) and len(set(x)) > threshold:
		if len(x) >= threshold and len(y) >= threshold:
			confusion_matrix = pd.crosstab(x,y)
			chi2 = ss.chi2_contingency(confusion_matrix)[0]
			n = confusion_matrix.sum().sum()
			phi2 = chi2/n
			r,k = confusion_matrix.shape
			phi2corr = max(0, phi2-((k-1)*(r-1))/(n-1))
			rcorr = r-((r-1)**2)/(n-1)
			kcorr = k-((k-1)**2)/(n-1)
			V = np.sqrt(phi2corr/min((kcorr-1),(rcorr-1))) #Cramer's V
			xlabels = set(x)
			ylabels = set(y)
			obs = {(xi, yi):len([i for i in range(len(x)) if x[i] == xi and y[i] == yi]) for xi in xlabels for yi in ylabels}
			exp = {(xi, yi):len(x)*len([i for i in range(len(x)) if x[i] == xi])*len(y)*len([i for i in range(len(y)) if y[i] == yi]) for xi in xlabels for yi in ylabels}
			cs2 = scipy.stats.chisquare(obs, exp)
			return({"Correlation":v, "p-value":cs2.p})
	return(None)		
		

def correlation_df(df, forbidden=["", None], threshold = 3, sigthresh=0.05, mode="sm"):
	cor = {}
	for i, col1 in enumerate(df.columns):
		for j in range(i):
			col2 = list(df.columns)[j]
			try:
				cor[(col1, col2)] = correlation_lists(df[col1], df[col2], mode=mode)
				if cor[(col1, col2)] == None:
					del cor[(col1, col2)]
			except:
				pass
	cor = {k:v for k,v in cor.items() if v != None}
	keys = [k for k in list(cor.keys()) if "p-value" in cor[k] and is_numeric(cor[k]["p-value"]) and not np.isnan(float(cor[k]["p-value"]))]
	returnval = multipletests([float(cor[k]["p-value"]) for k in keys])
	pvals_corrected = returnval[1]
	for i in range(len(keys)):
		cor[keys[i]]["p-value"] = pvals_corrected[i]
	cor = {k:v for k,v in cor.items() if k in keys and not np.isnan(v["p-value"]) and v["p-value"] < sigthresh}
	formatted_cor = {}
	for col1 in df.columns:
		formatted_cor[col1] = {}
		for col2 in df.columns:
			if col1 == col2:
				formatted_cor[col1][col2] = 1
			elif (col1, col2) in cor:
				formatted_cor[col1][col2] = cor[(col1, col2)]["Correlation"]
			elif (col2, col1) in cor:
				formatted_cor[col1][col2] = cor[(col2, col1)]["Correlation"]
	return(cor, formatted_cor)

def filter_corr(corr, analogues):
	for k,v in list(corr.items()):
		s1 = None
		s2 = None
		for i, l in enumerate(analogues):
			for j, x in enumerate(l):
				if x.lower() in str(k[0]).lower() and s1 == None:
					s1 = i
				if x.lower() in str(k[1]).lower() and s2 == None:
					s2 = i
				if s1 != None and s2 != None:
					break
			if s1 != None and s2 != None:
				break
		if s1 != None and s1 == s2:
			del(corr[k])
	return(corr)
