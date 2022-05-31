import argparse, re
import numpy as np
import os.path
import tools, compute_features

def compute_RSCU(codons, translation_file=None, quiet=True):

	if translation_file == None:
		translation = compute_features.translate
	else:
		translation = tools.two_column_csv(translation_file)
	try:
		codon_data = tools.parse_cut(codons)
	except:
		codon_data = tools.two_column_csv(codons)
	RSCU = {}
	W = {}

	for codon1 in codon_data.keys():
		count = 0
		numAA = 0
		max_val = 0.0
		for codon2 in codon_data.keys():
			if translation[codon1] == translation[codon2]:
				count += 1
				numAA += int(codon_data[codon2])
				if codon_data[codon2] > max_val:
					max_val = codon_data[codon2]
		try:
			RSCU[codon1] = float(codon_data[codon1])*count/numAA
		except ZeroDivisionError as e:
			RSCU[codon1] = np.NaN
	for codon1 in codon_data.keys():
		W[codon1] = RSCU[codon1]/max([RSCU[codon2] for codon2 in codon_data.keys() if translation[codon1] == translation[codon2]])

	if not quiet:
		with open(os.path.splitext(codons)[0] + "_RSCU.tsv", "w") as codon_out:
			for codon in RSCU.keys():
				codon_out.write(codon + "\t" + str(RSCU[codon]) + "\n")

		with open(os.path.splitext(codons)[0] + "_W.tsv", "w")as cp_out:
			for cp in W.keys():
				cp_out.write(cp + "\t" + str(W[cp]) + "\n")

	return(RSCU, W)
	

def compute_RSCPU(codonpairs, translation_file, quiet=True):
	translation = tools.two_column_csv(translation_file)
	cp_data = tools.two_column_csv(codonpairs)
	RSCPU = {}
	W_CP = {}

	for cp1 in cp_data.keys():
		count = 0
		numAAp = 0
		max_val = 0.0
		for cp2 in cp_data.keys():
			if translation[cp1[0:3]] == translation[cp2[0:3]] and translation[cp1[3:6]] == translation[cp2[3:6]]:
				count += 1
				numAAp += int(cp_data[cp2])
				if cp_data[cp2] > max_val:
					max_val = cp_data[cp2]
		try:
			RSCPU[cp1] = float(cp_data[cp1])*count/numAAp
		except ZeroDivisionError as e:
			RSCPU[cp1] = np.NaN

	for cp1 in cp_data.keys():
		W_CP[cp1] = RSCPU[cp1]/max([RSCPU[cp2] for cp2 in cp_data.keys() if translation[cp1[0:3]] == translation[cp2[0:3]] and translation[cp1[3:6]] == translation[cp2[3:6]]])

	if not quiet:
		with open(os.path.splitext(codonpairs)[0] + "_RSCPU.tsv", "w")as cp_out:
			for cp in RSCPU.keys():
				cp_out.write(cp + "\t" + str(RSCPU[cp]) + "\n")

		with open(os.path.splitext(codonpairs)[0] + "_W_CP.tsv", "w")as cp_out:
			for cp in W_CP.keys():
				cp_out.write(cp + "\t" + str(W_CP[cp]) + "\n")

	return(RSCPU, W_CP)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='This script will compute RSCU, RSCPU, W, and W (CP), based on codon and codon-pair usage data and a translation table.', add_help=True)
	parser.add_argument('codons', type=str, metavar='FILE', default="CUT.cut", help="location of the codon data")
	parser.add_argument('codonpairs', type=str, metavar='FILE', default="CPUT.tsv", help="location of the codon pairdata")
	parser.add_argument('translation', type=str, metavar='FILE', default="translation.csv", help="location of the translation table")
	parser.add_argument('--write', action=store_true, help="whether or not to write RSCU, RSCPU, and W information to the hard drive")
	args = parser.parse_args()

	compute_RSCU(args.codons, args.translation, quiet=args.write)
	compute_RSCPU(args.codonpairs, args.translation, quiet=args.write)
