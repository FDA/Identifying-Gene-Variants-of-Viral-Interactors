# David Holcomb
# 2018

import pandas as pd
from ev_model import CouplingsModel
import numpy as np
import ev_tools
import argparse, subprocess, re, csv
import os.path, sys
import align_muts_with_msa
sys.path.insert(0, "../")
import tools #used only for reading fasta file

#assumes plmc executable is accessible (in the PATH)
def run_plmc(model_params="", focus="", fasta="", eij_lambda=16.2, hi_lambda=0.01, iterations=500, stepsize=0.2, alphabet ="-ACDEFGHIKLMNPQRSTVWY", quiet=False):
	if fasta is "":
		raise RuntimeError("Cannot run PLMC without a multiple sequence alignment.")
	if focus is "":
		raise RuntimeError("Cannot run PLMC without a focus sequence.")
	if model_params == None or len(model_params) < 1:
		model_params = focus + ".model_params"
	if not quiet:
		print("Running PLMC on the MSA. This may take a while...")
	plmc_return = subprocess.Popen(["plmc", '-o', model_params, '-a', alphabet, '-c', os.path.join(os.path.dirname(model_params), focus) + ".ec", '-f', focus, '-le', str(eij_lambda), '-lh', str(hi_lambda), '-m', str(iterations), '-t', str(stepsize), '-g', fasta], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
	
	#nonzero status returned, (usually a segfault)
	if int(plmc_return) != 0:
		raise RuntimeError("Error (probably segfault) occurred while running PLMC")
	elif not quiet:
		print("PLMC is complete.")

#convert dataframe to a matrix/dictionary, assuming each line in the dataframe is one mutation
def dictize(df, alphabet="ACDEFGHIKLMNPQRSTVWY"):
	alphabet = [alphabet[i] for i in range(len(alphabet))]
	m = {}
	for i, row in df.iterrows():
		mut = row['mutant']
		val = float(row['effect_prediction_epistatic'])
		aa = str(row['subs'])
		pos = int(row['pos'])

		if pos not in m.keys():
			m[pos] = {}
			for key in alphabet:
				m[pos][key] = 0.0

		#add new value to array
		m[pos][aa] = val
	return m

def create_analysis(fasta="", focus="", model_params="", mutations="", skip_plmc=False, skip_align_muts=False, mode="matrix", alphabet=None, space="-", quiet=False):
	if alphabet == None or len(alphabet) < 1: #infer alphabet from sequences
		seqs = tools.read_fasta(fasta)
		alphabet = [seqs[k][i] for k in seqs for i in range(len(seqs[k]))]
		alphabet = "".join(set(alphabet))
		nonalpha = re.findall("(?!\w).", alphabet) #find non-alphanumeric (gap) characters
		nonalpha = list(set([x[i] for x in nonalpha for i in range(len(x))]))
		if space == None and len(nonalpha) < 1:
			space = "-"
		elif space not in nonalpha and len(nonalpha) > 0: #space character not specified
			space = nonalpha[0]
			nonalpha = nonalpha[1:]
		else:	#soace character specified
			nonalpha = [nonalpha[i] for i in range(len(nonalpha)) if nonalpha[i] != space]
		alphabet = "".join(sorted([alphabet[i] for i in range(len(alphabet)) if alphabet[i] not in nonalpha and alphabet[i] != space])) #remove "gap" characters
		alphabet = str(space) + alphabet + "".join(nonalpha)
	else: #remove sequences with characters outside of alphabet (may be ambiguity chars)
		space = alphabet[0]
		seqs = tools.read_fasta(fasta)
		for k,v in list(seqs.items()):
			if not re.match("[" + alphabet + "]+", v, re.IGNORECASE):
				del seqs[k]
	if model_params is "":
		model_params = focus + ".model_params"

	#run plmc
	if not skip_plmc:	
		run_plmc(model_params, focus, fasta, alphabet=alphabet, quiet=quiet)
	
	#align mutations spreadsheet with multiple sequence alignment to account for deletions and insertions ('-')
	if not skip_align_muts:
		align_muts_with_msa.align_with_msa(fasta, focus, mutations)
		mutations = os.path.join(os.path.dirname(fasta), "new_mutations.csv")
	else:
		mutations = mutations

	output = os.path.join(os.path.dirname(fasta), focus + "_outfile.csv")

	#generate coupling model
	c = CouplingsModel(model_params)
	if mode == "epistatic" or mode == "independent":
		data = pd.read_csv(mutations, sep=",", comment="#")
	
	#run EVmutation tools on coupling model, then convert data back to unaligned sequence
	if mode == 'epistatic':
		data_pred = ev_tools.predict_mutation_table(c, data, "effect_prediction_epistatic")
		data_pred.to_csv(output) #write data to output file
		align_muts_with_msa.align_with_orig(fasta, focus, output, outfile=output, pos_col=1) #realign positions with sequence (instead of MSA)
		#first column of output file is just an index, so ignore it (pos_col=1)
	elif mode == 'independent':
		c = c.to_independent_model()
		data_pred = ev_tools.predict_mutation_table(c, data, "effect_prediction_independent")
		data_pred.to_csv(output)
		align_muts_with_msa.align_with_orig(fasta, focus, output, outfile=output, pos_col=1) 
	elif mode == 'matrix':
		data_pred = ev_tools.single_mutant_matrix(c, output_column="effect_prediction_epistatic")
		data_pred = pd.DataFrame(dictize(data_pred, alphabet=alphabet[1:])).T
		data_pred.to_csv(output, sep=",")
		align_muts_with_msa.align_with_orig(fasta, focus, output, outfile=output, pos_col=0)
	data_pred = pd.read_csv(output, index_col=0)
	return(data_pred)
		
		#first column of output file is important (position), so keep it (pos_col=0)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='This script will run all the necessary subcomponents to compute EV scores, including PLMC, aligning mutations with the MSA, generating the coupling model, and outputting results. It is possible to skip the mutation alignment and PLMC if the model_params file from PLMC already exists and if the mutations are already aligned with the MSA. It is also possible to change the alphabet used, but remember to include the gap character in the alphabet. For best results, include as many arguments as possible.', add_help=True)
	parser.add_argument('--fasta', type=str, metavar='PATH', default="", help="path to multiple sequence alignment fasta file")
	parser.add_argument('--focus', type=str, default="", help="the focus sequence in the multiple sequence alignment")
	parser.add_argument('--model_params', type=str, metavar='PATH', default="", help="path to model_params file generated by plmc")
	parser.add_argument('--mutations', type=str, default="mutations.csv", metavar='PATH', help="path to CSV file containing mutants")
	parser.add_argument('--skip_plmc', action='store_true', help="whether or not to run PLMC to generate model parameters")
	parser.add_argument('--skip_align_muts', action='store_true', help="whether or not to align the mutations CSV file to the multiple sequence alignment")
	parser.add_argument('--mode', choices={'epistatic', 'independent', 'matrix'}, default="epistatic", help="which analysis to run (\'epistatic\' to compute the epistatic effects of the mutations in the spreadsheet, \'independent\' to compute the independent effects of the mutations in the spreadsheet, \'matrix\' to generate the matrix of values for all amino acid substitutions in all positions")
	parser.add_argument('--alphabet', type=str, help="the alphabet used in the MSA (-UAGC, -TAGC, -ACDEFGHIKLMNPQRSTVWY)", default=None)
	args = parser.parse_args()

	create_analysis(args.fasta, args.focus, args.model_params, args.mutations, args.skip_plmc, args.skip_align_muts, args.mode, args.alphabet)
