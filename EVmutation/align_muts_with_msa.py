# David Holcomb
# 2018

import pandas as pd
import argparse, subprocess, re, os.path
import warnings

indel_chars = ['-', '.', '*']

#read the msa for the focus sequence
def read_msa(fasta, focus):
	if focus == None or len(focus) <= 0:
		raise Exception("Cannot align mutations against the MSA without the name of the focus sequence.")

	#read the focus sequence from the MSA fasta
	wt_seq = ""
	focusflag = False
	with open(fasta, 'r') as f:
		for line in f:
			if len(line) > 1 and '#' not in line: #read line
				if '>' in line: #sequence name
					focusflag = (line[1:].strip() == focus) #indicates whether we're reading the focus sequence
				else: #sequence data
					if focusflag:
						wt_seq = wt_seq + line.strip()
	f.close()
	return wt_seq

def get_mutation_data(line):
	pos = 0
	aas = re.findall("[a-zA-Z]", line)
	nums = re.findall("\d+", line)
	for num in nums:
		if int(num) > pos:
			pos = int(num)
	return pos, aas

def read_muts(mutations):
	#read the mutations
	muts = []
	with open(mutations, 'r') as f:
		for line in f:
			if len(line) > 1 and '#' not in line and "mutant" not in line:
				line = re.findall("[a-zA-Z]\d+[a-zA-Z]", line)[0]
				pos, aas = get_mutation_data(line)
				if len(aas) < 2 or pos == 0:
					raise Exception("Mutation line formatted incorrectly: " + line)
				else:
					muts.append((pos, aas[0], aas[1]))
	f.close()
	return muts

#converts a list of positions/mutations that have been aligned to a MSA back to the original format
#note that there may be duplicate amino acid mutations due to synonymous mutations
def align_with_orig(fasta, focus, csv, outfile="", pos_col=0, header=True, verbose=False):
	if verbose:
		print("Reverting mutation positions for ease of use.")

	#read the focus sequence from the MSA fasta
	wt_seq = read_msa(fasta, focus)

	#create a list of position/mutation and row data
	l = []
	with open(csv, "r") as f:
		for line in f:
			parts = line.split(",")
			l.append(parts[pos_col:])
	f.close()

	#convert position into new position by subtracting number of '-' characters in sequence
	new_l = []
	for i, line in enumerate(l):
		if i == 0 and header: #append the header to the new list of lines without changing
			new_l.append(line)
		else:
			pos, aas = get_mutation_data(line[0])
			if pos == 0:
				continue
			offset = 0
			for c in indel_chars:
				offset = offset + int(wt_seq[:pos].count(c))
			new_pos = int(pos) - offset
			if len(aas) < 2:
				new_key = str(new_pos)
			else:
				new_key = aas[0] + str(new_pos) + aas[1]
			new_l.append([new_key] + line[1:])
	
	if outfile == None or len(outfile) < 1:
		outfile = os.path.join(os.path.dirname(csv), "outfile.csv")

	#print new list to output file
	with open(outfile, "w") as f:
		for line in new_l:
			for cell in line:
				f.write(str(cell).strip() + ",")
			f.write('\n')
	f.close()

#aligns a list of mutations to a MSA alignment (by accounting for '-' characters)
def align_with_msa(fasta, focus, mutations, verbose=False):
	
	if verbose:
		print("Aligning mutations with multiple sequence alignment...")

	#read the focus sequence from the MSA fasta
	wt_seq = read_msa(fasta, focus)

	#read the mutations
	muts = read_muts(mutations)

	#relates original mutations to aligned mutations
	correspondence = {}

	#sort the mutations, and adjust their positions to account for insertions/deletions ('-')
	sorted_muts = sorted(muts, key=lambda tup: tup[0])
	offset = 0	#number of '-' characters in sequence
	i = 0		#position in sequence
	for mut in sorted_muts:
		while i < mut[0] + offset:
			if wt_seq[i] in indel_chars:
				offset = offset + 1
			i = i + 1
		if wt_seq[i-1] is not mut[1]:
			warnings.warn("Mutations don't line up with focus sequence: \n" + str(i) + wt_seq[i] + "\n" + str(mut[0]) + mut[1])
		else:
			new_mut = (str(i), mut[1], mut[2])
			correspondence[mut] = new_mut

	#write out newly aligned mutations (maintains original ordering)
	with open(os.path.join(os.path.dirname(fasta), "new_mutations.csv"), "w") as f:
		f.write("mutant\n")
		for mut in muts:
			if mut in correspondence:
				new_mut = correspondence[mut]
				f.write(new_mut[1] + str(new_mut[0]) + new_mut[2] + "\n")
		f.close()

	#write out a correspondence between the old mutations and the newly aligned mutations (also maintains original ordering)
	#probably useful later
	with open(os.path.join(os.path.dirname(fasta), "muts_correspondence.csv"), "w") as f:
		for mut in muts:
			if mut in correspondence:
				f.write(mut[1] + str(mut[0]) + mut[2] + ",")
				new_mut = correspondence[mut]
				f.write(new_mut[1] + str(new_mut[0]) + new_mut[2] + "\n")
	f.close()
	
	print("Done.")

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='This script will align a list of mutations with a multiple sequence alignment, where positions may have changed due to insertions and deletions. It can also revert aligned positions to the original alignment (align with the original sequence).', add_help=True)
	parser.add_argument('fasta', type=str, metavar='PATH', default="", help="path to multiple sequence alignment fasta file")
	parser.add_argument('focus', type=str, default="", help="the focus sequence in the multiple sequence alignment")
	parser.add_argument('mutations', type=str, default="AA_muts.csv", metavar='PATH', help="path to CSV file containing mutants")
	parser.add_argument('--outfile', type=str, metavar='PATH', default="", help="path to the output file")
	parser.add_argument('--pos_col', type=int, default=0, help="index of the column containing the position (mutation information")
	parser.add_argument('--mode', choices={'align', 'unalign'}, default="align", help="whether to align your mutations with the MSA, or to unalign them (align them with the original sequence")
	args = parser.parse_args()
	
	if args.mode == "align":
		align_with_msa(args.fasta, args.focus, args.mutations)
	elif args.mode =="unalign":
		align_with_orig(args.fasta, args.focus, args.mutations, outfile=args.outfile, pos_col=args.pos_col)
	else:
		print("Unrecognized mode: please specify \'align\' to adjust your mutation positions to align them with the MSA, or \'unalign\' to return your mutation positions to be aligned with only the focus sequence.")
	
