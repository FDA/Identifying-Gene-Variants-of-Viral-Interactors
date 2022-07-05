import os, os.path, argparse, re, copy

#Input: directory containing all score files
#Output: A dictionary of collected scores for each file/structure
#Utility function to parse all Rosetta score files in any directory and organize the results into one file, sorted on the specified column
#If rm is set to True, it will remove all files except for the top five to save space
def parse_scores(directory, col="total_score", outfile="scores.tsv", prefix=None, rm=False):
	if outfile == "scores.tsv":
		outfile = os.path.join(directory, outfile)
	data = {}
	colnames = {}
	for f in os.listdir(directory):
		if f.endswith(".sc"):
			with open(os.path.join(directory, f), "r") as inf:
				lines = list(inf.readlines())
			lines = [lines[i].split() for i in range(len(lines))]
			try:
				if len(colnames) < 1:
					colnames = {k-1:lines[1][k] for k in range(1, len(lines[1]))}
				data[f.replace(".sc", "")] = {colnames[k]:lines[2][k+1] for k in colnames.keys()}
			except Exception as e:
				print(e)
	data = sorted(data.items(), key = lambda kv: float(kv[1][col]), reverse=False)
	with open(outfile, "w") as outf:
		outf.write("\t" + "\t".join(colnames.values()) + "\n")
		for k,v in data:
			outf.write(str(k) + "\t")
			outf.write("\t".join([v[colnames[k]] for k in colnames.keys()]) + "\n")
	#print([data[i][0] for i in range(5)])
	if rm:
		#print(data[0][0])
		if prefix == None or len(prefix) < 1:
			prefix=max(re.findall("(.+)\_\d+", data[0][0]))
		#print(prefix)
		top5 = [max(re.findall(prefix + "\_(\d+)", data[i][0])) for i in range(5)]
		#print(top5)
		for f in os.listdir(directory):
			if f.endswith(".sc"):
				try:
					index = max(re.findall(prefix + "\_(\d+)",os.path.splitext(os.path.basename(f))[0]))
					if index not in top5:
						os.remove(os.path.join(directory, f))
				except ValueError as e:
					pass
			elif f.endswith(".pdb.gz"):
				try:
					index = max(re.findall(prefix + prefix + "\_(\d+)\_", os.path.splitext(os.path.basename(f))[0]))
					if index not in top5:
						os.remove(os.path.join(directory, f))
				except ValueError as e:
					pass
	return(data)

#Can parse through a directory to organize and sort all Rosetta score files
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="This script will read through all score files, and find the lowest five scores.", add_help=True)
	parser.add_argument("directory", type=str, metavar="PATH", help="path to directory containing score files")
	parser.add_argument("--col", "-c", type=str, default="total_score", help="the number of the column to use (0-indexed)")
	parser.add_argument("--outfile", "-o", type=str, metavar="FILE", default="scores.tsv", help="file to output all data to")
	parser.add_argument("--prefix", type=str, default=None, help="prefix of files to be analyzed (useful when multiple Rosetta outputs stored in same directory)")
	parser.add_argument("--rm", action="store_true", default=False, help="delete files for all but top five")
	clargs = parser.parse_args()

	data = parse_scores(clargs.directory, col=clargs.col, outfile=clargs.outfile, prefix=clargs.prefix, rm=clargs.rm)
