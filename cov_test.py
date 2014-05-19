#! /usr/bin/env python

## Ramakrishna Sompallae, Univ of Iowa
## pipeline to analyse and annotate iontorrent data
## included the var qual information
## last edited 7/23/13 5:55 PM
## filter vcf files
## saves analysed files to PGM_CHPv2_results 

import csv
import fnmatch
import operator
import os
import shutil
import sys
import time
import numpy as np

def inadequate_cov(cov_file, final_file):
	in_file = cov_file
	#data = np.genfromtxt('compare.dat',dtype=None,delimiter=18, autostrip=True,names=True, usecols=range(COLUMNS), comments=None)
	data = np.genfromtxt(in_file, dtype=None, autostrip=True)

	subset = [d for d in data if d[6] < 250]

	out_file = final_file
	output = csv.writer(open(out_file, 'a'))

	output.writerow(["Inadequately covered amplicons (<250X)"])
	output.writerow(["Gene", "Exon", "Amplicon", "Coverage"])

	for s in subset:
#		print s[5]+"\t"+s[4]+"\t"+s[3]+"\t"+str(s[6])
		output.writerow ([s[5], s[4], s[3], s[6]])

#inadequate_cov(NGS_results+sname+"/"+sname+"_coverage_regions.txt", NGS_results+sname+"/"+sname+"_final_to_report.csv")
inadequate_cov(sys.argv[1], sys.argv[2])
	