#! /usr/bin/env python
import os
import sys
import csv

in_file = sys.argv[1]                   # input file: first argument

out_file = '_'.join([sys.argv[2], "prefiltered.csv"])
out_file1 = '_'.join([sys.argv[2], "final_to_report.csv"])

print out_file1


# output in csv format
output = csv.writer(open(out_file, 'wb'))
output1 = csv.writer(open(out_file1, 'wb'))

count = 0	
for line in csv.reader(open(in_file, 'rU'), skipinitialspace=True, dialect="excel"):
# first line to out file
#	print len(line)
	count=count+1
	hotspot = line[14]
	line[14] = line[37]
	line[37] = hotspot
	
	if "Chrom" in line[0]:
		output.writerow (line[:43])
		output1.writerow ([line[0], line[1], line[2], line[18], line[19], line[17], line[14], line[4], line[6], line[7], line[8], line[9], line[12], line[13], line[24], line[23], line[22], line[3], line[37], line[27], line[29], line[31], line[33], line[35], line[36]])
		continue	
	
	line[6] = line[41] #Ref
	line[7] = line[42] #Var
	output.writerow (line[:43])
	
# filter out calls with freq < 5, Qual < 25, coverage < 100 and intronic and synonymous SNVs
	
	if float(line[8]) < 5 or int(float(line[9])) < 25 or int(float(line[11])) < 100 or line[15] == "intronic" or line[17] == "synonymous SNV":
		continue
	
#filter out MAF greater than 0.1
	else:
		maf = line[23].replace('\"','') #strip "" 
#		print maf
		if line[23] == '' or float(maf) <= 0.1:
			output1.writerow ([line[0], line[1], line[2], line[18], line[19], line[17], line[14], line[4], line[6], line[7], line[8], line[9], line[12], line[13], line[24], line[23], line[22], line[3], line[37], line[27], line[29], line[31], line[33], line[35], line[36]])
		elif float(maf) > 0.1 and line[14] != '':
			output1.writerow ([line[0], line[1], line[2], line[18], line[19], line[17], line[14], line[4], line[6], line[7], line[8], line[9], line[12], line[13], line[24], line[23], line[22], line[3], line[37], line[27], line[29], line[31], line[33], line[35], line[36]])

output1.writerow("")
output1.writerow("")
output1.writerow(["Number of variants found by TVC 3.4.5:",str(count - 1)])
output1.writerow("")
output1.writerow("")
				
'''
# for all lines of csv input, load them in, filter, and output.
for (query_name, subject_name, score, expect) in csv.reader(open(in_file)):
    # filter HERE.
    if float(expect) < 1e-3:            # e.g. keep only good matches
        row = [query_name, subject_name, score, expect]
        output.writerow(row)
'''
