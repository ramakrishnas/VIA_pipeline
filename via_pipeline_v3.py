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

global annovar, amplicon_bed, coverage_dir, watchdir, hg19_index, CHPv1_bed, CHPv2_bed, CCP_bed, CP_NGa_bed, hotspot_bed, chiptype, runsdb, mut_test, temp_hs_bed

somatic_tests = ["CHPv1", "CHPv2", "CP_NGa", "CCP", "MLP", "SCP"]
germline_tests = ["SMCHD1", "DYSF"]
mut_test= "somatic"
run_dir = sys.argv[1]

sys.path.append('/Users/molpathuser1/')
#os.chdir('/Users/molpathuser1/')

#test_folder
#watchdir = "/Volumes/Macintosh HD/Users/molpathuser1/ngs_testfolder/"

watchdir = "/Volumes/Iontorrent/"



runsdb = "/Volumes/Iontorrent/VIA_pipeline/runsdb.txt"
tmp_dir = "/Volumes/Iontorrent/VIA_pipeline/tmp/"

###molpath1
annovar = "/Volumes/Iontorrent/annovar"

NGS_results = "/Volumes/Iontorrent/All_NGS_results/"

#NGS_results = "/Volumes/Iontorrent/NGS_results"

CHPv2_results = "/Volumes/Iontorrent/PGM_CHPv2_results/"

shared_drive = "/Volumes/pathology/Molecular\ Pathology/All_NGS_results/ "

CHPv1_bed = "/Volumes/Iontorrent/ampliseq_files/Cancer_hotspot_panel_1/HSMv12.1_reqions_NO_JAK2_NODUP.bed"
CHPv2_bed = "/Volumes/Iontorrent/ampliseq_files/Cancer_hotspot_panel_2/4475346_CHP2_designed_20120806.bed"
CCP_bed = "/Volumes/Iontorrent/ampliseq_files/Comprehensive_cancer_panel/CCP_Amplicons.bed"
CP_NGa_bed = "/Volumes/Iontorrent/ampliseq_files/Custom_panel/IAD34427_Designed_v2.bed"
MLP_bed = "/Volumes/Iontorrent/ampliseq_files/Melanoma_panel/IAD40174_30_Designed_genes.bed"
#SCP_bed = "/Volumes/Iontorrent/ampliseq_files/Small_cancer_panel/New_Small_Panel_St.bed"
SCP_bed = "/Volumes/Iontorrent/ampliseq_files/Small_cancer_panel/SCP_IAD54614_124_Designed.bed"
SMCHD1_bed = "/Volumes/Iontorrent/ampliseq_files/SMCHD1/SMCHD1.bed"
DYSF_bed = "/Volumes/Iontorrent/ampliseq_files/SMCHD1/DYSF.bed"

CHPv1_hs_bed = "/Volumes/Iontorrent/ampliseq_files/Cancer_hotspot_panel_1/HSMv12.1_hotspots_NO_JAK2.bed"
CHPv2_hs_bed = "/Volumes/Iontorrent/ampliseq_files/Cancer_hotspot_panel_2/4475346_CHP2_hotspots_20120927.bed"
CCP_hs_bed = "/Volumes/Iontorrent/ampliseq_files/Comprehensive_cancer_panel/CCP_hotspots_20121225.bed"
CP_NGa_hs_bed = "/Volumes/Iontorrent/ampliseq_files/Custom_panel/IAD34427_Submitted_v2.bed"
SCP_hs_bed = "/Volumes/Iontorrent/ampliseq_files/Cancer_hotspot_panel_2/4475346_CHP2_hotspots_20120927.bed"
temp_hs_bed = "/Volumes/Iontorrent/ampliseq_files/hs_temp.bed"

hg19_index = "/Volumes/Iontorrent/ampliseq_files/hg19.fasta.fai"

coverage_dir = "/Volumes/Iontorrent/amplicon_analysis/coveragebed"

####bioinfo core
#annovar = "/Users/ramakrishnasompallae/Work/Pathology/Mol_Path/Iontorrent/annovar"

#amplicon_bed = "/Users/ramakrishnasompallae/Work/Pathology/Mol_Path/Iontorrent/ampliseq_files/HSMv12.1_reqions_NO_JAK2_NODUP.bed"

#coverage_dir = "/Users/ramakrishnasompallae/Work/Pathology/Mol_Path/Iontorrent/amplicon_analysis/coveragebed"


##looks for vcf file under each barcode

'''
def get_folder_size(folder_path):
	print folder_path
	total_size = 0
	for filenames in os.listdir(folder_path):
		if filenames.endswith("PTRIM.bam") and:
		total_size += os.path.getsize(folder_path+"/"+filenames)
		print total_size
	return total_size
'''    
    
def vcf_rename(folder, num):
	barcodes = os.listdir(folder)
	for i in barcodes:
		if "IonXpress" in i:
			for filenames in os.listdir(folder+"/"+i):
				if filenames.endswith("PTRIM.bam") and os.path.getsize(folder+"/"+i+"/"+filenames) > 10000000:
					runname = "R_"+num
					samplename = "R"+num+"_BC"+i.split("_")[1]
					tempdir = folder+"/"+i	
					vcf_annovar(tempdir, samplename, runname)		

##coverage analysis of amplicon regions
def amp_coverage(path, folder, num):
	global amplicon_bed, hotspot_bed, samplename, NGS_results, shared_drive
	for files in os.listdir(os.path.join(path+folder)):
		tempfile = path+folder+"/"+files
		print tempfile
		if tempfile.endswith(".bam") and os.path.getsize(tempfile) > 10000000: # if bam file greater than 10mb
			samplename = "R"+num+"_BC"+(files.split("_")[1])
#			print files
#			print filesize
					
			# run 21 was not labelled with CHPv2			
			if "CHPv2" in folder or "CHP2" in folder or "MOL-21-316A" in folder or "MOL-60" in folder or "Mol-72_lin_097" in folder:
				chiptype = "CHPv2"
				amplicon_bed = CHPv2_bed
				hotspot_bed = CHPv2_hs_bed
				NGS_results = CHPv2_results
				shared_drive = "/Volumes/pathology/Molecular\ Pathology/Active\ Validated\ Tests/NGS\ CHPv2\ data/PGM\ CHPv2\ results/"
						
			elif "CP_NGa" in folder:
				chiptype = "CP_NGa"
				amplicon_bed = CP_NGa_bed
				hotspot_bed = CP_NGa_hs_bed

			elif "CCP" in folder:
				chiptype = "CCP"
				amplicon_bed = CCP_bed
				hotspot_bed = CCP_hs_bed

			elif "MLP" in folder:
				chiptype = "MLP"
				amplicon_bed = MLP_bed
				hotspot_bed = CCP_hs_bed
				
			elif "SCP" in folder:
				chiptype = "SCP"
				amplicon_bed = SCP_bed
				hotspot_bed = SCP_hs_bed	

			elif "SMCHD1" in folder:
				chiptype = "SMCHD1"
				amplicon_bed = SMCHD1_bed
				hotspot_bed = temp_hs_bed
				mut_test = "germline"

			elif "DYSF" in folder:
				chiptype = "DYSF"
				amplicon_bed = DYSF_bed
				hotspot_bed = temp_hs_bed
				mut_test = "germline"

			else:
				chiptype = "CHPv1", "CHPv2", "CP_NGa", "CCP", "MLP", "SCP"
				amplicon_bed = CHPv1_bed
				hotspot_bed = CHPv1_hs_bed
							
			print amplicon_bed
			os.system("mkdir "+annovar+"/R_"+num)
			os.system("mkdir "+NGS_results+samplename)
			os.system("cp "+tempfile+" "+NGS_results+samplename)
			os.system("cp "+tempfile+".bai "+NGS_results+samplename)
			os.system("/Users/molpathuser1/othertools/BEDTools-Version-2.16.2/bin/coverageBed -abam "+ tempfile +" -b "+ amplicon_bed +" > "+NGS_results+samplename+"/"+samplename+"_coverage_regions.txt")
'''			if os.path.exists(path+folder+"/plugin_out/createPdf_out/"):
				os.system("cp "+path+folder+"/plugin_out/createPdf_out/report.pdf "+NGS_results+"/"+samplename+"/"+samplename+"_report.pdf")
			else:
				os.system("cp "+path+folder+"/report.pdf "+NGS_results+"/"+samplename+"/R"+num+"_"+chiptype+"_report.pdf")
'''

def create_avinput_vcf(tempdir, varfile):
	ids = ['Zygosity', 'Het', 'Hom']
	out = open(tempdir+"/var.avinput", 'w')
	with open(varfile) as f:
		for line in f:
			parts = line.strip().split('\t')
			if parts[5] in ids:
				out.write(line),
	f.close()
	out.close()


def sort_csv(hsfile):
	in_file = hsfile
	out_file = in_file.replace(".csv", "")+"_temp.csv"
	out_file2 = in_file.replace(".csv", "")+"_sort.csv"

	with open(in_file, "rU") as f1:
		lines = csv.reader(f1)
		header = lines.next()
		thedata = list(lines)

	if "ann_out.genome_summary" in in_file:
		thedata.sort(key=operator.itemgetter(28,29))

	else:
		thedata.sort(key=operator.itemgetter(0,1))
	
	with open(out_file, 'w') as f2:
		writeit = csv.writer(f2)
		f2.write(str(header).strip("[]")+"\n")
		writeit.writerows(thedata)

	with open(out_file, 'rU') as infile, open(out_file2, 'wb') as outfile:
		for line in infile:
			outfile.write(line.replace('\r', ''))
    	os.remove(out_file)	

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

##Annotation for variants using Annovar			
def vcf_annovar(tempdir, sname, rname):
	ann_folder = annovar+"/"+rname
	#var_file = tempdir+"/variants.xls"
	vcf_file = tempdir+"/TSVC_variants.vcf"
	#create_avinput_vcf(tempdir, var_file)
	print amplicon_bed
	print os.path.abspath('')
	print os.system("echo $PATH > test.txt") #test
	
	os.system("/Volumes/Iontorrent/VIA_pipeline/vcf_parser_2.py "+vcf_file+" "+tmp_dir+"temp.xls "+amplicon_bed) # vcf parser to get read depth
	os.system("cut -f1,2,5,7-14 "+tmp_dir+"temp.xls > "+tmp_dir+"temp.vcf") # rearrange columns for annovar
	os.system("/Users/molpathuser1/othertools/annovar/convert2annovar.pl --format vcf4 --includeinfo "+tmp_dir+"temp.vcf > "+annovar+"/"+rname+"/"+sname+".avinput")

	if (mut_test == "somatic" and chiptype in somatic_tests):
		os.system("/Users/molpathuser1/othertools/annovar/summarize_annovar_cosmic.pl -buildver hg19 -ver1000g 1000g2012apr -verdbsnp 137 -veresp 6500 -remove "+annovar+"/"+rname+"/"+sname+".avinput /Users/molpathuser1/othertools/annovar/humandb -outfile "+annovar+"/"+rname+"/"+sname+"_ann_out")
	
	elif (mut_test == "germline" and chiptype in germline_tests):
		os.system("/Users/molpathuser1/othertools/annovar/summarize_annovar_cosmic.pl -buildver hg19 -ver1000g 1000g2012apr -verdbsnp 137 -veresp 6500 -remove "+annovar+"/"+rname+"/"+sname+".avinput /Users/molpathuser1/othertools/annovar/humandb -outfile "+annovar+"/"+rname+"/"+sname+"_ann_out")

#	os.system("table_annovar.pl "+annovar+"/"+rname+"/"+sname+".avinput humandb/ -buildver hg19 -out anno -protocol refGene,phastConsElements46way,genomicSuperDups,1000g2012apr_all,esp6500si_all,snp137,ljb_all,cosmic64 -operation g,r,r,f,f,f,f,f -nastring NA -csvout")

	ann_out_file = annovar+"/"+rname+"/"+sname+"_ann_out.genome_summary.csv"
	os.system("/Volumes/Iontorrent/VIA_pipeline/collate_variants.pl -l -G "+hg19_index+" -B "+ hotspot_bed+ " "+tmp_dir+"temp2.xls "+tmp_dir+"temp.xls")
	stringl = "awk -F'\t' '{ for (n=1; n <= NF; ++n) { if (n>1) printf(\",\"); printf(\"\\\"%s\\\"\",$n); } print \"\"; }' "+tmp_dir+"temp2.xls > "+annovar+"/"+rname+"/"+sname+"_hotspot.csv"
	os.system(stringl)

	num_vcf = sum(1 for line in open(ann_out_file)) ## from annotated csv file
	print "Number of variants in VCF file : "+ str(num_vcf - 1)
	sort_csv(ann_out_file)

	hs_file = annovar+"/"+rname+"/"+sname+"_hotspot.csv"
	sort_csv(hs_file)
	
	stringp = "paste -d \", \" "+ann_folder+"/"+sname+"_hotspot_sort.csv "+ann_folder+"/"+sname+"_ann_out.genome_summary_sort.csv > "+ann_folder+"/"+sname+"_tmp.csv"
	os.system(stringp)
	os.system("/Volumes/Iontorrent/VIA_pipeline/filter_annovar_v4.py "+ann_folder+"/"+sname+"_tmp.csv "+ann_folder+"/"+sname)
	
	os.system("cp "+ann_folder+"/"+sname+"_prefiltered.csv "+NGS_results+sname+"/")
	os.system("cp "+ann_folder+"/"+sname+"_final_to_report.csv "+NGS_results+sname+"/")
	os.system("cp "+vcf_file+" "+NGS_results+sname+"/")
	
	inadequate_cov(NGS_results+sname+"/"+sname+"_coverage_regions.txt", NGS_results+sname+"/"+sname+"_final_to_report.csv")
	print ann_folder+"/"+sname+"_final_to_report.csv"
	
	os.system("rm "+tmp_dir+"temp.*")

def analyze_data(dir):
	run = dir.split("_")[2][4:]
	print run
	runs.append(run)
	amp_coverage(watchdir, dir, run)
	if os.path.exists(watchdir+dir+"/plugin_out/variantCaller3.4.2_out"):
		tempdir = os.path.join(watchdir+dir+"/plugin_out/variantCaller3.4.2_out")
		samplenames = list()
		barcodes = os.listdir(tempdir)
		for i in barcodes:
			if ("IonXpress" in i) and os.path.exists(watchdir+dir+"/plugin_out/variantCaller3.4.2_out/"+i+"/TSVC_variants.vcf"):
				print i+" variants present"
				samplenames.append(i)
			else:
				sys.exit("Check the folder "+watchdir+dir+"/plugin_out/variantCaller3.4.2_out")
		vcf_rename(tempdir, run)

dirs = list()
runs = list()

localtime = time.asctime(time.localtime(time.time())) 
if os.path.exists(runsdb):
	oldcontents = open(runsdb, "r").read().splitlines()

print run_dir
if "_MOL-" in run_dir:
	print run_dir
	if run_dir in oldcontents:
		print run_dir+" is previously analysed"
	else:
		print "Now analyzing "+run_dir
		analyze_data(run_dir)
		f = open(runsdb, 'a')
		f.write(localtime+'\n')
		f.write(run_dir +'\n')
else:
	print "Please check the directory name."


#myruns = sorted(set(runs))
#print myruns
print localtime

#os.system("rsync -avz /Volumes/Iontorrent/NGS_results/ /Volumes/pathology/Molecular\ Pathology/Test\ Development/IonTorrent\ PGM/VALIDATION\ DATA/NGS_results/")

#os.system("rsync -avz /Volumes/Iontorrent/PGM_CHPv2_results/ /Volumes/pathology/Molecular\ Pathology/Active\ Validated\ Tests/NGS\ CHPv2\ data/PGM\ CHPv2\ results/")
print "copying files from "+NGS_results+" to "+shared_drive
os.system("rsync -avz "+NGS_results+" "+shared_drive)

