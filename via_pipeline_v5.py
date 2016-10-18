#! /usr/bin/env python

## Ramakrishna Sompallae, Univ of Iowa
## pipeline to analyse and annotate iontorrent data
## included the var qual information
## edited to analyse the results from 2nd PGM and SCP 20 gene panel
## last edited 4/1/15 5:53 PM
## filter vcf files
## saves analysed files to PGM_CHPv2_results 

import csv
import fnmatch
import operator
import glob
import os
import shutil
import sys
import time
import re
import numpy as np

global annovar, amplicon_bed, coverage_dir, watchdir, run_dir, var_caller, hg19_index, CHPv1_bed, CHPv2_bed, \
CCP_bed, CP_NGa_bed, hotspot_bed, chiptype, runsdb, sequencer, mut_test, fil_par, \
temp_hs_bed, var_caller, tvc, somatic_tests, germline_tests, bc_chip_match

somatic_tests = ["CHPv1", "CHPv2", "CP_NGa", "CCP", "MLP", "SCP", "AML", "MPCP", "TP53", "Myeloma", "Melanoma"]
germline_tests = ["SMCHD1", "DGP", "DYSF"]
mut_test = "somatic"

run_dir = sys.argv[1]

#bc_chip_match = sys.argv[2]

sys.path.append('/Users/molpathuser1/')

runsdb = "/Volumes/Iontorrent/VIA_pipeline/runsdb.txt"
tmp_dir = "/Volumes/Iontorrent/VIA_pipeline/tmp/"

###molpath1
'''
NGS_results = "/Volumes/Iontorrent/All_NGS_results/"

CHPv2_results = "/Volumes/Iontorrent/PGM_CHPv2_results/"

SCP_results = "/Volumes/Iontorrent/All_NGS_results/SCP_results/"
MPCP_results = "/Volumes/Iontorrent/All_NGS_results/MPCP_results/"

shared_drive = "/Volumes/pathology/Molecular\ Pathology/All_NGS_results/"
'''

CHPv1_bed = "/Volumes/Iontorrent/ampliseq_files/Cancer_hotspot_panel_1/HSMv12.1_reqions_NO_JAK2_NODUP.bed"
CHPv2_bed = "/Volumes/Iontorrent/ampliseq_files/Cancer_hotspot_panel_2/CHP2/CHP2.20131001.designed.bed"
#"/Volumes/Iontorrent/ampliseq_files/Cancer_hotspot_panel_2/4475346_CHP2_designed_20120806.bed"#old file
CCP_bed = "/Volumes/Iontorrent/ampliseq_files/Comprehensive_cancer_panel/CCP_Amplicons.bed"
CP_NGa_bed = "/Volumes/Iontorrent/ampliseq_files/Custom_panel/IAD34427_Designed_v2.bed"
#MLP_bed = "/Volumes/Iontorrent/ampliseq_files/Melanoma_panel/MLP_IAD40174_30_Designed_genes.bed" Previous/old bed file
SCP_bed = "/Volumes/Iontorrent/ampliseq_files/Small_cancer_panel/SCP_IAD54614_124_Designed.bed"
SCP20_bed="/Volumes/Iontorrent/ampliseq_files/Small_cancer_panel/SCP20_IAD74300_182_modified.bed"
#MPCP_bed = "/Volumes/Iontorrent/ampliseq_files/molpath_panel/MPCP23_IAD77734_236_modified.bed"
MPCP_bed = "/Volumes/Iontorrent/ampliseq_files/molpath_panel/MPCP.bed"

MPCP_ARMS_bed = "/Volumes/Iontorrent/ampliseq_files/molpath_panel/MPCP_ARMS.bed"

AML_bed = "/Volumes/Iontorrent/ampliseq_files/AML_panel/AML_IAD63562_182_Designed_npm1.bed"

SMCHD1_bed = "/Volumes/Iontorrent/ampliseq_files/SMCHD1_DYSF/SMCHD1.bed"
DYSF_bed = "/Volumes/Iontorrent/ampliseq_files/SMCHD1_DYSF/DYSF_IAD49541_124_Designed_mod.bed"
DGP_bed = "/Volumes/Iontorrent/ampliseq_files/Dystroglycanopathy_panel/DGP_IAD73267_204_Designed.bed"
TP53_bed = "/Volumes/Iontorrent/ampliseq_files/TP53/TP53.20140108.designed.bed"
myeloma_bed = "/Volumes/Iontorrent/ampliseq_files/Myeloma_panel/Myeloma_IAD87360_238_modified.bed"
melanoma_bed = "/Volumes/Iontorrent/ampliseq_files/Melanoma_panel/IAD77906_231_MLP.bed"


hg19_index = "/Volumes/Iontorrent/ampliseq_files/hg19.fasta.fai"

coverage_dir = "/Volumes/Iontorrent/amplicon_analysis/coveragebed"

NGS_results = "/Volumes/Iontorrent/NGS_V5_Results/" #Outdir_local
shared_drive = "/Volumes/pathology/Molecular\ Pathology/Active\ Validated\ Tests/NGS\ V5\ data/S5\ Results/" #Outdir on S-drive


#shared_drive = "/opt/molpath_sdrive/Active\ Validated\ Tests/NGS\ V5\ data/Test\ S5\ Results/" #Test Outdir on S-drive

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

print len(sys.argv)

def find_var_Dir(directory, var_dir):
#	print directory
#	print var_dir
	os.chdir(directory)
	dirs = {}
	for dir in glob.glob(var_caller+'*'):
		print dir
#		if os.path.isdir(dir) and var_dir in dir:
		dirs[dir] = os.path.getctime(dir)

	all_vardirs = sorted(dirs.iteritems(), key=operator.itemgetter(1))
#	print lister
	return all_vardirs[-1][0] #get last item [-1] and folder name [0]
    

def somatic_filter(in_file):
	with open(in_file,"rb") as inf:
		out_file = ''.join([in_file.strip("prefiltered.tsv"), "final_to_report.csv"]) #final to report
		out_file2 = ''.join([in_file.strip("prefiltered.tsv"), "filtered.csv"]) #filtered to report
		outf = csv.writer(open(out_file, 'wb'))#final to report
		outf2 = csv.writer(open(out_file2, 'wb'))  #filtered

		for line in inf:
			parts = line.strip().split('\t',)
			
			hgvsc = ''
			hgvsp = ''
			exon = ''

			if parts[0] == 'CHROM':
				outf2.writerow(parts)
#	   			'\t'.join([str(x) for x in list])
#	   			outf.write(parts[0]+'\t'+parts[3]+'\t'+parts[16]+'\tRefseqID\t'+parts[21]+'\tHGVSvariant\t'+parts[4:6]+'\t'+parts[8:10]+'\t'+parts[14]+'\n')
				outf.writerow([parts[0], parts[1], parts[16], parts[21], 'HGVSvariant', parts[4], 'RefseqID', parts[5], parts[6], parts[8], parts[9], parts[10], parts[14]])
	   		
	   		elif float(parts[4]) >= 0.04 and int(parts[8]) >= 100 and parts[14] != 'synonymous_variant' and parts[14] != 'intron_variant':
		   		parts[4] = float(parts[4])*100
		   		
				if len(parts[21].split('/')) > 1:
					exon = parts[21].split('/')[0]
					parts[21] = exon+'('+parts[21].split('/')[1]+')'

				else:
					exon = parts[21]
		   		
	  			if len(parts[23].split(':')) > 1:
	  				hgvsc = parts[23].split(':')[1]
	  			else:
	  				hgvsc = parts[23]
	  				
	  			if len(parts[24].split(':')) > 1:
	  				hgvsp = parts[24].split(':')[1]
	  			else:
	  				hgvsp = parts[24]
	  			
	  			if hgvsc != '' and hgvsp != '':
	  				hgvs = hgvsc+':'+hgvsp
	  			elif hgvsc != '' and hgvsp == '':
	  				hgvs = hgvsc
	  			else:
	  				hgvs = ''

#	  			print hgvs	
	  			alt_freq = parts[4]

	 			if len(parts) < 52:
		   			outf2.writerow(parts)
#		   			outf.write(', '.join([parts[0], parts[1], parts[16], exon, hgvs,  "%.2f" %alt_freq, parts[19], parts[5].replace(",", "/"), parts[6], parts[8], parts[9], parts[10], parts[14], '\n']))
					outf.writerow([parts[0], parts[1], parts[16], exon, hgvs, "%.2f" %alt_freq, parts[19], parts[5], parts[6], parts[8], parts[9], parts[10], parts[14]])

				elif len(parts) > 52 and parts[52] == '':
		   			outf2.writerow(parts)
					outf.writerow([parts[0], parts[1], parts[16], exon, hgvs, "%.2f" %alt_freq, parts[19], parts[5], parts[6], parts[8], parts[9], parts[10], parts[14]])

				elif len(parts) > 52 and parts[52] != '' and float(parts[52].split(':',1)[1]) <= 0.1:
		   			outf2.writerow(parts)
					outf.writerow([parts[0], parts[1], parts[16], exon, hgvs, "%.2f" %alt_freq, parts[19], parts[5], parts[6], parts[8], parts[9], parts[10], parts[14]])

		


def coverage_filter(in_file):
	with open(in_file,"rb") as inf, open(''.join([in_file.strip("_cov_regs.txt"), "_low_cov_regs.csv"]),"wb") as outf:
	#	out_csv = csv.writer(outf, delimiter=',', quoting=csv.QUOTE_ALL)
		header = ['Chr', 'start', 'stop', 'amplicon', 'exon', 'Genes', 'Depth', 'targetbp', 'target_covered', 'target_portion']
		outf.write(', '.join([i for i in header]))
		ll = 0
		for line in inf:
			parts = line.strip().split('\t',)
			if	int(parts[6]) <= 250:
				ll += 1
				print parts[6]
				print "gret"
				outf.write('\n')
				outf.write(', '.join([i for i in parts]))
		return ll
			    
##coverage analysis of amplicon regions

def analyse_sample(barcode, folder, chip, sname): 
	global amplicon_bed, hotspot_bed, samplename, NGS_results, shared_drive, sequencer, chiptype, mut_test
	vardir = folder
	print barcode
	print folder
	print chip
	print len(chip)
	print sname	
	
	vcf_file = os.path.join(vardir+"/"+barcode+"/TSVC_variants.vcf.gz")	
	hs_calls = os.path.join(vardir+"/"+barcode+"/alleles.xls")	

	if chip.lower() == 'chpv2':
		chiptype = "CHPv2"
		amplicon_bed = CHPv2_bed
#		NGS_results = CHPv2_results
#		shared_drive = "/Volumes/pathology/Molecular\ Pathology/Active\ Validated\ Tests/NGS\ CHPv2\ data/PGM\ CHPv2\ results/"

	elif chip.lower() == 'scp':
		chiptype = "SCP"
		amplicon_bed = SCP_bed
#		NGS_results = SCP_results
#		shared_drive = "/Volumes/pathology/Molecular\ Pathology/Active\ Validated\ Tests/NGS\ SCP\ data/PGM\ SCP\ results/"
 
			#Kitchen sink / molpath custom gene panel
	elif chip.lower() == 'mpcp':
		chiptype = "MPCP"
		amplicon_bed = MPCP_bed
#		NGS_results = "/Volumes/Iontorrent/PGM_MPCP_results/"
#		shared_drive = "/Volumes/pathology/Molecular\ Pathology/Active\ Validated\ Tests/NGS\ MPCP\ data/PGM\ MPCP\ results/"
	
	elif chip.lower() == 'mpcp_arms' or chip.lower() == 'mpcparms':
		chiptype = "MPCP_ARMS"
		amplicon_bed = MPCP_ARMS_bed
#		NGS_results = "/Volumes/Iontorrent/PGM_MPCP_results/"
#		shared_drive = "/Volumes/pathology/Molecular\ Pathology/Active\ Validated\ Tests/NGS\ MPCP\ data/PGM\ MPCP\ results/"
		
	elif chip.lower() == 'aml':
		chiptype = "AML"
		amplicon_bed = AML_bed
#		NGS_results = "/Volumes/Iontorrent/PGM_AML_results/"
#		shared_drive = "/Volumes/pathology/Molecular\ Pathology/Active\ Validated\ Tests/NGS\ AML\ data/PGM\ AML\ results/"

	elif chip.lower() == 'smchd1':
		chiptype = "SMCHD1"
		amplicon_bed = SMCHD1_bed
#		NGS_results = "/Volumes/Iontorrent/PGM_SMCHD1_results/"
#		shared_drive = "/Volumes/pathology/Molecular\ Pathology/Active\ Validated\ Tests/NGS\ SMCHD1\ data/PGM\ SMCHD1\ results/"
		mut_test = "germline"

	elif chip.lower() == 'dysf':
		chiptype = "DYSF"
		amplicon_bed = DYSF_bed
#		NGS_results = "/Volumes/Iontorrent/PGM_DYSF_results/"
#		shared_drive = "/Volumes/pathology/Molecular\ Pathology/Active\ Validated\ Tests/NGS\ DYSF\ data/PGM\ DYSF\ results/"
		mut_test = "germline"

	elif chip.lower() == 'dgp':
		chiptype = "DGP"
		amplicon_bed = DGP_bed
#		NGS_results = "/Volumes/Iontorrent/PGM_DGP_validation/"
#		shared_drive = "/Volumes/pathology/Molecular\ Pathology/All_NGS_results/DGP_validation/"
		mut_test = "germline"

	elif chip.lower() == 'tp53':
		chiptype = "TP53"
		amplicon_bed = TP53_bed
#		NGS_results = "/Volumes/Iontorrent/PGM_TP53_validation/"
#		shared_drive = "/Volumes/pathology/Molecular\ Pathology/All_NGS_results/TP53_validation/"

	elif chip.lower() == 'myeloma':
		chiptype = "Myeloma"
		amplicon_bed = myeloma_bed
#		NGS_results = "/Volumes/Iontorrent/PGM_myeloma_validation/"		
#		shared_drive = "/Volumes/pathology/Molecular\ Pathology/All_NGS_results/myeloma_validation/"

	elif chip.lower() == 'melanoma':
		chiptype = "Melanoma"
		amplicon_bed = melanoma_bed

	print chip
# test print
	print amplicon_bed

	print NGS_results
	if not os.path.exists(NGS_results):
		os.makedirs(NGS_results)


	tvcf = NGS_results+sname+'/TSVC_variants.vcf.gz'
	
	os.system("mkdir "+NGS_results+sname)

	
	if os.path.exists(vardir.split("plugin_out")[0]+barcode+"_rawlib.bam"):
		bamfile = vardir.split("plugin_out")[0]+barcode+"_rawlib.bam"
		
	elif os.path.exists(vardir+"/"+barcode+"_rawlib.realigned.bam"):
		bamfile = vardir+"/"+barcode+"_rawlib.realigned.bam"
	
	elif os.path.exists(vardir+"/"+barcode+"/"+barcode+"_rawlib_processed.bam"):
		bamfile = vardir+"/"+barcode+"/"+barcode+"_rawlib_processed.bam"
		
	else:
		print 'Please check the folder for BAM files'
		
	os.system("cp "+bamfile+"* "+NGS_results+sname)
#	tvcf = NGS_results+sname+'/'+sname+vcf_file.split('/TSVC')[1]

	os.system("cp "+vcf_file+"* "+NGS_results+sname)
	
	if os.path.exists(vardir+"/"+barcode+"_parameters.json"):
		os.system("cp "+vardir+"/"+barcode+"_parameters.json "+NGS_results+sname)
	else:
		os.system("cp "+vardir+"/local_parameters.json "+NGS_results+sname)

	os.system("cp "+hs_calls+"* "+NGS_results+sname)
	
	os.system("/Users/molpathuser1/othertools/BEDTools-Version-2.16.2/bin/coverageBed -abam "+ bamfile +" -b "+ amplicon_bed +" > "+NGS_results+sname+"/"+sname+"_cov_regs.txt")

	os.system("run_vep.sh "+tvcf)
	os.system("mv "+tvcf.strip(".vcf.gz")+"_prefiltered.tsv "+NGS_results+sname+"/"+sname+"_prefiltered.tsv")

	lowl = coverage_filter(os.path.join(NGS_results+sname+"/"+sname+"_cov_regs.txt"))
	
	print lowl
	
	somatic_filter(os.path.join(NGS_results+sname+"/"+sname+"_prefiltered.tsv"))

#	stringp = "paste -d \"\\t\" "+NGS_results+"/"+sname+"/"+sname+"_final_to_report.tsv "+NGS_results+"/"+sname+"/"+sname+"_filtered.tsv | cut -f1-13,45,63-92 >"+NGS_results+"/"+sname+"/"+sname+"_temp_report.tsv"
#	os.system(stringp)
#	os.system("mv "+os.path.join(NGS_results+sname+"/"+sname+"_temp_report.tsv")+" "+os.path.join(NGS_results+sname+"/"+sname+"_final_to_report.tsv"))
#	print stringp
	
	with open(os.path.join(NGS_results+sname+"/"+sname+"_final_to_report.csv"),"a") as repf:
		repf.write("\n\n\n\nPanel given for this barcode "+chip+'\n')
		repf.write("Bed file used in the analysis "+amplicon_bed+'\n')
		repf.write("Total number of variants found _____\n")
		repf.write("There are %d low Covered amplicons with depth <250:\n"%lowl)
	os.system("cat "+os.path.join(NGS_results+sname+"/"+sname+"_low_cov_regs.csv")+" >> "+NGS_results+sname+"/"+sname+"_final_to_report.csv")

	print "copying files from "+NGS_results+" to "+shared_drive
	print NGS_results
	print shared_drive
	os.system("rsync -avz "+NGS_results+" "+shared_drive)

'''
	if "AML" in folder:
		for files in os.listdir(os.path.join(path+folder+"/plugin_out/downloads/")):
			tempfile2 = path+folder+"/plugin_out/downloads/"+files
			print tempfile
			if "IonXpress" in tempfile and tempfile.endswith(".fastq") and os.path.getsize(tempfile) > 5000000: # if bam file greater than 10mb
				os.system("usearch8.osx -search_oligodb "+tempfile+" -db /Volumes/Pathology/Molecular\ Pathology/Active\ Validated\ Tests/Oncology/AML\ Panel/FLT3-ITD/flt3_ex14_oligo.fa -strand both -matched "+NGS_results+samplename+"/FLT3-ITD/"+samplename+"_flt3_exon14.fa")
'''


def analyze_data(dir):
	global watchdir, var_caller, vardir, molid

	run = dir.split("_")[2][4:]
	print run
	runname = "R"+run
			
	print runname
	molid = ""	
	runs.append(run)
	watchdir = "/Volumes/Iontorrent/"
	sequencer = ""	

	if not os.path.exists(os.path.join(watchdir+dir)) and "SN2" in dir:		
		watchdir = "/Volumes/Iontorrent-2/"
		sequencer = "-SN2"

	elif not os.path.exists(os.path.join(watchdir+dir)) and "_S5XL-" in dir:
		watchdir = "/Volumes/Iontorrent-2/"
		run2 = run.split("-")[2]
		run=run2
		if "_S5XL-00642" in dir:
			sequencer = "_S5XLA"

		elif "_S5XL-00700" in dir:
			sequencer = "_S5XLB"
		
	bc = []
	chip = []
	try:
		with open(os.path.join(watchdir+dir+"/barcode_test.csv"), 'rb+') as bc_file:
			lines = (line for line in bc_file if 'IonXpress' in line)
			print lines
			for line in lines:
				print line
				bc.append(line.strip().split(",")[0])
				chip.append(line.strip().split(",")[1])

	except IOError as e:
		print "Barcode file, barcode_test.csv not found" #Does not exist OR no read permissions
		sys.exit(1)

	print "Number of samples in the run :"
	print len(bc)
	print bc	

	var_caller = "variantCaller_out"
	pluginsdir = os.path.join(watchdir, dir, "plugin_out")
	vardir = os.path.join(pluginsdir+"/"+find_var_Dir(pluginsdir, var_caller))
	print "vardir is "+vardir

	count = 0
	for i in bc:
		print i
		
		print vardir+"/"+i+"/TSVC_variants.vcf"
		
		print os.path.exists(vardir+"/"+i+"/TSVC_variants.vcf")
		
		if os.path.exists(watchdir+"/"+run_dir+"/"+i+"_rawlib.bam"):
			raw_bam = watchdir+"/"+run_dir+"/"+i+"_rawlib.bam"
			bam_size = os.path.getsize(watchdir+"/"+run_dir+"/"+i+"_rawlib.bam")
		
		elif os.path.exists(vardir+"/"+i+"_rawlib.realigned.bam"):
			raw_bam = vardir+"/"+i+"_rawlib.realigned.bam"
			bam_size = os.path.getsize(vardir+"/"+i+"_rawlib.realigned.bam")
		
		elif os.path.exists(vardir+"/"+i+"/"+i+"_rawlib_processed.bam"):
			raw_bam = vardir+"/"+i+"/"+i+"_rawlib_processed.bam"
			bam_size = os.path.getsize(vardir+"/"+i+"/"+i+"_rawlib_processed.bam")
	
		else:
			print "BAM file not found, Please check and rerun the analysis"
			
		if os.path.exists(vardir+"/"+i+"/TSVC_variants.vcf") and bam_size > 10000000:
			print "bam and variant files present for "+i
			print "panel is "+chip[count]
			with open(os.path.join(vardir+"/"+i+"/TSVC_variants.vcf"),"rb") as vcff:
				for line in vcff:
					if '#CHROM' in line and 'POS' in line and 'ID' in line:
						molid = line.strip().split('\t',)[9]
						molid = molid.replace (" ", "-")
						molid = molid.split (".",)
						print len(molid)
						if len(molid)>1:
							print "Sample name as per VCF file is "+molid[0]
							molid = molid[0]
						else :
							print "Sample name as per VCF file is "+molid
							molid = molid



			panel = chip[count] 
			samplename = "R"+run+"-"+panel+"_BC"+i.split("_")[1]+"_"+molid+sequencer			
			analyse_sample(i, vardir, panel, samplename)
			
		else:
			print "Please check the folder "+vardir
		count+=1	
		
		
	print os.path.abspath('')
	print os.system("echo $PATH > test.txt") #test

dirs = list()
runs = list()

localtime = time.asctime(time.localtime(time.time())) 
if os.path.exists(runsdb):
	oldcontents = open(runsdb, "r").read().splitlines()


#print run_dir

if "_MOL-" or "_SN2-" or "_S5XL-" in run_dir:
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



'''	#	if "SMCHD1" in dir:
#		all_vardirs = [d for d in os.listdir(pluginsdir) if 'variantCaller_out' in d]
#	print all_vardirs
	
	if len(all_vardirs) > 1:
		vardir = max([os.path.join(pluginsdir,all_vardirs) for all_vardirs in os.listdir(pluginsdir)], key=os.path.getmtime)
		print "vardir is "+vardir
		vardir = sorted(dirs.iteritems(), key=operator.itemgetter(1))
		print "vardir is "+vardir
	else:
		vardir = os.path.join(pluginsdir,all_vardirs[0])
		
	print "vardir is "+vardir '''	
	
	#if os.path.exists(watchdir+dir+"/plugin_out/variantCaller3.4.2_out"):
	#	tempdir = os.path.join(watchdir+dir+"/plugin_out/variantCaller3.4.2_out")
