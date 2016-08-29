#!/bin/bash

# new via pipeline
# ~/Myscripts/via_pipeline_new_temp.sh /Volumes/Iontorrent/TSv5.0_validation/TVC_5.0/S5XL-AmpliSeq_Controls/R_2016_05_17_11_24_33_user_S5XL-00642-11-MPCP

echo "Analyzing VCf file"

echo $1
bcftools norm -m -both -f /Users/molpathuser1/othertools/cpipe/hg19/hg19.fasta $1 -o ${1/.vcf.gz/_norm.vcf.gz};
bcftools view -i 'FORMAT/GT=="0/1" || FORMAT/GT=="1/0" || FORMAT/GT=="1/1"' ${1/.vcf.gz/_norm.vcf.gz} -o ${1/.vcf.gz/_norm_filt.vcf.gz};
cd /Users/molpathuser1/othertools/ensembl-tools-release-84/scripts/variant_effect_predictor
perl /Users/molpathuser1/othertools/ensembl-tools-release-84/scripts/variant_effect_predictor/variant_effect_predictor.pl -i ${1/.vcf.gz/_norm_filt.vcf.gz} \
	--offline --cache --dir_cache ~/othertools/vep/databases --refseq \
	--fasta ~/othertools/vep/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
	--everything --check_existing --no_escape\
	--dir_plugins ~/othertools/vep/databases/Plugins --plugin dbNSFP,/Users/molpathuser1/othertools/vep/databases/dbNSFP.gz \
	--force_overwrite --vcf -o STDOUT -fork 4 --flag_pick_allele > ${1/.vcf.gz/_norm_filt_vep_out.vcf}

grep "^#" ${1/.vcf.gz/_norm_filt_vep_out.vcf} > ${1/.vcf.gz/_norm_filt_vep_out_NM.vcf}
cut -f2 /Volumes/Iontorrent/VIA_pipeline/gene_refSeq.txt | sort | uniq | while read LINE ; 
	do 
	echo ${LINE/.?/}; \
	perl /Users/molpathuser1/othertools/ensembl-tools-release-84/scripts/variant_effect_predictor/filter_vep.pl -filter "Feature match ${LINE/.?/}" \
	--only_matched -i ${1/.vcf.gz/_norm_filt_vep_out.vcf} | grep -v "^#" >> ${1/.vcf.gz/_norm_filt_vep_out_NM.vcf}; 
	done
	
echo $'CHROM\tPOS\tREF\tALT\tAllele Freq\tQUAL\tCosmic ID\tTYPE\tDepth\tRef Cov\tAlt Cov\tAlt Fwd Cov\tAlt Rev Cov\tAllele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\t\Amino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tPICK\tVARIANT_CLASS\tSYMBOL_SOURCE\tHGNC_ID\tCANONICAL\tTSL\tAPPRIS\tCCDS\tENSP\tSWISSPROT\tTREMBL\tUNIPARC\tREFSEQ_MATCH\tGENE_PHENO\tSIFT\tPolyPhen\tDOMAINS\tHGVS_OFFSET\tGMAF\tAFR_MAF\tAMR_MAF\tEAS_MAF\tEUR_MAF\tSAS_MAF\tAA_MAF\tEA_MAF\tExAC_MAF\tExAC_Adj_MAF\tExAC_AFR_MAF\tExAC_AMR_MAF\tExAC_EAS_MAF\tExAC_FIN_MAF\tExAC_NFE_MAF\tExAC_OTH_MAF\tExAC_SAS_MAF\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tMOTIF_NAME\tMOTIF_POS\tHIGH_INF_POS\tMOTIF_SCORE_CHANGE' > ${1/.vcf.gz/_norm_filt_vep_out_NM.tsv}
bcftools query   -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%QUAL\t%OID\t%TYPE\t%DP\t%RO\t%AO\t%SAF\t%SAR\t%CSQ\n' ${1/.vcf.gz/_norm_filt_vep_out_NM.vcf} | awk 'gsub(/\|/,"\t") gsub(/&/,", ")' >> ${1/.vcf.gz/_norm_filt_vep_out_NM.tsv}
grep "^CHROM" ${1/.vcf.gz/_norm_filt_vep_out_NM.tsv} > ${1/.vcf.gz/_prefiltered.tsv}
grep -v "^CHROM" ${1/.vcf.gz/_norm_filt_vep_out_NM.tsv} | sort | uniq >> ${1/.vcf.gz/_prefiltered.tsv}
rm ${1/TSVC_variants.vcf.gz/*_norm*vcf*}

# Key files - variant_effect_predictor.pl, hg19.fasta, gene_refSeq
