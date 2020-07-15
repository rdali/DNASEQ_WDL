#-------------------------------------------------------------------------------
# DnaSeq WDL workflow
# Version1 based on GenPipes 3.1.5-beta
# Created on: 2020-06-18
#-------------------------------------------------------------------------------


import "dnaseq_tasks.wdl"
import "dnaseq_readsetScatter.wdl" as readsetScatter
import "dnaseq_gatkCaller.wdl" as gatkCallerScatter


#-------------------------------------------------------------------------------
# Workflow
#-------------------------------------------------------------------------------

workflow DnaSeq {

	## samples:

	Array[Object] samples

	## Species VARS:
	String SPECIES
	String ASSEMBLY
	Array[String] CHRS
	Array[String] CHR_EXCLUDE
	File GENOME_FASTA
	String BWA_INDX
	File GENOME_DICT
	File DB_SNP
	File DB_SNP_COMMON
	File DB_SFP
	File G1000
	File HAPMAP
	File OMNI
	File GNOMAD
	File IlluEXCLUSION


	## ENV VARS:
	String TMPDIR

	## MODULE VARS
	String MOD_PICARD
	String PICARD_HOME
	String MOD_JAVA
	String MOD_TRIMMOMATIC
	String TRIMMOMATIC_JAR
	String MOD_BWA
	String MOD_SAMTOOLS
	String MOD_SAMBAMBA
	String MOD_BVATOOLS
	String BVATOOLS_JAR
	String MOD_GATK
	String GATK_JAR
	String MOD_R
	String MOD_QUALIMAP
	String MOD_FASTQC
	String MOD_HTSLIB
	String MOD_VT
	String MOD_VCFTOOLS
	String MOD_TABIX
	String MOD_SNPEFF
	String SNPEFF_HOME
	String MOD_GEMINI
	String MOD_PYTHON
	String MOD_MUGQICTOOLS
	String PYTHON_TOOLS
	String MOD_MULTIQC

	## loop over samples then readsets

	scatter (sample in samples) {



	String IN_BAM = "/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/alignment/S2/S2.sorted.dup.recal.bam"
	String IN_VCF_G = "/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/alignment/S2/S2.hc.g.vcf.gz"

	


	#call combine_gvcf


	#output {


#	}

}

#	call dnaseq_tasks.merge_and_call_combined_gvcf_merges {
		
#		input:
#		VCFS = ,
#		GENOME_FASTA = GENOME_FASTA,
#	
#		TMPDIR = TMPDIR,
#		MOD_JAVA = MOD_JAVA,
#		MOD_GATK = MOD_GATK,
#		GATK_JAR = GATK_JAR

#	}



}













