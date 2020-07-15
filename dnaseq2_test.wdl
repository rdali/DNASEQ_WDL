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


	

	call gatkCallerScatter.gatkCallerScatter {

		input:
		SAMPLE = sample.sample,
		#IN_BAM = recalibration.OUT_BAM,
		IN_BAM = IN_BAM,
		INTERVALS = CHRS,
		CHR_EXCLUDE = CHR_EXCLUDE,
		GENOME_FASTA = GENOME_FASTA,
		G1000 = G1000,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

		#output:
#      Array[File] OUT_VCFs
#      Array[File] OUT_TBIs


	}

	call dnaseq_tasks.merge_and_call_individual_gvcf_merge {

		input:
		SAMPLE = sample.sample,
		IN_VCFS_INTRVL = gatkCallerScatter.OUT_VCFs,
		GENOME_FASTA = GENOME_FASTA,
		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR


	}

	call dnaseq_tasks.merge_and_call_individual_gvcf_calls {

		input:
		SAMPLE = sample.sample,
		IN_VCF_G = merge_and_call_individual_gvcf_merge.OUT_VCF_G,

		GENOME_FASTA = GENOME_FASTA,
	
		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

	}

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













