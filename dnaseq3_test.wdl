#-------------------------------------------------------------------------------
# DnaSeq WDL workflow
# Version1 based on GenPipes 3.1.5-beta
# Created on: 2020-06-18
#-------------------------------------------------------------------------------


import "genpipes_tasks.wdl"


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


	String IN_BAM = "/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/alignment/S2/S2.sorted.dup.recal.bam"
	String IN_VCF = "/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/variants/allSamples.hc.vcf.gz"
	String IN_VQSR = "/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/variants/allSamples.hc.vqsr.vcf"
	String IN_ZIPPED = "/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.gz"
	String IN_STATS = "/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.stats.csv"







}













