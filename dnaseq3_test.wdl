#-------------------------------------------------------------------------------
# DnaSeq WDL workflow
# Version1 based on GenPipes 3.1.5-beta
# Created on: 2020-06-18
#-------------------------------------------------------------------------------


import "dnaseq_tasks.wdl"


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
	String IN_VCF = ""

	call dnaseq_tasks.merge_and_call_combined_gvcf_calls {

		input:
		#IN_VCF_G = merge_and_call_combined_gvcf_merges.OUT_VCF_G,
		IN_VCF_G = IN_VCF,

		GENOME_FASTA = GENOME_FASTA,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

	}

	call dnaseq_tasks.variant_recalibrator_prep {

		input:
		IN_VCF_ALL = merge_and_call_combined_gvcf_calls.OUT_VCF_ALL,

		GENOME_FASTA = GENOME_FASTA,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR,
		MOD_R = MOD_R


	}

	call dnaseq_tasks.variant_recalibrator_exec {

		input:
		IN_VCF_ALL = merge_and_call_combined_gvcf_calls.OUT_VCF_ALL,
		IN_SNP_RCL = variant_recalibrator_prep.OUT_SNP_RCL,
		IN_SNP_TRNCH = variant_recalibrator_prep.OUT_SNP_TRNCH,
		IN_SNP_R = variant_recalibrator_prep.OUT_SNP_R,
		IN_INDL_RCL = variant_recalibrator_prep.OUT_INDL_RCL,
		IN_INDL_TRNCH = variant_recalibrator_prep.OUT_INDL_TRNCH,
		IN_INDL_R = variant_recalibrator_prep.OUT_INDL_R,

		GENOME_FASTA = GENOME_FASTA,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

	}

	call dnaseq_tasks.haplotype_caller_decompose_and_normalize {

		input:
		IN_VQSR = variant_recalibrator_exec.OUT_VQSR,

		GENOME_FASTA = GENOME_FASTA,

		MOD_HTSLIB = MOD_HTSLIB,
		MOD_VT = MOD_VT

	}


	call dnaseq_tasks.haplotype_caller_flag_mappability {

		input:
		IN_VT = haplotype_caller_decompose_and_normalize.OUT_VT,

		IlluEXCLUSION = IlluEXCLUSION,

		MOD_VCFTOOLS = MOD_VCFTOOLS,
		MOD_TABIX = MOD_TABIX,
		MOD_HTSLIB = MOD_HTSLIB

	}


	call dnaseq_tasks.haplotype_caller_snp_id_annotation {

		input:
		IN_MIL = haplotype_caller_flag_mappability.OUT_MIL,

		DB_SNP = DB_SNP,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_SNPEFF = MOD_SNPEFF,
		SNPEFF_HOME = SNPEFF_HOME,
		MOD_HTSLIB = MOD_HTSLIB,
		MOD_TABIX = MOD_TABIX

	}


	call dnaseq_tasks.haplotype_caller_snp_effect {

		input:
		IN_SNPID = haplotype_caller_snp_id_annotation.OUT_SNPID,

		ASSEMBLY = ASSEMBLY,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_SNPEFF = MOD_SNPEFF,
		SNPEFF_HOME = SNPEFF_HOME,
		MOD_HTSLIB = MOD_HTSLIB,
		MOD_TABIX = MOD_TABIX

	}

	call dnaseq_tasks.haplotype_caller_dbnsfp_annotation {

		input:
		IN_ZIPPED = haplotype_caller_snp_effect.OUT_ZIPPED,

		DB_SFP = DB_SFP,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_SNPEFF = MOD_SNPEFF,
		SNPEFF_HOME = SNPEFF_HOME,
		MOD_HTSLIB = MOD_HTSLIB,
		MOD_TABIX = MOD_TABIX

	}

	call dnaseq_tasks.haplotype_caller_gemini_annotations {

		input:
		IN_ZIPPED = haplotype_caller_dbnsfp_annotation.OUT_ZIPPED,

		TMPDIR = TMPDIR,
		MOD_GEMINI = MOD_GEMINI,
		MOD_HTSLIB = MOD_HTSLIB

	}

	call dnaseq_tasks.haplotype_caller_metrics_vcf_stats {

		input:
		IN_ZIPPED = haplotype_caller_dbnsfp_annotation.OUT_ZIPPED,
		IN_STATS = haplotype_caller_snp_effect.OUT_STATS,

		GENOME_DICT = GENOME_DICT,

		MOD_PYTHON = MOD_PYTHON,
		PYTHON_TOOLS = PYTHON_TOOLS,
		MOD_MUGQICTOOLS = MOD_MUGQICTOOLS

	}

	call dnaseq_tasks.run_multiqc {

		input:
		## add files from every step
		METRICS = [haplotype_caller_metrics_vcf_stats.OUT_CHGRATE],
		MOD_MULTIQC = MOD_MULTIQC

	}



}













