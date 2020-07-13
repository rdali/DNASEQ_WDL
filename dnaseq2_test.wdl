#-------------------------------------------------------------------------------
# DnaSeq WDL workflow
# Version1 based on GenPipes 3.1.5-beta
# Created on: 2020-06-18
#-------------------------------------------------------------------------------


import "dnaseq_tasks.wdl"
import "dnaseq_readsetScatter.wdl" as readsetScatter


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

		call readsetScatter.readsetScatter {

			input:
			SAMPLE = sample.sample,
			READSETS = sample.readsets,

			BWA_INDX = BWA_INDX,

			TMPDIR = TMPDIR,
			MOD_JAVA = MOD_JAVA,
			MOD_PICARD = MOD_PICARD,
			PICARD_HOME = PICARD_HOME, 
			MOD_TRIMMOMATIC = MOD_TRIMMOMATIC,
			TRIMMOMATIC_JAR = TRIMMOMATIC_JAR,
			MOD_JAVA = MOD_JAVA,
			MOD_BWA = MOD_BWA

		}
	

		call dnaseq_tasks.fix_mate_by_coordinate {

			input:
			SAMPLE = sample.sample,
			#IN_BAM = sambamba_merge_sam_files_gatk.OUT_MERGED_BAM,
			IN_BAM="/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/alignment/S1/S1.realigned.sorted.bam",

			TMPDIR = TMPDIR,
			MOD_JAVA = MOD_JAVA,
			MOD_BVATOOLS = MOD_BVATOOLS,
			BVATOOLS_JAR = BVATOOLS_JAR,
			MOD_PICARD = MOD_PICARD,
			PICARD_HOME = PICARD_HOME

		}

		call dnaseq_tasks.picard_mark_duplicates {

			input:
			SAMPLE = sample.sample,
			IN_BAM = fix_mate_by_coordinate.OUT_BAM,

			MOD_JAVA = MOD_JAVA, 
			MOD_PICARD = MOD_PICARD,
			PICARD_HOME = PICARD_HOME,
			TMPDIR = TMPDIR

		}

		call dnaseq_tasks.recalibration_report {

			input:
			SAMPLE = sample.sample,
			IN_BAM = picard_mark_duplicates.OUT_BAM,
			GENOME_FASTA = GENOME_FASTA,
			KNOWNSITES = [DB_SNP, G1000, GNOMAD],

			MOD_JAVA = MOD_JAVA,
			MOD_GATK = MOD_GATK,
			GATK_JAR = GATK_JAR,
			TMPDIR = TMPDIR

		}

		call dnaseq_tasks.recalibration {

			input:
			SAMPLE = sample.sample,
			IN_BAM = picard_mark_duplicates.OUT_BAM,

			GENOME_FASTA = GENOME_FASTA,
			IN_CLB_RPT  = recalibration_report.OUT_CLB_RPT,

			MOD_JAVA = MOD_JAVA,
			MOD_GATK = MOD_GATK,
			GATK_JAR = GATK_JAR,
			MOD_SAMTOOLS = MOD_SAMTOOLS,
			MOD_SAMBAMBA = MOD_SAMBAMBA,
			TMPDIR = TMPDIR

		}

		call dnaseq_tasks.metrics_dna_picard_metrics_main {

			input:
			SAMPLE = sample.sample,
			IN_BAM = picard_mark_duplicates.OUT_BAM,

			GENOME_FASTA = GENOME_FASTA,

			MOD_JAVA = MOD_JAVA,
			MOD_PICARD = MOD_PICARD,
			PICARD_HOME = PICARD_HOME,
			MOD_R = MOD_R,
			TMPDIR = TMPDIR

		}


		call dnaseq_tasks.metrics_dna_picard_metrics_oxog {

			input:
			SAMPLE = sample.sample,
			IN_BAM = picard_mark_duplicates.OUT_BAM,

			DB_SNP = DB_SNP,
			GENOME_FASTA = GENOME_FASTA,

			MOD_JAVA = MOD_JAVA,
			MOD_PICARD = MOD_PICARD,
			PICARD_HOME = PICARD_HOME,
			MOD_R = MOD_R,
			TMPDIR = TMPDIR

		}


		call dnaseq_tasks.metrics_dna_picard_metrics_biasQc {

			input:
			SAMPLE = sample.sample,
			IN_BAM = picard_mark_duplicates.OUT_BAM,

			GENOME_FASTA = GENOME_FASTA,

			MOD_JAVA = MOD_JAVA,
			MOD_PICARD = MOD_PICARD,
			PICARD_HOME = PICARD_HOME,
			MOD_R = MOD_R,
			TMPDIR = TMPDIR

		}


		call dnaseq_tasks.metrics_dna_sample_qualimap {

			input:
			SAMPLE = sample.sample,
			IN_BAM = picard_mark_duplicates.OUT_BAM,

			SPECIES = SPECIES,

			MOD_JAVA = MOD_JAVA,
			MOD_QUALIMAP = MOD_QUALIMAP

		}


		call dnaseq_tasks.metrics_dna_sambamba_flagstat {

			input:
			SAMPLE = sample.sample,
			IN_BAM = picard_mark_duplicates.OUT_BAM,

			MOD_SAMBAMBA = MOD_SAMBAMBA

		}


		call dnaseq_tasks.metrics_dna_fastqc {

			input:
			SAMPLE = sample.sample,
			IN_BAM = picard_mark_duplicates.OUT_BAM,			
			ADAPTER1 = readsetScatter.ADPTER1[0],
			ADAPTER2 = readsetScatter.ADPTER2[0],

			MOD_FASTQC = MOD_FASTQC,
			MOD_JAVA = MOD_JAVA

		}


		call dnaseq_tasks.gatk_callable_loci	{

			input:
			SAMPLE = sample.sample,
			IN_BAM = recalibration.OUT_BAM,

			GENOME_FASTA = GENOME_FASTA,

			MOD_JAVA = MOD_JAVA,
			MOD_GATK = MOD_GATK,
			GATK_JAR = GATK_JAR,
			TMPDIR = TMPDIR

		}


		call dnaseq_tasks.extract_common_snp_freq {

			input:
			SAMPLE = sample.sample,
			IN_BAM = recalibration.OUT_BAM,

			DB_SNP_COMMON = DB_SNP_COMMON,

			MOD_JAVA = MOD_JAVA,
			MOD_BVATOOLS = MOD_BVATOOLS,
			BVATOOLS_JAR = BVATOOLS_JAR

		}	


		call dnaseq_tasks.baf_plot {

			input:
			SAMPLE = sample.sample,

			IN_COMMONSNPS = extract_common_snp_freq.OUT_COMMONSNPS,

			DB_SNP_COMMON = DB_SNP_COMMON,
			GENOME_DICT = GENOME_DICT,

			MOD_JAVA = MOD_JAVA,
			MOD_BVATOOLS = MOD_BVATOOLS,
			BVATOOLS_JAR = BVATOOLS_JAR

		}


	}
}













