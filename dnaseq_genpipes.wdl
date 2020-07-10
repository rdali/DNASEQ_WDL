#-------------------------------------------------------------------------------
# DnaSeq WDL workflow
# Version1 based on GenPipes 3.1.5-beta
# Created on: 2020-06-18
#-------------------------------------------------------------------------------


import "dnaseq_tasks.wdl"
import "dnaseq_readsetScatter.wdl" as readsetScatter
import "dnaseq_gatkRealinger.wdl" as gatkRealingerScatter


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

		if (length(readsetScatter.OUT_BAMs) > 1){

			call dnaseq_tasks.sambamba_merge_sam_files {

				input:

				SAMPLE = sample.sample,
				PREFIX = ".sorted.bam",
				IN_BAMS = readsetScatter.OUT_BAMs,

				MOD_SAMTOOLS = MOD_SAMTOOLS,
				MOD_SAMBAMBA = MOD_SAMBAMBA

			}
		}


		call gatkRealingerScatter.gatkRealingerScatter {

			input:
			SAMPLE = sample.sample,
			IN_BAM = select_first([sambamba_merge_sam_files.OUT_MERGED_BAM, if (length(readsetScatter.OUT_BAMs) == 1) then readsetScatter.OUT_BAMs[0] else "" ]),
			INTERVALS = CHRS,
			GENOME_FASTA = GENOME_FASTA,
			G1000 = G1000,
			TMPDIR = TMPDIR,
			MOD_JAVA = MOD_JAVA,
			MOD_GATK = MOD_GATK,
			GATK_JAR = GATK_JAR

		#output:
#      Array[File] OUT_INTRVLs = concat_intvls.OUT
#      Array[File] OUT_BAMs = concat_bams.OUT


		}


		call dnaseq_tasks.sambamba_merge_sam_files as sambamba_merge_sam_files_gatk {

			input:
			SAMPLE = sample.sample,
			PREFIX = ".realigned.sorted.bam",
			IN_BAMS = gatkRealingerScatter.OUT_BAMs,

			MOD_SAMTOOLS = MOD_SAMTOOLS,
			MOD_SAMBAMBA = MOD_SAMBAMBA

		}


		call dnaseq_tasks.fix_mate_by_coordinate {

			input:
			SAMPLE = sample.sample,
			IN_BAM = sambamba_merge_sam_files_gatk.OUT_MERGED_BAM,

			TMPDIR = TMPDIR,
			MOD_JAVA = MOD_JAVA,
			MOD_BVATOOLS = MOD_BVATOOLS,
			BVATOOLS_JAR = BVATOOLS_JAR,
			MOD_PICARD = MOD_PICARD,
			PICARD_HOME = PICARD_HOME

		}

	}
}













