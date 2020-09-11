#-------------------------------------------------------------------------------
# DnaSeq WDL workflow
# Version1 based on GenPipes 3.1.5-beta
# Created on: 2020-07-18
#-------------------------------------------------------------------------------



import "genpipes_tasks.wdl"
import "dnaseq_readsetScatter.wdl" as readsetScatter
import "dnaseq_gatkRealinger.wdl" as gatkRealingerScatter
import "dnaseq_gatkCaller.wdl" as gatkCallerScatter
import "dnaseq_combineGVCF.wdl" as combineGVCFScatter



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

			call genpipes_tasks.sambamba_merge_sam_files {

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

		}


		call genpipes_tasks.sambamba_merge_sam_files as sambamba_merge_sam_files_gatk {

			input:
			SAMPLE = sample.sample,
			PREFIX = ".realigned.sorted.bam",
			IN_BAMS = gatkRealingerScatter.OUT_BAMs,

			MOD_SAMTOOLS = MOD_SAMTOOLS,
			MOD_SAMBAMBA = MOD_SAMBAMBA

		}


		call genpipes_tasks.fix_mate_by_coordinate {

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


		call genpipes_tasks.picard_mark_duplicates {

			input:
			SAMPLE = sample.sample,
			IN_BAM = fix_mate_by_coordinate.OUT_BAM,

			MOD_JAVA = MOD_JAVA, 
			MOD_PICARD = MOD_PICARD,
			PICARD_HOME = PICARD_HOME,
			TMPDIR = TMPDIR

		}


		call genpipes_tasks.recalibration_report {

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


		call genpipes_tasks.recalibration {

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


		
		call genpipes_tasks.metrics_dna_picard_metrics_main {

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


		call genpipes_tasks.metrics_dna_picard_metrics_oxog {

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


		call genpipes_tasks.metrics_dna_picard_metrics_biasQc {

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


		call genpipes_tasks.metrics_dna_sample_qualimap {

			input:
			SAMPLE = sample.sample,
			IN_BAM = picard_mark_duplicates.OUT_BAM,

			SPECIES = SPECIES,

			MOD_JAVA = MOD_JAVA,
			MOD_QUALIMAP = MOD_QUALIMAP

		}



		call genpipes_tasks.metrics_dna_sambamba_flagstat {

			input:
			SAMPLE = sample.sample,
			IN_BAM = picard_mark_duplicates.OUT_BAM,

			MOD_SAMBAMBA = MOD_SAMBAMBA

		}


		call genpipes_tasks.metrics_dna_fastqc {

			input:
			SAMPLE = sample.sample,
			IN_BAM = picard_mark_duplicates.OUT_BAM,
			ADAPTER1 = readsetScatter.ADPTER1[0],
			ADAPTER2 = readsetScatter.ADPTER2[0],

			MOD_FASTQC = MOD_FASTQC,
			MOD_JAVA = MOD_JAVA

		}


		call genpipes_tasks.gatk_callable_loci	{

			input:
			SAMPLE = sample.sample,
			IN_BAM = recalibration.OUT_BAM,
			
			GENOME_FASTA = GENOME_FASTA,

			MOD_JAVA = MOD_JAVA,
			MOD_GATK = MOD_GATK,
			GATK_JAR = GATK_JAR,
			TMPDIR = TMPDIR

		}


		call genpipes_tasks.extract_common_snp_freq {

			input:
			SAMPLE = sample.sample,
			IN_BAM = recalibration.OUT_BAM,
			
			DB_SNP_COMMON = DB_SNP_COMMON,

			MOD_JAVA = MOD_JAVA,
			MOD_BVATOOLS = MOD_BVATOOLS,
			BVATOOLS_JAR = BVATOOLS_JAR

		}	


		call genpipes_tasks.baf_plot {

			input:
			SAMPLE = sample.sample,

			IN_COMMONSNPS = extract_common_snp_freq.OUT_COMMONSNPS,

			DB_SNP_COMMON = DB_SNP_COMMON,
			GENOME_DICT = GENOME_DICT,
			CHR_EXCLUDE = CHR_EXCLUDE,

			MOD_JAVA = MOD_JAVA,
			MOD_BVATOOLS = MOD_BVATOOLS,
			BVATOOLS_JAR = BVATOOLS_JAR

		}


		call genpipes_tasks.cram {

			input:
			SAMPLE = sample.sample,
			IN_BAM = recalibration.OUT_BAM,

			GENOME_FASTA = GENOME_FASTA,

			MOD_SAMTOOLS = MOD_SAMTOOLS

		}


		call gatkCallerScatter.gatkCallerScatter {

			input:
			SAMPLE = sample.sample,
			IN_BAM = recalibration.OUT_BAM,
			INTERVALS = CHRS,
			CHR_EXCLUDE = CHR_EXCLUDE,
			GENOME_FASTA = GENOME_FASTA,
			G1000 = G1000,

			TMPDIR = TMPDIR,
			MOD_JAVA = MOD_JAVA,
			MOD_GATK = MOD_GATK,
			GATK_JAR = GATK_JAR

		}


		call genpipes_tasks.merge_and_call_individual_gvcf_merge {

			input:
			SAMPLE = sample.sample,
			IN_VCFS_INTRVL = gatkCallerScatter.OUT_VCFs,
			GENOME_FASTA = GENOME_FASTA,
			TMPDIR = TMPDIR,
			MOD_JAVA = MOD_JAVA,
			MOD_GATK = MOD_GATK,
			GATK_JAR = GATK_JAR


		}


		call genpipes_tasks.merge_and_call_individual_gvcf_calls {

			input:
			SAMPLE = sample.sample,
			IN_VCF_G = merge_and_call_individual_gvcf_merge.OUT_VCF_G,
			IN_VCF_G_INDEX = merge_and_call_individual_gvcf_merge.OUT_VCF_G_INDEX,

			GENOME_FASTA = GENOME_FASTA,
		
			TMPDIR = TMPDIR,
			MOD_JAVA = MOD_JAVA,
			MOD_GATK = MOD_GATK,
			GATK_JAR = GATK_JAR

		}


	}


	call combineGVCFScatter.combineGVCFScatter {

		input:
		IN_VCFS_G = merge_and_call_individual_gvcf_merge.OUT_VCF_G,
		IN_VCF_G_INDEX = merge_and_call_individual_gvcf_merge.OUT_VCF_G_INDEX,

		GENOME_FASTA = GENOME_FASTA,
		INTERVALS = CHRS,
		CHR_EXCLUDE = CHR_EXCLUDE,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

	}


	call genpipes_tasks.merge_combined_gvcf {
		
		input:
		VCFS = combineGVCFScatter.OUT_GVCFs,
		VCF_INDEXs = combineGVCFScatter.OUT_GVCF_INDEXs,

		GENOME_FASTA = GENOME_FASTA,
	
		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

	}



	call genpipes_tasks.call_combined_gvcf {

		input:
		IN_VCF_G = merge_combined_gvcf.OUT_VCF_G,
		IN_VCF_G_INDEX = merge_combined_gvcf.OUT_VCF_G_INDEX,

		GENOME_FASTA = GENOME_FASTA,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

	}


	call genpipes_tasks.variant_recalibrator_prep {

		input:
		IN_VCF_ALL = call_combined_gvcf.OUT_VCF_ALL,
		IN_VCF_ALL_INDEX = call_combined_gvcf.OUT_VCF_ALL_INDEX,

		GENOME_FASTA = GENOME_FASTA,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR,
		MOD_R = MOD_R


	}


	call genpipes_tasks.variant_recalibrator_exec {

		input:
		IN_VCF_ALL = call_combined_gvcf.OUT_VCF_ALL,
		IN_VCF_ALL_INDEX = call_combined_gvcf.OUT_VCF_ALL_INDEX,
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


	call genpipes_tasks.haplotype_caller_decompose_and_normalize {

		input:
		IN_VQSR = variant_recalibrator_exec.OUT_VQSR,
		IN_VQSR_INDEX = variant_recalibrator_exec.OUT_VQSR_INDEX,

		GENOME_FASTA = GENOME_FASTA,

		MOD_HTSLIB = MOD_HTSLIB,
		MOD_VT = MOD_VT

	}


	call genpipes_tasks.haplotype_caller_flag_mappability {

		input:
		IN_VT = haplotype_caller_decompose_and_normalize.OUT_VT,
		IN_VT_INDEX = haplotype_caller_decompose_and_normalize.OUT_VT_INDEX,

		IlluEXCLUSION = IlluEXCLUSION,

		MOD_VCFTOOLS = MOD_VCFTOOLS,
		MOD_TABIX = MOD_TABIX,
		MOD_HTSLIB = MOD_HTSLIB

	}


	call genpipes_tasks.haplotype_caller_snp_id_annotation {

		input:
		IN_MIL = haplotype_caller_flag_mappability.OUT_MIL,
		IN_MIL_INDEX = haplotype_caller_flag_mappability.OUT_MIL_INDEX,

		DB_SNP = DB_SNP,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_SNPEFF = MOD_SNPEFF,
		SNPEFF_HOME = SNPEFF_HOME,
		MOD_HTSLIB = MOD_HTSLIB,
		MOD_TABIX = MOD_TABIX

	}


	call genpipes_tasks.haplotype_caller_snp_effect {

		input:
		IN_SNPID = haplotype_caller_snp_id_annotation.OUT_SNPID,
		IN_SNPID_INDEX = haplotype_caller_snp_id_annotation.OUT_SNPID_INDEX,

		ASSEMBLY = ASSEMBLY,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_SNPEFF = MOD_SNPEFF,
		SNPEFF_HOME = SNPEFF_HOME,
		MOD_HTSLIB = MOD_HTSLIB,
		MOD_TABIX = MOD_TABIX

	}


	call genpipes_tasks.haplotype_caller_dbnsfp_annotation {

		input:
		IN_ZIPPED = haplotype_caller_snp_effect.OUT_ZIPPED,
		IN_ZIPPED_INDEX = haplotype_caller_snp_effect.OUT_ZIPPED_INDEX,

		DB_SFP = DB_SFP,

		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_SNPEFF = MOD_SNPEFF,
		SNPEFF_HOME = SNPEFF_HOME,
		MOD_HTSLIB = MOD_HTSLIB,
		MOD_TABIX = MOD_TABIX

	}


	call genpipes_tasks.haplotype_caller_gemini_annotations {

		input:
		IN_ZIPPED = haplotype_caller_dbnsfp_annotation.OUT_ZIPPED,
		IN_ZIPPED_INDEX = haplotype_caller_dbnsfp_annotation.OUT_ZIPPED_INDEX,

		TMPDIR = TMPDIR,
		MOD_GEMINI = MOD_GEMINI,
		MOD_HTSLIB = MOD_HTSLIB

	}


	call genpipes_tasks.haplotype_caller_metrics_vcf_stats {

		input:
		IN_ZIPPED = haplotype_caller_dbnsfp_annotation.OUT_ZIPPED,
		IN_ZIPPED_INDEX = haplotype_caller_dbnsfp_annotation.OUT_ZIPPED_INDEX,
		IN_STATS = haplotype_caller_snp_effect.OUT_STATS,

		GENOME_DICT = GENOME_DICT,

		MOD_PYTHON = MOD_PYTHON,
		PYTHON_TOOLS = PYTHON_TOOLS,
		MOD_MUGQICTOOLS = MOD_MUGQICTOOLS

	}


	call genpipes_tasks.run_multiqc {

		input:
		## add files from every step
		METRICS = [haplotype_caller_metrics_vcf_stats.OUT_CHGRATE],
		MOD_MULTIQC = MOD_MULTIQC

	}


}













