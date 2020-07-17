#-------------------------------------------------------------------------------
# DnaSeq WDL workflow
# Version1 based on GenPipes 3.1.5-beta
# Created on: 2020-07-18
#-------------------------------------------------------------------------------

import "genpipes_tasks.wdl"


workflow gatkCallerScatter {

 
	String SAMPLE
	File IN_BAM
	Array[String] INTERVALS
	Array[String] CHR_EXCLUDE

	File GENOME_FASTA
	File G1000

	String TMPDIR
	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR



	call genpipes_tasks.array_extend as concat_exclude {

    	input:
    	list1 = CHR_EXCLUDE,
    	list2 = INTERVALS
    	
    }


  scatter(chr in INTERVALS){

    call genpipes_tasks.gatk_haplotype_caller as gatk_haplotype_caller_INCLD {

        input:
		SAMPLE = SAMPLE,
		IN_BAM = IN_BAM,
		INTERVALS = [chr],
		INTERVAL_NAME = chr,
		GENOME_FASTA = GENOME_FASTA,
		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

    }

  }

    call genpipes_tasks.gatk_haplotype_caller as gatk_haplotype_caller_EXCLD {

        input:
		SAMPLE = SAMPLE,
		IN_BAM = IN_BAM,
		CHR_EXCLUDE = concat_exclude.OUT,
		INTERVAL_NAME = "others",
		GENOME_FASTA = GENOME_FASTA,
		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

    }



    call genpipes_tasks.array_extend_file as concat_vcf {

    	input:
    	list1 = gatk_haplotype_caller_INCLD.OUT_VCF_INTRVL,
    	list2 = [gatk_haplotype_caller_EXCLD.OUT_VCF_INTRVL]

    }

    call genpipes_tasks.array_extend_file as concat_tbi {

    	input:
    	list1 = gatk_haplotype_caller_INCLD.OUT_TBI_INTRVL,
    	list2 = [gatk_haplotype_caller_EXCLD.OUT_TBI_INTRVL]
    	
    }

    output {

      Array[File] OUT_VCFs = concat_vcf.OUT
      Array[File] OUT_TBIs = concat_tbi.OUT

    }


}