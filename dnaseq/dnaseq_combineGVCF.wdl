#-------------------------------------------------------------------------------
# DnaSeq WDL workflow
# Version1 based on GenPipes 3.1.5-beta
# Created on: 2020-07-18
#-------------------------------------------------------------------------------

import "genpipes_tasks.wdl"


workflow combineGVCFScatter {

	Array[File] IN_VCFS_G
	Array[File] IN_VCF_G_INDEX

	String GENOME_FASTA
	Array[String] INTERVALS
	Array[String] CHR_EXCLUDE

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR

	String TMPDIR 



	call genpipes_tasks.array_extend as concat_exclude {

    	input:
    	list1 = CHR_EXCLUDE,
    	list2 = INTERVALS
    	
    }


  scatter(chr in INTERVALS){

    call genpipes_tasks.combine_gvcf as combine_gvcf_INCLD {

        input:
		IN_VCFS_G = IN_VCFS_G,
		IN_VCF_G_INDEX = IN_VCF_G_INDEX,
		GENOME_FASTA = GENOME_FASTA,
		INTERVALS = [chr],
		INTERVAL_NAME = chr,
		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

    }

  }

    call genpipes_tasks.combine_gvcf as combine_gvcf_EXCLD {

        input:
        IN_VCFS_G = IN_VCFS_G,
        IN_VCF_G_INDEX = IN_VCF_G_INDEX,
		GENOME_FASTA = GENOME_FASTA,
		CHR_EXCLUDE = concat_exclude.OUT,
		INTERVAL_NAME = "others",
		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

    }



    call genpipes_tasks.array_extend_file as concat_gvcf {

    	input:
    	list1 = combine_gvcf_INCLD.OUT_VCF_G,
    	list2 = [combine_gvcf_EXCLD.OUT_VCF_G]

    }

    call genpipes_tasks.array_extend_file as concat_gvcf_index {

    	input:
    	list1 = combine_gvcf_INCLD.OUT_VCF_G_INDEX,
    	list2 = [combine_gvcf_EXCLD.OUT_VCF_G_INDEX]

    }

    output {

      Array[File] OUT_GVCFs = concat_gvcf.OUT
      Array[File] OUT_GVCF_INDEXs = concat_gvcf_index.OUT

    }


}