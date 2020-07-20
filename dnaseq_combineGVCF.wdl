#-------------------------------------------------------------------------------
# DnaSeq WDL workflow
# Version1 based on GenPipes 3.1.5-beta
# Created on: 2020-07-18
#-------------------------------------------------------------------------------

import "genpipes_tasks.wdl"


workflow combineGVCFScatter {

	Array[File] IN_VCFS_G

	String GENOME_FASTA
	Array[String] INTERVALS
	Array[String] CHR_EXCLUDE

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR

	String TMPDIR 



  scatter(chr in INTERVALS){

    call genpipes_tasks.combine_gvcf as combine_gvcf_INCLD {

        input:
		IN_VCFS_G = IN_VCFS_G,
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
		GENOME_FASTA = GENOME_FASTA,
		CHR_EXCLUDE = CHR_EXCLUDE,
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


    output {

      Array[File] OUT_GVCFs = concat_gvcf.OUT

    }


}