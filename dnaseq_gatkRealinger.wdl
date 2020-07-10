#-------------------------------------------------------------------------------
# DnaSeq WDL workflow
# Version1 based on GenPipes 3.1.5-beta
# Created on: 2020-06-18
#-------------------------------------------------------------------------------

import "dnaseq_tasks.wdl"


workflow gatkRealingerScatter {

 
	String SAMPLE
	File IN_BAM
	Array[String] INTERVALS

	File GENOME_FASTA
	File G1000

	String TMPDIR
	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR




  scatter(chr in INTERVALS){

    call dnaseq_tasks.gatk_indel_realigner as gatk_indel_realigner_INCLD {

        input:
		SAMPLE = SAMPLE,
		IN_BAM = IN_BAM,
		INTERVALS = [chr],
		INTERVAL_NAME = chr,
		GENOME_FASTA = GENOME_FASTA,
		G1000 = G1000,
		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

    }

  }

    call dnaseq_tasks.gatk_indel_realigner as gatk_indel_realigner_EXCLD {

        input:
		SAMPLE = SAMPLE,
		IN_BAM = IN_BAM,
		EXCLUDE = INTERVALS,
		INTERVAL_NAME = "others",
		GENOME_FASTA = GENOME_FASTA,
		G1000 = G1000,
		TMPDIR = TMPDIR,
		MOD_JAVA = MOD_JAVA,
		MOD_GATK = MOD_GATK,
		GATK_JAR = GATK_JAR

    }



    call dnaseq_tasks.concat_arrays as concat_intvls {

    	input:
    	A = gatk_indel_realigner_INCLD.OUT_INTRVL,
    	B = [gatk_indel_realigner_EXCLD.OUT_INTRVL]

    }

    call dnaseq_tasks.concat_arrays as concat_bams {

    	input:
    	A = gatk_indel_realigner_INCLD.OUT_BAM,
    	B = [gatk_indel_realigner_EXCLD.OUT_BAM]
    	
    }

    output {

      Array[File] OUT_INTRVLs = concat_intvls.OUT
      Array[File] OUT_BAMs = concat_bams.OUT

    }


}