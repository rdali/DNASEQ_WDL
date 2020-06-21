#-------------------------------------------------------------------------------
# DnaSeq WDL workflow
# Version: 3.1.5-beta
# Created on: 2020-06-18
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Workflow
#-------------------------------------------------------------------------------

workflow DnaSeq {

	## samples:
	File ReadsetFile
    Array[Array[File]] ReadsetFileContent = read_tsv(ReadsetFile)

    ## Species VARS:
    String SPECIES
    String ASSEMBLY
    File GENOME_FASTA
    File GENOME_DICT
    File DB_SNP
    File DB_SNP_COMMON\
    File DB_SFP
    File IlluEXCLUSION

	## ENV VARS:
	String OUT_DIR
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

	## loop over readset

    scatter (ReadsetFileLine in ReadsetFileContent) {
	    
	    String READSET = ReadsetFileLine[1]

		call picard_sam_to_fastq {
			
			input: 
			READSET = READSET,
			IN_BAM = ReadsetFileLine[12],
			TMPDIR = TMPDIR,
			MOD_JAVA = MOD_JAVA,
			MOD_PICARD = MOD_PICARD,
			PICARD_HOME = PICARD_HOME

		}

		call trimmomatic {
			input: 
			READSET = READSET,
			IN_FQ1 = picard_sam_to_fastq.OUT_FQ1, ## or ReadsetFileLine[10]
			IN_FQ2 = picard_sam_to_fastq.OUT_FQ2, ## or ReadsetFileLine[11]
			ADAPTER1 = ReadsetFileLine[6],
			ADAPTER2 = ReadsetFileLine[7],
			MOD_JAVA = MOD_JAVA,
			MOD_TRIMMOMATIC = MOD_TRIMMOMATIC,
			TRIMMOMATIC_JAR = TRIMMOMATIC_JAR

		}

		call bwa_mem_picard_sort_sam {
			input:
			READSET = READSET,
			IN_FQ1 = trimmomatic.OUT_FQ1_TRIM,
			IN_FQ2 = trimmomatic.OUT_FQ2_TRIM,
			MOD_JAVA = MOD_JAVA,
			MOD_BWA = MOD_BWA,
			MOD_PICARD = MOD_PICARD,
			PICARD_HOME = PICARD_HOME,
			TMPDIR = TMPDIR
		}

		## IF 1 readset, symlink else merge
		call sambamba_merge_sam_files {
			input:

		}

	}
	
}






#-------------------------------------------------------------------------------
# Tasks
#-------------------------------------------------------------------------------



task picard_sam_to_fastq {

	String READSET
	String IN_BAM
	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	String MOD_JAVA
	String MOD_PICARD
	String PICARD_HOME

	command <<<

module purge && \
module load ${MOD_JAVA} ${MOD_PICARD} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${PICARD_HOME}/picard.jar SamToFastq \
VALIDATION_STRINGENCY=LENIENT \
INPUT=${IN_BAM} \
FASTQ=${READSET}.pair1.fastq.gz \
SECOND_END_FASTQ=${READSET}.pair2.fastq.gz
	>>>

	output {
	File OUT_FQ1 = "${READSET}.pair1.fastq.gz"
	File OUT_FQ2 = "${READSET}.pair1.fastq.gz"
	}
}


task trimmomatic {

	String READSET
	File IN_FQ1
	File IN_FQ2
	String ADAPTER1
	String ADAPTER2	
	String MOD_JAVA
	String MOD_TRIMMOMATIC
	String TRIMMOMATIC_JAR
	Int THREADS
	Int BUFFER
	String RAM
	Int TRAILING
	Int MINLEN


command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_TRIMMOMATIC} && \
`cat > ${READSET}.trim.adapters.fa << END
>Prefix/1
${ADAPTER1}
>Prefix/2
${ADAPTER2}
END
` && \
java -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${TRIMMOMATIC_JAR} PE \
  -threads ${THREADS} \
  -phred33 \
  ${IN_FQ1} \
  ${IN_FQ2} \
  ${READSET}.trim.pair1.fastq.gz \
  ${READSET}.trim.single1.fastq.gz \
  ${READSET}.trim.pair2.fastq.gz \
  ${READSET}.trim.single2.fastq.gz \
  ILLUMINACLIP:${READSET}.adapters.fa:2:30:15 \
  TRAILING:${TRAILING} \
  MINLEN:${MINLEN} \
  2> ${READSET}.trim.log
	>>>

output {
	File OUT_FQ1_TRIM = "${READSET}.trim.pair1.fastq.gz"
	File OUT_FQ1_TRIM_SINGLE = "${READSET}.trim.single1.fastq.gz"
	File OUT_FQ2_TRIM = "${READSET}.trim.pair2.fastq.gz"
	File OUT_FQ2_TRIM_SINGLE = "${READSET}.trim.single2.fastq.gz"
	File OUT_TRIM_LOG = "${READSET}.trim.log"
	
	}
}


task bwa_mem_picard_sort_sam {

	String READSET
	File IN_FQ1
	File IN_FQ2
	File REF_FASTA
	String BAM_HEADER

	String MOD_JAVA
	String MOD_BWA
	String MOD_PICARD

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	String PICARD_HOME
	Int MAX_REC

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_BWA} ${MOD_PICARD} && \
bwa mem  \
  -M -t 15 \
  -R ${BAM_HEADER} \
  ${REF_FASTA} \
  ${IN_FQ1} \
  ${IN_FQ2} | \
 java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${PICARD_HOME}/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=${READSET}.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=${MAX_REC}
	>>>

output {
	File OUT_BAM="${READSET}.sorted.bam"
	File OUT_BAI="${READSET}.sorted.bai"
	
	}
}




task sambamba_merge_sam_files {

    String SAMPLE
    String PREFIX
    Array[File] IN_BAMS

    String MOD_SAMTOOLS
	String MOD_SAMBAMBA

command <<<
module purge && \
module load ${MOD_SAMTOOLS} ${MOD_SAMBAMBA} && \
sambamba merge -t 7 \
${SAMPLE}.sorted.bam \
(${sep=" " IN_BAMS})
	>>>

output {
	File OUT_MERGED_BAM="${SAMPLE}${PREFIX}"
	
	}
}



##gatk_indel_realigner


task fix_mate_by_coordinate {

	String SAMPLE
	File IN_BAM
	
	String MOD_JAVA
	String MOD_BVATOOLS
	String BVATOOLS_JAR
	String MOD_PICARD
	String PICARD_HOME

	Int THREADS
	Int BUFFER1
	Int BUFFER2
	Int MAX_REC
	String TMPDIR

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_BVATOOLS} ${MOD_PICARD} && \
java -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${BVATOOLS_JAR} \
  groupfixmate \
  --level 1 \
  --bam ${IN_BAM} \
  --out ${SAMPLE}.matefixed.sorted.tmp.bam && \
 java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER2} -Xmx${RAM} -jar ${PICARD_HOME}/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${TMPDIR} \
 INPUT=${SAMPLE}.matefixed.sorted.tmp.bam \
 OUTPUT=${SAMPLE}.matefixed.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=${MAX_REC}
	>>>

output {
	File OUT_TMP_BAM="${SAMPLE}.matefixed.sorted.tmp.bam"
	File OUT_BAM="${SAMPLE}.matefixed.sorted.bam"
	File OUT_BAI="${SAMPLE}.matefixed.sorted.bai"
	
	}
}



task picard_mark_duplicates {

	String SAMPLE
	File IN_BAM

	String MOD_JAVA
	String MOD_PICARD
	String PICARD_HOME
	String TMPDIR

	Int THREADS
	Int BUFFER
	String RAM
	Int MAX_REC


	
command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_PICARD} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${PICARD_HOME}/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${TMPDIR} \
 INPUT=${IN_BAM} \
 OUTPUT=${SAMPLE}.sorted.dup.bam \
 METRICS_FILE=${SAMPLE}.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=${MAX_REC}
	>>>

output {
	File OUT_BAM="${SAMPLE}.sorted.dup.bam"
	File OUT_BAI="${SAMPLE}.sorted.dup.bai"
	File OUT_METRICS="${SAMPLE}.sorted.dup.metrics"
	}
}


## TO COMPLETE
task recalibration_report {

	String SAMPLE
	File IN_BAM
	File REF_FASTA

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM



command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_GATK} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type BaseRecalibrator --bqsrBAQGapOpenPenalty 30 \
  -nt 1 --num_cpu_threads_per_data_thread 12 \
  --input_file ${IN_BAM} \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa  \
  --knownSites /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.dbSNP142.vcf.gz \
  --knownSites /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz \
  --knownSites /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/S1.sorted.dup.recalibration_report.grp
	>>>

output {
	String OUT1="alignment/S1/S1.sorted.dup.recalibration_report.grp"
	
	}
}


task recalibration {

	String IN1="alignment/S1/S1.sorted.dup.bam"
	String IN2="alignment/S1/S1.sorted.dup.recalibration_report.grp"
	

command <<<
	module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 mugqic/samtools/1.4.1 mugqic/sambamba/0.6.6 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=4194304 -Xmx24G -jar $GATK_JAR \
  --analysis_type PrintReads --generate_md5 \
  -nt 1 --num_cpu_threads_per_data_thread 5 \
  --input_file alignment/S1/S1.sorted.dup.bam \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --BQSR alignment/S1/S1.sorted.dup.recalibration_report.grp \
  --out alignment/S1/S1.sorted.dup.recal.bam && \
sambamba index -t 10 \
  alignment/S1/S1.sorted.dup.recal.bam \
  alignment/S1/S1.sorted.dup.recal.bam.bai
	>>>

output {
	String OUT1="alignment/S1/S1.sorted.dup.recal.bam"
	String OUT2="alignment/S1/S1.sorted.dup.recal.bai"
	String OUT3="alignment/S1/S1.sorted.dup.recal.bam.bai"
	
	}
}

#####


task metrics_dna_picard_metrics_main {

	String SAMPLE
	File IN_BAM

	File GENOME_FASTA

	String MOD_JAVA
	String MOD_PICARD
	String PICARD_HOME
	String MOD_R

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int MAX_REC

	

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_PICARD} ${MOD_R} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${PICARD_HOME}/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${TMPDIR} \
 REFERENCE_SEQUENCE=${GENOME_FASTA} \
 INPUT=${IN_BAM} \
 OUTPUT=${SAMPLE}.picard_metrics.all.metrics \
 MAX_RECORDS_IN_RAM=${MAX_REC}
	>>>

output {

	File OUT_ALIGN="${SAMPLE}.picard_metrics.all.metrics.alignment_summary_metrics"
	File OUT_INSERT="${SAMPLE}.picard_metrics.all.metrics.insert_size_metrics"
	File OUT_INSERT_PDF="${SAMPLE}.picard_metrics.all.metrics.insert_size_histogram.pdf"
	File OUT_CYCLE="${SAMPLE}.picard_metrics.all.metrics.quality_by_cycle_metrics"
	File OUT_CYCLE_PDF="${SAMPLE}.picard_metrics.all.metrics.quality_by_cycle.pdf"
	File OUT_DIST="${SAMPLE}.picard_metrics.all.metrics.quality_distribution_metrics"
	File OUT_DIST_PDF="${SAMPLE}.picard_metrics.all.metrics.quality_distribution.pdf"
	
	}
}



task metrics_dna_picard_metrics_oxog {

	String SAMPLE
	File IN_BAM
	
	File DB_SNP
	File GENOME_FASTA

	String MOD_JAVA
	String MOD_PICARD
	String PICARD_HOME
	String MOD_R

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int MAX_REC


command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_PICARD} ${MOD_R} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx6G -jar ${PICARD_HOME}/picard.jar CollectOxoGMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=${TMPDIR} \
 INPUT=${IN_BAM} \
 OUTPUT=${SAMPLE}.picard_metrics.oxog_metrics.txt \
 DB_SNP=${DB_SNP} \
 REFERENCE_SEQUENCE=${GENOME_FASTA} \
 MAX_RECORDS_IN_RAM=4000000
	>>>

output {

	File OUT_OXOG="${SAMPLE}.picard_metrics.oxog_metrics.txt"
	
	}
}



task metrics_dna_picard_metrics_biasQc {

	String SAMPLE
	File IN_BAM

	File GENOME_FASTA

	String MOD_JAVA
	String MOD_PICARD
	String PICARD_HOME
	String MOD_R

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int MAX_REC

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_PICARD} ${MOD_R} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${PICARD_HOME}/picard.jar CollectGcBiasMetrics \
 VALIDATION_STRINGENCY=SILENT ALSO_IGNORE_DUPLICATES=TRUE \
 TMP_DIR=${TMPDIR} \
 INPUT=${IN_BAM} \
 OUTPUT=${SAMPLE}.picard_metrics.qcbias_metrics.txt \
 CHART=${SAMPLE}.picard_metrics.qcbias_metrics.pdf \
 SUMMARY_OUTPUT=${SAMPLE}.picard_metrics.qcbias_summary_metrics.txt \
 REFERENCE_SEQUENCE=${GENOME_FASTA} \
 MAX_RECORDS_IN_RAM=${MAX_REC}
	>>>

output {

	File OUT_BIAS="${SAMPLE}.picard_metrics.qcbias_metrics.txt"
	
	}
}


task metrics_dna_sample_qualimap {

	String SAMPLE
	File IN_BAM

	String SPECIES

	String MOD_JAVA
	String MOD_QUALIMAP

	String RAM


command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_QUALIMAP} && \
qualimap bamqc -nt 11 -gd ${SPECIES} \
  -bam ${IN_BAM} -outdir ${SAMPLE} \
  --java-mem-size=${RAM}
	>>>

output {

	File OUT_QUALIMAP="${SAMPLE}.qualimap.genome_results.txt"
	
	}
}


task metrics_dna_sambamba_flagstat {

	String SAMPLE
	File IN_BAM
	
	String MOD_SAMBAMBA

command <<<
module purge && \
module load ${MOD_SAMBAMBA} && \
sambamba flagstat -t 5 \
  ${IN_BAM} \
  > ${SAMPLE}.flagstat
	>>>

output {
	File OUT_FLAGSTAT="${SAMPLE}.flagstat"
	
	}
}


task metrics_dna_fastqc {

	String SAMPLE
	File IN_BAM
	
	String ADAPTER1
	String ADAPTER2

	String MOD_FASTQC
	String MOD_JAVA

command <<<
module purge && \
module load ${MOD_FASTQC} ${MOD_JAVA} && \
`cat > adapter.tsv << END
>Adapter1	${ADAPTER1}
>Adapter2	${ADAPTER2}
END` && \
fastqc \
  -o ${SAMPLE}.fastqc \
  -t 3 \
  -a adapter.tsv \
  -f bam \
  ${IN_BAM}
	>>>

output {
	File OUT_FASTQC="${SAMPLE}.sorted.dup_fastqc.zip"
	
	}
}


task gatk_callable_loci {

	String SAMPLE
	File IN_BAM

	File GENOME_FASTA
	
	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int MINDEPTH
	Int MAXDEPTH


command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_GATK} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type CallableLoci -dt none --minDepth ${MINDEPTH} --maxDepth ${MAXDEPTH} --minDepthForLowMAPQ ${MINDEPTH} --minMappingQuality ${MINDEPTH} --minBaseQuality ${MINMQ} \
  --input_file ${IN_BAM} \
  --reference_sequence ${GENOME_FASTA} \
  --summary ${SAMPLE}.callable.summary.txt \
  --out ${SAMPLE}.callable.bed
	>>>

output {

	File OUT_BED="${SAMPLE}.callable.bed"
	File OUT_SUMMARY="${SAMPLE}.callable.summary.txt"
	
	}
}


task extract_common_snp_freq {

	String SAMPLE
	File IN_BAM

	File DB_SNP_COMMON

	String MOD_JAVA
	String MOD_BVATOOLS
	String BVATOOLS_JAR
	
	Int THREADS
	Int BUFFER
	String RAM


command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_BVATOOLS} && \
java -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${BVATOOLS_JAR} \
  basefreq \
  --pos ${DB_SNP_COMMON} \
  --bam ${IN_BAM} \
  --out ${SAMPLE}.commonSNPs.alleleFreq.csv
	>>>

output {

	File OUT_COMMONSNPS="${SAMPLE}.commonSNPs.alleleFreq.csv"
	
	}
}



task baf_plot {

	String SAMPLE
	File IN_COMMONSNPS
	
	File DB_SNP_COMMON
	File GENOME_DICT
	String CHR_EXCLUDE	

	String MOD_JAVA
	String MOD_BVATOOLS
	String BVATOOLS_JAR
	
	Int THREADS
	Int BUFFER
	String RAM

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_BVATOOLS} && \
java -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${BVATOOLS_JAR} \
  ratiobaf --plot --maxDepth ${MAXDEPTH}  --exclude ${CHR_EXCLUDE} \
  --refdict ${GENOME_DICT} \
  --snppos ${DB_SNP_COMMON} \
  --basefreq ${IN_COMMONSNPS} \
  --prefix ${SAMPLE}.ratioBAF
	>>>

output {

	File OUT_BAF="${SAMPLE}.ratioBAF.png"
	
	}
}


## TO DO
#call gatk_haplotype_caller	#scatters on sample and Chr
#call merge_and_call_individual_gvcf_merge	## gathers on Chr to get samples
#call merge_and_call_individual_gvcf_calls	
#call combine_gvcf_1_JOB_ID
#call combine_gvcf_2_JOB_ID	
#call combine_gvcf_3_JOB_ID	
#call combine_gvcf_4_JOB_ID	

task merge_and_call_combined_gvcf_merges {

	Array[File] VCFS

	File GENOME_FASTA

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM	

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_GATK} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -cp ${GATK_JAR} \
  org.broadinstitute.gatk.tools.CatVariants  \
  --reference ${GENOME_FASTA} \
  (${sep="--variant " VCFS})
  --outputFile allSamples.hc.g.vcf.gz
	>>>

output {

	File OUT_VCF_G="allSamples.hc.g.vcf.gz"
	
	}
}



task merge_and_call_combined_gvcf_calls {

	File IN_VCF_G

	File GENOME_FASTA

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM	

	

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_GATK} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type GenotypeGVCFs --useNewAFCalculator -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence ${GENOME_FASTA} \
  --variant ${IN_VCF_G} \
  --out allSamples.hc.vcf.gz
	>>>

output {
	File OUT_VCF_ALL="allSamples.hc.vcf.gz"
	
	}
}



task variant_recalibrator_prep {

	String IN1="variants/allSamples.hc.vcf.gz"
	

command <<<
	module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 mugqic/R_Bioconductor/3.4.1_3.5 && \
mkdir -p variants && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx24G -jar $GATK_JAR \
  --analysis_type VariantRecalibrator -nt 11 \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  -input variants/allSamples.hc.vcf.gz \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/annotations/hapmap_3.3.b37.sites.vcf.gz -resource:omni,known=false,training=true,truth=false,prior=12.0 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/annotations/1000G_omni2.5.b37.vcf.gz -resource:1000G,known=false,training=true,truth=false,prior=10.0 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/annotations/1000G_phase1.snps.high_confidence.b37.vcf.gz -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.dbSNP142.vcf.gz -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an DP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -mode SNP \
  --recal_file variants/allSamples.hc.snps.recal \
  --tranches_file variants/allSamples.hc.snps.tranches  \
  --rscript_file variants/allSamples.hc.snps.R && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx24G -jar $GATK_JAR \
  --analysis_type VariantRecalibrator -nt 11 \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  -input variants/allSamples.hc.vcf.gz \
  -resource:mills,known=false,training=true,truth=true,prior=12.0 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.dbSNP142.vcf.gz -an QD -an DP -an FS -an ReadPosRankSum -an MQRankSum -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -mode INDEL \
  --recal_file variants/allSamples.hc.indels.recal \
  --tranches_file variants/allSamples.hc.indels.tranches  \
  --rscript_file variants/allSamples.hc.indels.R
	>>>

output {
	String OUT1="variants/allSamples.hc.snps.recal"
	String OUT2="variants/allSamples.hc.snps.tranches"
	String OUT3="variants/allSamples.hc.snps.R"
	String OUT4="variants/allSamples.hc.indels.recal"
	String OUT5="variants/allSamples.hc.indels.tranches"
	String OUT6="variants/allSamples.hc.indels.R"
	
	}
}

task variant_recalibrator_exec {

	String IN1="variants/allSamples.hc.vcf.gz"
	String IN2="variants/allSamples.hc.snps.recal"
	String IN3="variants/allSamples.hc.snps.tranches"
	String IN4="variants/allSamples.hc.indels.recal"
	String IN5="variants/allSamples.hc.indels.tranches"
	

command <<<
	module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p variants && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx24G -jar $GATK_JAR \
  --analysis_type ApplyRecalibration -nt 11 \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  -input variants/allSamples.hc.vcf.gz \
  --ts_filter_level 99.5 -mode SNP \
  --tranches_file variants/allSamples.hc.snps.tranches \
  --recal_file variants/allSamples.hc.snps.recal \
  --out variants/allSamples.hc.snps_raw_indels.vqsr.vcf.gz && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx24G -jar $GATK_JAR \
  --analysis_type ApplyRecalibration -nt 11 \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  -input variants/allSamples.hc.snps_raw_indels.vqsr.vcf.gz \
  --ts_filter_level 99.0 -mode INDEL \
  --tranches_file variants/allSamples.hc.indels.tranches \
  --recal_file variants/allSamples.hc.indels.recal \
  --out variants/allSamples.hc.vqsr.vcf
	>>>

output {
	String OUT1="variants/allSamples.hc.snps_raw_indels.vqsr.vcf.gz"
	String OUT2="variants/allSamples.hc.vqsr.vcf"
	
	}
}

##


task haplotype_caller_decompose_and_normalize {

	File IN_VQSR
	File GENOME_FASTA
	

command <<<
module purge && \
module load ${MOD_HTSLIB} ${MOD_VT} && \
zless variants/${IN_VQSR} | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r ${GENOME_FASTA} -  \
          | \
bgzip -cf \
 > \
allSamples.hc.vqsr.vt.vcf.gz && tabix -pvcf allSamples.hc.vqsr.vt.vcf.gz   \

	>>>

output {

	String OUT_VT="allSamples.hc.vqsr.vt.vcf.gz"
	
	}
}



task haplotype_caller_flag_mappability {

	File IN_VT

	File IlluEXCLUSION
	
	String MOD_VCFTOOLS
	String MOD_TABIX
	String MOD_HTSLIB

command <<<
module purge && \
module load ${MOD_VCFTOOLS} ${MOD_TABIX} ${MOD_HTSLIB} && \
vcf-annotate \
  -d key=INFO,ID=MIL,Number=1,Type=String,Description='Mappability annotation. 300IS 40SD 1SHI. HC = too high coverage (>400), LC = too low coverage (<50), MQ = too low mean mapQ (<20), ND = no data at the position' \
  -c CHROM,FROM,TO,INFO/MIL \
  -a ${IlluEXCLUSION} \
  ${IN_VT} | \
bgzip -cf \
 > \
allSamples.hc.vqsr.vt.mil.vcf.gz && tabix -pvcf allSamples.hc.vqsr.vt.mil.vcf.gz \
        
	>>>

output {

	File OUT_MIL="allSamples.hc.vqsr.vt.mil.vcf.gz"
	
	}
}



task haplotype_caller_snp_id_annotation {

	File IN_MIL

	File DB_SNP

	String MOD_JAVA
	String MOD_SNPEFF
	String SNPEFF_HOME
	String MOD_HTSLIB
	String MOD_TABIX

	String TMPDIR
	Int THREADS
	String RAM
	

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_SNPEFF} ${MOD_HTSLIB} ${MOD_TABIX} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Xmx${RAM} -jar ${SNPEFF_HOME}/SnpSift.jar annotate \
  ${DB_SNP} \
  ${IN_MIL} | \
bgzip -cf \
 > \
allSamples.hc.vqsr.vt.mil.snpId.vcf.gz && tabix -pvcf allSamples.hc.vqsr.vt.mil.snpId.vcf.gz \
        
	>>>

output {

	File OUT_SNPID="allSamples.hc.vqsr.vt.mil.snpId.vcf.gz"
	
	}
}



task haplotype_caller_snp_effect {

	File IN_SNPID

	String ASSEMBLY
	
	String MOD_JAVA
	String MOD_SNPEFF
	String SNPEFF_HOME
	String MOD_HTSLIB
	String MOD_TABIX

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM


command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_SNPEFF} ${MOD_HTSLIB} ${MOD_TABIX} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Xmx${RAM} -jar ${SNPEFF_HOME}/snpEff.jar eff -lof \
   \
  -c ${SNPEFF_HOME}/snpEff.config \
  -i vcf \
  -o vcf \
  -csvStats allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.stats.csv \
  -stats allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.stats.html \
  ${ASSEMBLY} \
  ${IN_SNPID} > allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf && \
bgzip -cf \
 \
 allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf > \
 allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.gz && tabix -pvcf allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.gz \
        
	>>>

output {
	File OUT_SNPEFF="allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf"
	File OUT_STATS="allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.stats.csv"
	File OUT_ZIPPED="allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.gz"
	
	}
}



task haplotype_caller_dbnsfp_annotation {

	File IN_ZIPPED

	File DB_SFP

	String MOD_JAVA
	String MOD_SNPEFF
	String SNPEFF_HOME
	String MOD_HTSLIB
	String MOD_TABIX

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_SNPEFF} ${MOD_HTSLIB} ${MOD_TABIX} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Xmx${RAM} -jar ${SNPEFF_HOME}/SnpSift.jar dbnsfp \
  -v -db ${DB_SFP} \
  ${IN_ZIPPED} \
  > allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf && \
bgzip -cf \
 \
 allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf > \
 allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf.gz && tabix -pvcf allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf.gz \
        
	>>>

output {

	File OUT_SFP="allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf"
	File OUT_ZIPPED="allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf.gz"
	
	}
}



task haplotype_caller_gemini_annotations {

	File IN_ZIPPED

	String MOD_GEMINI
	String MOD_HTSLIB

	String TMPDIR
	Int CPU

command <<<
module purge && \
module load ${MOD_GEMINI} ${MOD_HTSLIB} && \
gemini load -v ${IN_ZIPPED} \
  -t snpEff --cores ${CPU} --save-info-string \
  --tempdir ${TMPDIR} \
  allSamples.gemini.db
	>>>

output {

	File OUT_GEMINI="allSamples.gemini.db"
	
	}
}



task haplotype_caller_metrics_vcf_stats {

	File IN_ZIPPED
	File IN_STATS

	File GENOME_DICT

	String MOD_PYTHON
	String PYTHON_TOOLS
	String MOD_MUGQICTOOLS

command <<<
module purge && \
module load ${MOD_PYTHON} ${MOD_MUGQICTOOLS} && \
python ${PYTHON_TOOLS}/vcfStats.py \
  -v ${IN_ZIPPED} \
  -d ${GENOME_DICT} \
  -o allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf.part_changeRate.tsv \
  -f ${IN_STATS}
	>>>

output {

	String OUT_CHGRATE="allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf.part_changeRate.tsv"
	
	}
}



task run_multiqc {
	
	Array[File] METRICS
	## Need to make it run at the very end	

command <<<
module purge && \
module load ${MOD_MULTIQC} && \
multiqc .
	>>>

output {

	File OUT_MULTIQC="multiqc_report"
	
	}
}




task cram_output {

	String SAMPLE
	File IN_BAM
	
	File GENOME_FASTA

	String MOD_SAMTOOLS

command <<<
module purge && \
module load ${MOD_SAMTOOLS} && \
samtools view -h -T ${GENOME_FASTA} -C \
  ${IN_BAM} \
  > ${SAMPLE}.sorted.dup.recal.cram
	>>>

output {

	File OUT_CRAM="${SAMPLE}.sorted.dup.recal.cram"
	
	}
}











