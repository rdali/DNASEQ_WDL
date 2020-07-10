#-------------------------------------------------------------------------------
# DnaSeq WDL workflow
# Version1 based on GenPipes 3.1.5-beta
# Created on: 2020-06-18
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# DnaSeq Tasks
#-------------------------------------------------------------------------------



task picard_sam_to_fastq {

	String READSET
	File IN_BAM

	String MOD_JAVA
	String MOD_PICARD
	String PICARD_HOME

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM


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
	String BWA_INDX  ## must be String so that bwa can find the index
	String BAM_HEADER

	String MOD_JAVA
	String MOD_BWA
	String MOD_PICARD
	String PICARD_HOME

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int MAX_REC

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_BWA} ${MOD_PICARD} && \
bwa mem \
  -M -t 15 \
  -R "${BAM_HEADER}" \
  ${BWA_INDX} \
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
${sep=" " IN_BAMS}
	>>>

output {

	File OUT_MERGED_BAM="${SAMPLE}${PREFIX}"

	}
}



task gatk_indel_realigner {

	String SAMPLE
	File IN_BAM
	Array[String]? INTERVALS
	Array[String]? EXCLUDE
	String INTERVAL_NAME

	String GENOME_FASTA
	File G1000

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int NT
	Int NCT
	Int MAX_REC


command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_GATK} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type RealignerTargetCreator -nct ${NT} -nt ${NCT} \
  --reference_sequence ${GENOME_FASTA} \
  --input_file ${IN_BAM} \
  --known ${G1000} \
  --out ${SAMPLE}.sorted.realigned.${INTERVAL_NAME}.intervals \
  ${true=' --intervals ' false='' defined(INTERVALS)}${sep=' --intervals ' INTERVALS} \
  ${true=' --excludeIntervals ' false='' defined(EXCLUDE)}${sep=' --excludeIntervals ' EXCLUDE}  && \
	java -Djava.io.tmpdir=${TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type IndelRealigner -nt ${NT} -nct ${NT} \
  --reference_sequence ${GENOME_FASTA} \
  --input_file ${IN_BAM} \
  --targetIntervals ${SAMPLE}.sorted.realigned.${INTERVAL_NAME}.intervals \
  --knownAlleles ${G1000} \
  --out ${SAMPLE}.sorted.realigned.${INTERVAL_NAME}.bam \
  ${true=' --intervals ' false='' defined(INTERVALS)}${sep=' --intervals ' INTERVALS} \
  ${true=' --excludeIntervals ' false='' defined(EXCLUDE)}${sep=' --excludeIntervals ' EXCLUDE} \
  --maxReadsInMemory ${MAX_REC}
	>>>

output {

	File OUT_INTRVL="${SAMPLE}.sorted.realigned.${INTERVAL_NAME}.intervals"
	File OUT_BAM="${SAMPLE}.sorted.realigned.${INTERVAL_NAME}.bam"

	}
}



task fix_mate_by_coordinate {

	String SAMPLE
	File IN_BAM

	String MOD_JAVA
	String MOD_BVATOOLS
	String BVATOOLS_JAR
	String MOD_PICARD
	String PICARD_HOME

	String TMPDIR
	Int THREADS
	Int BUFFER1
	Int BUFFER2
	String RAM
	Int MAX_REC

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_BVATOOLS} ${MOD_PICARD} && \
java -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER1} -Xmx${RAM} -jar ${BVATOOLS_JAR} \
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


task recalibration_report {

	String SAMPLE
	File IN_BAM
	String GENOME_FASTA
	Array[File] KNOWNSITES

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int BQSR
	Int CPUTHREADS


##  KNOWNSITES = [DB_SNP, G1000, GNOMAD]

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_GATK} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type BaseRecalibrator --bqsrBAQGapOpenPenalty ${BQSR} \
  -nt 1 --num_cpu_threads_per_data_thread ${CPUTHREADS} \
  --input_file ${IN_BAM} \
  --reference_sequence ${GENOME_FASTA} \
  ${sep=" --knownSites " KNOWNSITES} \
  --out ${SAMPLE}.sorted.dup.recalibration_report.grp
	>>>

output {

	File OUT_CLB_RPT="${SAMPLE}.sorted.dup.recalibration_report.grp"

	}
}


task recalibration {

	String SAMPLE
	File IN_BAM

	String GENOME_FASTA
	File IN_CLB_RPT

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR
	String MOD_SAMTOOLS
	String MOD_SAMBAMBA

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int CPUTHREADS

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_GATK} ${MOD_SAMTOOLS} ${MOD_SAMBAMBA} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type PrintReads --generate_md5 \
  -nt 1 --num_cpu_threads_per_data_thread ${CPUTHREADS} \
  --input_file ${IN_BAM} \
  --reference_sequence ${GENOME_FASTA} \
  --BQSR ${IN_CLB_RPT} \
  --out ${SAMPLE}.sorted.dup.recal.bam && \
sambamba index -t 10 \
  ${SAMPLE}.sorted.dup.recal.bam \
  ${SAMPLE}.sorted.dup.recal.bam.bai
	>>>

output {

	File OUT_BAM="${SAMPLE}.sorted.dup.recal.bam"
	File OUT_BAI="${SAMPLE}.sorted.dup.recal.bai"
	File OUT_BAMBAI="${SAMPLE}.sorted.dup.recal.bam.bai"

	}
}


task metrics_dna_picard_metrics_main {

	String SAMPLE
	File IN_BAM

	String GENOME_FASTA

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
	String GENOME_FASTA

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
 VALIDATION_STRINGENCY=SILENT \
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

	String GENOME_FASTA

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

	String GENOME_FASTA

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int MINDEPTH
	Int MAXDEPTH
	Int MINMQ


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
	Int MAXDEPTH

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


task gatk_haplotype_caller {

	String SAMPLE
	File IN_BAM
	Array[String] CHR_EXCLUDE

	String GENOME_FASTA
	String INTERVAL_NAME
	Array[String] INTERVALS

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int NCT

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_GATK} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct ${NCT} -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence ${GENOME_FASTA} \
  --input_file ${IN_BAM} \
  --out ${SAMPLE}.${INTERVAL_NAME}.hc.g.vcf.gz \
  ${sep=" --intervals " INTERVALS} \
  ${sep=" --excludeIntervals " CHR_EXCLUDE}
	>>>

output {

	File OUT_VCF_INTRVL="${SAMPLE}.${INTERVAL_NAME}.hc.g.vcf.gz"
	File OUT_TBI_INTRVL="${SAMPLE}.${INTERVAL_NAME}.hc.g.vcf.gz.tbi"

	}
}



task merge_and_call_individual_gvcf_merge {

	String SAMPLE
	Array[File] IN_VCFS_INTRVL

	String GENOME_FASTA

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
  ${sep=" --variant " IN_VCFS_INTRVL} \
  --outputFile ${SAMPLE}.hc.g.vcf.gz
	>>>

output {

	File OUT_VCF_G="${SAMPLE}.hc.g.vcf.gz"

	}
}


task merge_and_call_individual_gvcf_calls {

	String SAMPLE
	File IN_VCF_G

	String GENOME_FASTA
	String INTERVAL

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int NCT

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_GATK} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type GenotypeGVCFs --useNewAFCalculator -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence ${GENOME_FASTA} \
  --variant ${IN_VCF_G} \
  --out ${SAMPLE}.hc.vcf.gz
	>>>

output {

	File OUT_VCF="${SAMPLE}.hc.vcf.gz"

	}
}

task combine_gvcf {

	Array[File] IN_VCFS_G
	Array[String] CHR_EXCLUDE

	String GENOME_FASTA
	String INTERVAL

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
java -Djava.io.tmpdir=${TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type CombineGVCFs  \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence ${GENOME_FASTA} \
	${sep=" --variant " IN_VCFS_G} \
  --out allSamples.${INTERVAL}.hc.g.vcf.bgz \
	${sep=" --excludeIntervals " CHR_EXCLUDE}
	>>>

output {

	File OUT_VCF_G="allSamples.${INTERVAL}.hc.g.vcf.bgz"

	}
}



task merge_and_call_combined_gvcf_merges {

	Array[File] VCFS

	String GENOME_FASTA

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
  org.broadinstitute.gatk.tools.CatVariants \
  --reference ${GENOME_FASTA} \
  ${sep=" --variant " VCFS} \
  --outputFile allSamples.hc.g.vcf.gz
	>>>

output {

	File OUT_VCF_G="allSamples.hc.g.vcf.gz"

	}
}



task merge_and_call_combined_gvcf_calls {

	File IN_VCF_G

	String GENOME_FASTA

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

	File IN_VCF_ALL

	String GENOME_FASTA

	String RSRC_SNP
	String RSCR_INDL

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR
	String MOD_R

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int NT


command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_GATK} ${MOD_R} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type VariantRecalibrator -nt ${NT} \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence ${GENOME_FASTA} \
  -input ${IN_VCF_ALL} \
  ${RSRC_SNP} \
  --recal_file allSamples.hc.snps.recal \
  --tranches_file allSamples.hc.snps.tranches \
  --rscript_file allSamples.hc.snps.R && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type VariantRecalibrator -nt ${NT} \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence ${GENOME_FASTA} \
  -input ${IN_VCF_ALL} \
  ${RSRC_SNP} \
  --recal_file allSamples.hc.indels.recal \
  --tranches_file allSamples.hc.indels.tranches \
  --rscript_file allSamples.hc.indels.R
	>>>

output {

	File OUT_SNP_RCL="allSamples.hc.snps.recal"
	File OUT_SNP_TRNCH="allSamples.hc.snps.tranches"
	File OUT_SNP_R="allSamples.hc.snps.R"
	File OUT_INDL_RCL="allSamples.hc.indels.recal"
	File OUT_INDL_TRNCH="allSamples.hc.indels.tranches"
	File OUT_INDL_R="allSamples.hc.indels.R"

	}
}

task variant_recalibrator_exec {

	File IN_VCF_ALL
	File IN_SNP_RCL
	File IN_SNP_TRNCH
	File IN_SNP_R
	File IN_INDL_RCL
	File IN_INDL_TRNCH
	File IN_INDL_R

	String GENOME_FASTA

	String MOD_JAVA
	String MOD_GATK
	String GATK_JAR

	String TMPDIR
	Int THREADS
	Int BUFFER
	String RAM
	Int NT
	Float FILT_SNP
	Float FILT_INDL

command <<<
module purge && \
module load ${MOD_JAVA} ${MOD_GATK} && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type ApplyRecalibration -nt ${NT} \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence ${GENOME_FASTA} \
  -input ${IN_VCF_ALL} \
  --ts_filter_level ${FILT_SNP} -mode SNP \
  --tranches_file ${IN_SNP_TRNCH} \
  --recal_file ${IN_SNP_RCL} \
  --out allSamples.hc.snps_raw_indels.vqsr.vcf.gz && \
java -Djava.io.tmpdir=${TMPDIR} -XX:ParallelGCThreads=${THREADS} -Dsamjdk.buffer_size=${BUFFER} -Xmx${RAM} -jar ${GATK_JAR} \
  --analysis_type ApplyRecalibration -nt ${NT} \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence {$GENOME_FASTA} \
  -input allSamples.hc.snps_raw_indels.vqsr.vcf.gz \
  --ts_filter_level ${FILT_INDL} -mode INDEL \
  --tranches_file ${IN_INDL_TRNCH} \
  --recal_file ${IN_INDL_RCL} \
  --out allSamples.hc.vqsr.vcf
	>>>

output {

	File OUT_VQSR_RAW="allSamples.hc.snps_raw_indels.vqsr.vcf.gz"
	File OUT_VQSR="allSamples.hc.vqsr.vcf"

	}
}

##


task haplotype_caller_decompose_and_normalize {

	File IN_VQSR
	String GENOME_FASTA

	String MOD_HTSLIB
	String MOD_VT


command <<<
module purge && \
module load ${MOD_HTSLIB} ${MOD_VT} && \
zless variants/${IN_VQSR} | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r ${GENOME_FASTA} - \
		  | \
bgzip -cf \
 > \
allSamples.hc.vqsr.vt.vcf.gz && tabix -pvcf allSamples.hc.vqsr.vt.vcf.gz \

	>>>

output {

	File OUT_VT="allSamples.hc.vqsr.vt.vcf.gz"

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

	File OUT_CHGRATE="allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf.part_changeRate.tsv"

	}
}



task run_multiqc {

	Array[File] METRICS
	## Need to make it run at the very end
	String MOD_MULTIQC

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

	String GENOME_FASTA

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



#-------------------------------------------------------------------------------
# General Tasks
#-------------------------------------------------------------------------------



task concat_arrays {

    Array[File] A
    Array[File] B

    command {

        cat write_lines(A)
        cat write_lines(B)

    }

    output {

        Array[File] OUT = read_lines(stdout())

    }
}
