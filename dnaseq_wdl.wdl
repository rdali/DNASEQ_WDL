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

	## ENV VARS
	String OUT_DIR
	String TMPDIR

	## MODULE VARS
	String PICARD_HOME
	String MOD_JAVA
	String MOD_PICARD
	String MOD_TRIMMOMATIC
	String TRIMMOMATIC_JAR
	String MOD_BWA

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
	String OUT_BAM="${READSET}.sorted.bam"
	String OUT_BAI="${READSET}.sorted.bai"
	
	}
}




















