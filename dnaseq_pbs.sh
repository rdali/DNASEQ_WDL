#!/bin/bash
# Exit immediately on error

set -eu -o pipefail

#-------------------------------------------------------------------------------
# DnaSeq PBS/TORQUE Job Submission Bash script
# Version: 3.1.5-beta
# Created on: 2020-06-23T14:56:21
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   sym_link_fastq: 1 job
#   trimmomatic: 1 job
#   merge_trimmomatic_stats: 1 job
#   skewer_trimming: 1 job
#   bwa_mem_picard_sort_sam: 1 job
#   sambamba_merge_sam_files: 1 job
#   gatk_indel_realigner: 23 jobs
#   sambamba_merge_realigned: 1 job
#   fix_mate_by_coordinate: 1 job
#   picard_mark_duplicates: 1 job
#   recalibration: 2 jobs
#   sym_link_final_bam: 1 job
#   metrics_dna_picard_metrics: 3 jobs
#   metrics_dna_sample_qualimap: 1 job
#   metrics_dna_sambamba_flagstat: 1 job
#   metrics_dna_fastqc: 1 job
#   picard_calculate_hs_metrics: 0 job... skipping
#   gatk_callable_loci: 1 job
#   extract_common_snp_freq: 1 job
#   baf_plot: 1 job
#   gatk_haplotype_caller: 23 jobs
#   merge_and_call_individual_gvcf: 2 jobs
#   combine_gvcf: 4 jobs
#   merge_and_call_combined_gvcf: 2 jobs
#   variant_recalibrator: 2 jobs
#   haplotype_caller_decompose_and_normalize: 1 job
#   haplotype_caller_flag_mappability: 1 job
#   haplotype_caller_snp_id_annotation: 1 job
#   haplotype_caller_snp_effect: 1 job
#   haplotype_caller_dbnsfp_annotation: 1 job
#   haplotype_caller_gemini_annotations: 1 job
#   haplotype_caller_metrics_vcf_stats: 1 job
#   run_multiqc: 1 job
#   cram_output: 1 job
#   TOTAL: 86 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/DnaSeq_job_list_$TIMESTAMP
export CONFIG_FILES="/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/repos/genpipes//pipelines/dnaseq/dnaseq.base.ini,/home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/repos/genpipes//pipelines/dnaseq/dnaseq.beluga.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: sym_link_fastq
#-------------------------------------------------------------------------------
STEP=sym_link_fastq

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: sym_link_fastq_1_JOB_ID: sym_link_fastq.pair_end.S1_R1
#-------------------------------------------------------------------------------
JOB_NAME=sym_link_fastq.pair_end.S1_R1
JOB_DEPENDENCIES=
JOB_DONE=job_output/sym_link_fastq/sym_link_fastq.pair_end.S1_R1.6e3c72ee976c933bb12c8df0d6fb23a6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sym_link_fastq.pair_end.S1_R1.6e3c72ee976c933bb12c8df0d6fb23a6.mugqic.done' > $COMMAND
mkdir -p /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/deliverables/S1/wgs/raw_reads && \       
ln -sf \
  /home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/rawReads/NA24385_WE5_R1.fastq.gz \
  /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/deliverables/S1/wgs/raw_reads/NA24385_WE5_R1.fastq.gz && \
mkdir -p /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/deliverables/S1/wgs/raw_reads && \       
ln -sf \
  /home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/rawReads/NA24385_WE5_R2.fastq.gz \
  /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/deliverables/S1/wgs/raw_reads/NA24385_WE5_R2.fastq.gz
sym_link_fastq.pair_end.S1_R1.6e3c72ee976c933bb12c8df0d6fb23a6.mugqic.done
chmod 755 $COMMAND

sym_link_fastq_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=3:00:00 --mem-per-cpu=4000M -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$sym_link_fastq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.S1_R1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.S1_R1
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.S1_R1.3d457a2421e33e42d17a0f4fc4d47e38.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.S1_R1.3d457a2421e33e42d17a0f4fc4d47e38.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/S1 && \
`cat > trim/S1/S1_R1.trim.adapters.fa << END
>Prefix/1
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>Prefix/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
END
` && \
java -XX:ParallelGCThreads=5 -Dsamjdk.buffer_size=1048576 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 5 \
  -phred33 \
  /home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/rawReads/NA24385_WE5_R1.fastq.gz \
  /home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/rawReads/NA24385_WE5_R2.fastq.gz \
  trim/S1/S1_R1.trim.pair1.fastq.gz \
  trim/S1/S1_R1.trim.single1.fastq.gz \
  trim/S1/S1_R1.trim.pair2.fastq.gz \
  trim/S1/S1_R1.trim.single2.fastq.gz \
  ILLUMINACLIP:trim/S1/S1_R1.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/S1/S1_R1.trim.log
trimmomatic.S1_R1.3d457a2421e33e42d17a0f4fc4d47e38.mugqic.done
chmod 755 $COMMAND

trimmomatic_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu=4775M -n 5 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.ec4e5a3d1b6210705d8014a42a608328.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'merge_trimmomatic_stats.ec4e5a3d1b6210705d8014a42a608328.mugqic.done' > $COMMAND
module purge && \
module load mugqic/perl/5.22.1 mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Paired Reads #	Surviving Paired Reads #	Surviving Paired Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/S1/S1_R1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/S1	S1_R1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"	" '{OFS="	"; if (NR==1) {if ($2=="Raw Paired Reads #") {paired=1};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$2=$2*2; $3=$3*2}; raw[$1]+=$2; surviving[$1]+=$3}}END{for (sample in raw){print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"} else {print $1, $2, sprintf("%\47d", $3), sprintf("%\47d", $4), sprintf("%.1f", $5)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/repos/genpipes/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/repos/genpipes/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=30 \
  --variable min_length=50 \
  --variable read_type=Paired \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats.ec4e5a3d1b6210705d8014a42a608328.mugqic.done
chmod 755 $COMMAND

merge_trimmomatic_stats_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu=4000M -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: skewer_trimming
#-------------------------------------------------------------------------------
STEP=skewer_trimming

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: skewer_trimming_1_JOB_ID: skewer_trimming.S1_R1
#-------------------------------------------------------------------------------
JOB_NAME=skewer_trimming.S1_R1
JOB_DEPENDENCIES=
JOB_DONE=job_output/skewer_trimming/skewer_trimming.S1_R1.178735a8db51f5402f3cd347c9f6d9ef.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'skewer_trimming.S1_R1.178735a8db51f5402f3cd347c9f6d9ef.mugqic.done' > $COMMAND
module purge && \
module load mugqic/skewer/0.2.2 && \
mkdir -p trim/S1 && \
`cat > trim/S1/adapter.tsv << END
>Adapter1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>Adapter2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
END` && \
$SKEWER_HOME/./skewer --threads 8 --min 25 -q 25 --compress -f sanger \
  -x trim/S1/adapter.tsv \
  /home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/rawReads/NA24385_WE5_R1.fastq.gz \
  /home/rdali/projects/rrg-bourqueg-ad/rdali/C3G/projects/dnaseqTestSet/rawReads/NA24385_WE5_R2.fastq.gz \
  -o trim/S1/S1_R1 && \
ln -s -f \
  S1_R1-trimmed-pair1.fastq.gz \
  trim/S1/S1_R1.trim.pair1.fastq.gz && \
ls trim/S1 && \
ln -s -f \
  S1_R1-trimmed-pair2.fastq.gz \
  trim/S1/S1_R1.trim.pair2.fastq.gz && \
ls trim/S1
skewer_trimming.S1_R1.178735a8db51f5402f3cd347c9f6d9ef.mugqic.done
chmod 755 $COMMAND

skewer_trimming_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:00 --mem=12G -N 1 -n 10 | grep "[0-9]" | cut -d\  -f4)
echo "$skewer_trimming_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: bwa_mem_picard_sort_sam
#-------------------------------------------------------------------------------
STEP=bwa_mem_picard_sort_sam

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.S1_R1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.S1_R1
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$skewer_trimming_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.S1_R1.1923c32720abc3030577b1b7b7def9aa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'bwa_mem_picard_sort_sam.S1_R1.1923c32720abc3030577b1b7b7def9aa.mugqic.done' > $COMMAND
module purge && \
module load mugqic/bwa/0.7.15 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.9.0 && \
mkdir -p alignment/S1/S1_R1 && \
bwa mem  \
  -M -t 15 \
  -R '@RG	ID:S1_R1	SM:S1	LB:6	PU:run1_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/bwa_index/Homo_sapiens.GRCh37.fa \
  trim/S1/S1_R1.trim.pair1.fastq.gz \
  trim/S1/S1_R1.trim.pair2.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx16G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/S1/S1_R1/S1_R1.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.S1_R1.1923c32720abc3030577b1b7b7def9aa.mugqic.done
chmod 755 $COMMAND

bwa_mem_picard_sort_sam_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=96:00:00 --mem=60G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: sambamba_merge_sam_files
#-------------------------------------------------------------------------------
STEP=sambamba_merge_sam_files

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: sambamba_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.S1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.S1
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/sambamba_merge_sam_files/symlink_readset_sample_bam.S1.5899c1222b9ce6fa1bb342bcc2d99610.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'symlink_readset_sample_bam.S1.5899c1222b9ce6fa1bb342bcc2d99610.mugqic.done' > $COMMAND
mkdir -p alignment/S1 && \
ln -s -f S1_R1/S1_R1.sorted.bam alignment/S1/S1.sorted.bam && \
ln -s -f S1_R1/S1_R1.sorted.bai alignment/S1/S1.sorted.bai && sleep 180
symlink_readset_sample_bam.S1.5899c1222b9ce6fa1bb342bcc2d99610.mugqic.done
chmod 755 $COMMAND

sambamba_merge_sam_files_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu=4000M -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: gatk_indel_realigner
#-------------------------------------------------------------------------------
STEP=gatk_indel_realigner

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_1_JOB_ID: gatk_indel_realigner.S1.0
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.0
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.0.0db12ccf89ad68b1a157f96f40c3b859.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.0.0db12ccf89ad68b1a157f96f40c3b859.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.0.intervals \
  --intervals 1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.0.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.0.bam \
  --intervals 1 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.0.0db12ccf89ad68b1a157f96f40c3b859.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_2_JOB_ID: gatk_indel_realigner.S1.1
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.1
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.1.f110f75b83dc4b04e080cbbc0fa0f3e3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.1.f110f75b83dc4b04e080cbbc0fa0f3e3.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.1.intervals \
  --intervals 2 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.1.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.1.bam \
  --intervals 2 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.1.f110f75b83dc4b04e080cbbc0fa0f3e3.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_2_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_3_JOB_ID: gatk_indel_realigner.S1.2
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.2
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.2.1240f101aa9312661622095eea1be88e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.2.1240f101aa9312661622095eea1be88e.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.2.intervals \
  --intervals 3 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.2.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.2.bam \
  --intervals 3 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.2.1240f101aa9312661622095eea1be88e.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_3_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_4_JOB_ID: gatk_indel_realigner.S1.3
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.3
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.3.ceb600d6c10fc2be092098be6d3f6b10.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.3.ceb600d6c10fc2be092098be6d3f6b10.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.3.intervals \
  --intervals 4 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.3.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.3.bam \
  --intervals 4 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.3.ceb600d6c10fc2be092098be6d3f6b10.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_4_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_5_JOB_ID: gatk_indel_realigner.S1.4
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.4
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.4.7d7e3e0548b63833bdd081dde56c2e82.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.4.7d7e3e0548b63833bdd081dde56c2e82.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.4.intervals \
  --intervals 5 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.4.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.4.bam \
  --intervals 5 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.4.7d7e3e0548b63833bdd081dde56c2e82.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_5_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_6_JOB_ID: gatk_indel_realigner.S1.5
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.5
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.5.04440147960896f995989e27fcc30bf3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.5.04440147960896f995989e27fcc30bf3.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.5.intervals \
  --intervals 6 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.5.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.5.bam \
  --intervals 6 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.5.04440147960896f995989e27fcc30bf3.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_6_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_7_JOB_ID: gatk_indel_realigner.S1.6
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.6
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.6.aa90938f867f41cc452c244e7b67f4ee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.6.aa90938f867f41cc452c244e7b67f4ee.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.6.intervals \
  --intervals 7 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.6.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.6.bam \
  --intervals 7 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.6.aa90938f867f41cc452c244e7b67f4ee.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_7_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_8_JOB_ID: gatk_indel_realigner.S1.7
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.7
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.7.b60b925bdf25401d43f37002a61a1d7b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.7.b60b925bdf25401d43f37002a61a1d7b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.7.intervals \
  --intervals 8 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.7.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.7.bam \
  --intervals 8 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.7.b60b925bdf25401d43f37002a61a1d7b.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_8_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_9_JOB_ID: gatk_indel_realigner.S1.8
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.8
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.8.7c0eb37192402e3cbae839072aeafd6d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.8.7c0eb37192402e3cbae839072aeafd6d.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.8.intervals \
  --intervals 9 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.8.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.8.bam \
  --intervals 9 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.8.7c0eb37192402e3cbae839072aeafd6d.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_9_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_10_JOB_ID: gatk_indel_realigner.S1.9
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.9
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.9.93f73b4b652a204c3f456c96d586697b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.9.93f73b4b652a204c3f456c96d586697b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.9.intervals \
  --intervals 10 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.9.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.9.bam \
  --intervals 10 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.9.93f73b4b652a204c3f456c96d586697b.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_10_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_11_JOB_ID: gatk_indel_realigner.S1.10
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.10
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.10.b0fa7566a6c6ac3f79ae2b622f6b8084.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.10.b0fa7566a6c6ac3f79ae2b622f6b8084.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.10.intervals \
  --intervals 11 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.10.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.10.bam \
  --intervals 11 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.10.b0fa7566a6c6ac3f79ae2b622f6b8084.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_11_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_12_JOB_ID: gatk_indel_realigner.S1.11
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.11
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.11.555a0a3ab3d716f4305926b936446968.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.11.555a0a3ab3d716f4305926b936446968.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.11.intervals \
  --intervals 12 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.11.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.11.bam \
  --intervals 12 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.11.555a0a3ab3d716f4305926b936446968.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_12_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_13_JOB_ID: gatk_indel_realigner.S1.12
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.12
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.12.11a84595d4dcefacdff679953104b2b2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.12.11a84595d4dcefacdff679953104b2b2.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.12.intervals \
  --intervals 13 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.12.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.12.bam \
  --intervals 13 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.12.11a84595d4dcefacdff679953104b2b2.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_13_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_14_JOB_ID: gatk_indel_realigner.S1.13
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.13
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.13.ebbdaa713ccf98cba7b1292c3c0db010.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.13.ebbdaa713ccf98cba7b1292c3c0db010.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.13.intervals \
  --intervals 14 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.13.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.13.bam \
  --intervals 14 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.13.ebbdaa713ccf98cba7b1292c3c0db010.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_14_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_15_JOB_ID: gatk_indel_realigner.S1.14
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.14
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.14.41ed77d54477a84dba2e42224f58bf36.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.14.41ed77d54477a84dba2e42224f58bf36.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.14.intervals \
  --intervals 15 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.14.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.14.bam \
  --intervals 15 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.14.41ed77d54477a84dba2e42224f58bf36.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_15_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_16_JOB_ID: gatk_indel_realigner.S1.15
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.15
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.15.6a98c389ad2196ceece77ed69d354069.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.15.6a98c389ad2196ceece77ed69d354069.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.15.intervals \
  --intervals 16 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.15.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.15.bam \
  --intervals 16 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.15.6a98c389ad2196ceece77ed69d354069.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_16_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_17_JOB_ID: gatk_indel_realigner.S1.16
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.16
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.16.508d6fd2c04017c1a3f9205685413f93.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.16.508d6fd2c04017c1a3f9205685413f93.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.16.intervals \
  --intervals 17 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.16.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.16.bam \
  --intervals 17 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.16.508d6fd2c04017c1a3f9205685413f93.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_17_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_18_JOB_ID: gatk_indel_realigner.S1.17
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.17
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.17.828a1aa30739796ae6ede1e60b2463b6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.17.828a1aa30739796ae6ede1e60b2463b6.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.17.intervals \
  --intervals 18 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.17.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.17.bam \
  --intervals 18 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.17.828a1aa30739796ae6ede1e60b2463b6.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_18_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_19_JOB_ID: gatk_indel_realigner.S1.18
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.18
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.18.4a9fb4aaddf4e3ec97f0d220287a3e0f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.18.4a9fb4aaddf4e3ec97f0d220287a3e0f.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.18.intervals \
  --intervals 19 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.18.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.18.bam \
  --intervals 19 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.18.4a9fb4aaddf4e3ec97f0d220287a3e0f.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_19_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_20_JOB_ID: gatk_indel_realigner.S1.19
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.19
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.19.ea04d964e8c5164579a0c6960f6c2ce6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.19.ea04d964e8c5164579a0c6960f6c2ce6.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.19.intervals \
  --intervals 20 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.19.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.19.bam \
  --intervals 20 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.19.ea04d964e8c5164579a0c6960f6c2ce6.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_20_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_21_JOB_ID: gatk_indel_realigner.S1.20
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.20
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.20.13d28fcf96ad656080ef6fd9374f0774.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.20.13d28fcf96ad656080ef6fd9374f0774.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.20.intervals \
  --intervals 21 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.20.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.20.bam \
  --intervals 21 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.20.13d28fcf96ad656080ef6fd9374f0774.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_21_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_22_JOB_ID: gatk_indel_realigner.S1.21
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.21
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.21.e6e1bbc47200cdcdd45869feb7375671.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.21.e6e1bbc47200cdcdd45869feb7375671.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.21.intervals \
  --intervals 22 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.21.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.21.bam \
  --intervals 22 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.21.e6e1bbc47200cdcdd45869feb7375671.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_22_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_indel_realigner_23_JOB_ID: gatk_indel_realigner.S1.others
#-------------------------------------------------------------------------------
JOB_NAME=gatk_indel_realigner.S1.others
JOB_DEPENDENCIES=$sambamba_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/gatk_indel_realigner/gatk_indel_realigner.S1.others.b2890a02a597ad7b543f57006d4ee2b1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_indel_realigner.S1.others.b2890a02a597ad7b543f57006d4ee2b1.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/realign && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type RealignerTargetCreator -nct 1 -nt 3 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --known /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/realign/S1.sorted.realigned.others.intervals \
  --excludeIntervals 1 \
  --excludeIntervals 2 \
  --excludeIntervals 3 \
  --excludeIntervals 4 \
  --excludeIntervals 5 \
  --excludeIntervals 6 \
  --excludeIntervals 7 \
  --excludeIntervals 8 \
  --excludeIntervals 9 \
  --excludeIntervals 10 \
  --excludeIntervals 11 \
  --excludeIntervals 12 \
  --excludeIntervals 13 \
  --excludeIntervals 14 \
  --excludeIntervals 15 \
  --excludeIntervals 16 \
  --excludeIntervals 17 \
  --excludeIntervals 18 \
  --excludeIntervals 19 \
  --excludeIntervals 20 \
  --excludeIntervals 21 \
  --excludeIntervals 22 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx6G -jar $GATK_JAR \
  --analysis_type IndelRealigner -nt 1 -nct 1 \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --input_file alignment/S1/S1.sorted.bam \
   \
  --targetIntervals alignment/S1/realign/S1.sorted.realigned.others.intervals \
  --knownAlleles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
   \
  --out alignment/S1/realign/S1.sorted.realigned.others.bam \
  --excludeIntervals 1 \
  --excludeIntervals 2 \
  --excludeIntervals 3 \
  --excludeIntervals 4 \
  --excludeIntervals 5 \
  --excludeIntervals 6 \
  --excludeIntervals 7 \
  --excludeIntervals 8 \
  --excludeIntervals 9 \
  --excludeIntervals 10 \
  --excludeIntervals 11 \
  --excludeIntervals 12 \
  --excludeIntervals 13 \
  --excludeIntervals 14 \
  --excludeIntervals 15 \
  --excludeIntervals 16 \
  --excludeIntervals 17 \
  --excludeIntervals 18 \
  --excludeIntervals 19 \
  --excludeIntervals 20 \
  --excludeIntervals 21 \
  --excludeIntervals 22 \
  --maxReadsInMemory 750000
gatk_indel_realigner.S1.others.b2890a02a597ad7b543f57006d4ee2b1.mugqic.done
chmod 755 $COMMAND

gatk_indel_realigner_23_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_indel_realigner_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: sambamba_merge_realigned
#-------------------------------------------------------------------------------
STEP=sambamba_merge_realigned

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: sambamba_merge_realigned_1_JOB_ID: sambamba_merge_realigned.S1
#-------------------------------------------------------------------------------
JOB_NAME=sambamba_merge_realigned.S1
JOB_DEPENDENCIES=$gatk_indel_realigner_1_JOB_ID:$gatk_indel_realigner_2_JOB_ID:$gatk_indel_realigner_3_JOB_ID:$gatk_indel_realigner_4_JOB_ID:$gatk_indel_realigner_5_JOB_ID:$gatk_indel_realigner_6_JOB_ID:$gatk_indel_realigner_7_JOB_ID:$gatk_indel_realigner_8_JOB_ID:$gatk_indel_realigner_9_JOB_ID:$gatk_indel_realigner_10_JOB_ID:$gatk_indel_realigner_11_JOB_ID:$gatk_indel_realigner_12_JOB_ID:$gatk_indel_realigner_13_JOB_ID:$gatk_indel_realigner_14_JOB_ID:$gatk_indel_realigner_15_JOB_ID:$gatk_indel_realigner_16_JOB_ID:$gatk_indel_realigner_17_JOB_ID:$gatk_indel_realigner_18_JOB_ID:$gatk_indel_realigner_19_JOB_ID:$gatk_indel_realigner_20_JOB_ID:$gatk_indel_realigner_21_JOB_ID:$gatk_indel_realigner_22_JOB_ID:$gatk_indel_realigner_23_JOB_ID
JOB_DONE=job_output/sambamba_merge_realigned/sambamba_merge_realigned.S1.3933b2a0407cb3d3f9be6633e6ad863e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sambamba_merge_realigned.S1.3933b2a0407cb3d3f9be6633e6ad863e.mugqic.done' > $COMMAND
module purge && \
module load mugqic/samtools/1.4.1 mugqic/sambamba/0.6.6 && \
sambamba merge -t 7 \
  alignment/S1/S1.realigned.sorted.bam \
   \
  alignment/S1/realign/S1.sorted.realigned.0.bam \
  alignment/S1/realign/S1.sorted.realigned.1.bam \
  alignment/S1/realign/S1.sorted.realigned.2.bam \
  alignment/S1/realign/S1.sorted.realigned.3.bam \
  alignment/S1/realign/S1.sorted.realigned.4.bam \
  alignment/S1/realign/S1.sorted.realigned.5.bam \
  alignment/S1/realign/S1.sorted.realigned.6.bam \
  alignment/S1/realign/S1.sorted.realigned.7.bam \
  alignment/S1/realign/S1.sorted.realigned.8.bam \
  alignment/S1/realign/S1.sorted.realigned.9.bam \
  alignment/S1/realign/S1.sorted.realigned.10.bam \
  alignment/S1/realign/S1.sorted.realigned.11.bam \
  alignment/S1/realign/S1.sorted.realigned.12.bam \
  alignment/S1/realign/S1.sorted.realigned.13.bam \
  alignment/S1/realign/S1.sorted.realigned.14.bam \
  alignment/S1/realign/S1.sorted.realigned.15.bam \
  alignment/S1/realign/S1.sorted.realigned.16.bam \
  alignment/S1/realign/S1.sorted.realigned.17.bam \
  alignment/S1/realign/S1.sorted.realigned.18.bam \
  alignment/S1/realign/S1.sorted.realigned.19.bam \
  alignment/S1/realign/S1.sorted.realigned.20.bam \
  alignment/S1/realign/S1.sorted.realigned.21.bam \
  alignment/S1/realign/S1.sorted.realigned.others.bam
sambamba_merge_realigned.S1.3933b2a0407cb3d3f9be6633e6ad863e.mugqic.done
chmod 755 $COMMAND

sambamba_merge_realigned_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=32G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sambamba_merge_realigned_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: fix_mate_by_coordinate
#-------------------------------------------------------------------------------
STEP=fix_mate_by_coordinate

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: fix_mate_by_coordinate_1_JOB_ID: fix_mate_by_coordinate.S1
#-------------------------------------------------------------------------------
JOB_NAME=fix_mate_by_coordinate.S1
JOB_DEPENDENCIES=$sambamba_merge_realigned_1_JOB_ID
JOB_DONE=job_output/fix_mate_by_coordinate/fix_mate_by_coordinate.S1.ae0e9d4cda123882ca6d9c45cd38adea.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'fix_mate_by_coordinate.S1.ae0e9d4cda123882ca6d9c45cd38adea.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/picard/2.9.0 && \
java -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx16G -jar $BVATOOLS_JAR \
  groupfixmate \
  --level 1 \
  --bam alignment/S1/S1.realigned.sorted.bam \
  --out alignment/S1/S1.matefixed.sorted.tmp.bam && \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx16G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/S1/S1.matefixed.sorted.tmp.bam \
 OUTPUT=alignment/S1/S1.matefixed.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=3750000
fix_mate_by_coordinate.S1.ae0e9d4cda123882ca6d9c45cd38adea.mugqic.done
chmod 755 $COMMAND

fix_mate_by_coordinate_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=71:00:00 --mem=52G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$fix_mate_by_coordinate_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.S1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.S1
JOB_DEPENDENCIES=$fix_mate_by_coordinate_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.S1.ea2add593a70d94af6b32ccb2d0cae20.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_mark_duplicates.S1.ea2add593a70d94af6b32ccb2d0cae20.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=2 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/S1/S1.matefixed.sorted.bam \
 OUTPUT=alignment/S1/S1.sorted.dup.bam \
 METRICS_FILE=alignment/S1/S1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=2500000 
picard_mark_duplicates.S1.ea2add593a70d94af6b32ccb2d0cae20.mugqic.done
chmod 755 $COMMAND

picard_mark_duplicates_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=71:00:00 --mem=12G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: recalibration
#-------------------------------------------------------------------------------
STEP=recalibration

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: recalibration_1_JOB_ID: gatk_base_recalibrator.S1
#-------------------------------------------------------------------------------
JOB_NAME=gatk_base_recalibrator.S1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/recalibration/gatk_base_recalibrator.S1.94810674f0767f2364907e78bf4d7df1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_base_recalibrator.S1.94810674f0767f2364907e78bf4d7df1.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx55G -jar $GATK_JAR \
  --analysis_type BaseRecalibrator --bqsrBAQGapOpenPenalty 30 \
  -nt 1 --num_cpu_threads_per_data_thread 12 \
  --input_file alignment/S1/S1.sorted.dup.bam \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa  \
  --knownSites /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.dbSNP142.vcf.gz \
  --knownSites /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.gnomad.exomes.r2.0.1.sites.no-VEP.nohist.tidy.vcf.gz \
  --knownSites /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Mills_and_1000G_gold_standard.indels.vcf.gz \
  --out alignment/S1/S1.sorted.dup.recalibration_report.grp
gatk_base_recalibrator.S1.94810674f0767f2364907e78bf4d7df1.mugqic.done
chmod 755 $COMMAND

recalibration_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=60G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$recalibration_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: recalibration_2_JOB_ID: gatk_print_reads.S1
#-------------------------------------------------------------------------------
JOB_NAME=gatk_print_reads.S1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$recalibration_1_JOB_ID
JOB_DONE=job_output/recalibration/gatk_print_reads.S1.45ee3ac365ea54f41784225cfb91107e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_print_reads.S1.45ee3ac365ea54f41784225cfb91107e.mugqic.done' > $COMMAND
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
gatk_print_reads.S1.45ee3ac365ea54f41784225cfb91107e.mugqic.done
chmod 755 $COMMAND

recalibration_2_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=96:00:00 --mem=60G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$recalibration_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: sym_link_final_bam
#-------------------------------------------------------------------------------
STEP=sym_link_final_bam

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: sym_link_final_bam_1_JOB_ID: sym_link_final_bam.S1
#-------------------------------------------------------------------------------
JOB_NAME=sym_link_final_bam.S1
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/sym_link_final_bam/sym_link_final_bam.S1.3d956516c2df893d086b3ad9d3e31754.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'sym_link_final_bam.S1.3d956516c2df893d086b3ad9d3e31754.mugqic.done' > $COMMAND
md5sum /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/alignment/S1/S1.sorted.dup.recal.bam > /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/alignment/S1/S1.sorted.dup.recal.bam.md5 && \
mkdir -p /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/deliverables/S1/wgs/alignment && \       
ln -sf \
  /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/alignment/S1/S1.sorted.dup.recal.bam \
  /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/deliverables/S1/wgs/alignment/S1.sorted.dup.recal.bam && \
mkdir -p /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/deliverables/S1/wgs/alignment && \       
ln -sf \
  /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/alignment/S1/S1.sorted.dup.recal.bai \
  /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/deliverables/S1/wgs/alignment/S1.sorted.dup.recal.bai && \
mkdir -p /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/deliverables/S1/wgs/alignment && \       
ln -sf \
  /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/alignment/S1/S1.sorted.dup.recal.bam.md5 \
  /lustre03/project/6007512/rdali/C3G/projects/dnaseqTestSet/deliverables/S1/wgs/alignment/S1.sorted.dup.recal.bam.md5
sym_link_final_bam.S1.3d956516c2df893d086b3ad9d3e31754.mugqic.done
chmod 755 $COMMAND

sym_link_final_bam_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu=4000M -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$sym_link_final_bam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: metrics_dna_picard_metrics
#-------------------------------------------------------------------------------
STEP=metrics_dna_picard_metrics

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: metrics_dna_picard_metrics_1_JOB_ID: picard_collect_multiple_metrics.S1
#-------------------------------------------------------------------------------
JOB_NAME=picard_collect_multiple_metrics.S1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/metrics_dna_picard_metrics/picard_collect_multiple_metrics.S1.94d4b5202bd3c1f5eb7c6eeaa84c1634.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_collect_multiple_metrics.S1.94d4b5202bd3c1f5eb7c6eeaa84c1634.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.9.0 mugqic/R_Bioconductor/3.4.1_3.5 && \
mkdir -p metrics/dna/S1/picard_metrics && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx6G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
 INPUT=alignment/S1/S1.sorted.dup.bam \
 OUTPUT=metrics/dna/S1/picard_metrics.all.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_collect_multiple_metrics.S1.94d4b5202bd3c1f5eb7c6eeaa84c1634.mugqic.done
chmod 755 $COMMAND

metrics_dna_picard_metrics_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_dna_picard_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: metrics_dna_picard_metrics_2_JOB_ID: picard_collect_oxog_metrics.S1
#-------------------------------------------------------------------------------
JOB_NAME=picard_collect_oxog_metrics.S1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/metrics_dna_picard_metrics/picard_collect_oxog_metrics.S1.2821ab427935ffa4df62c825072f61a9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_collect_oxog_metrics.S1.2821ab427935ffa4df62c825072f61a9.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.9.0 mugqic/R_Bioconductor/3.4.1_3.5 && \
mkdir -p metrics/dna/S1/picard_metrics && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx6G -jar $PICARD_HOME/picard.jar CollectOxoGMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/S1/S1.sorted.dup.bam \
 OUTPUT=metrics/dna/S1/picard_metrics.oxog_metrics.txt \
 DB_SNP=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.dbSNP142.vcf.gz \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
 MAX_RECORDS_IN_RAM=4000000
picard_collect_oxog_metrics.S1.2821ab427935ffa4df62c825072f61a9.mugqic.done
chmod 755 $COMMAND

metrics_dna_picard_metrics_2_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=8G -N 1 -n 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_dna_picard_metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: metrics_dna_picard_metrics_3_JOB_ID: picard_collect_gcbias_metrics.S1
#-------------------------------------------------------------------------------
JOB_NAME=picard_collect_gcbias_metrics.S1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/metrics_dna_picard_metrics/picard_collect_gcbias_metrics.S1.57e90f0735b5dcbb5e3d8a54081b54c3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_collect_gcbias_metrics.S1.57e90f0735b5dcbb5e3d8a54081b54c3.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.9.0 mugqic/R_Bioconductor/3.4.1_3.5 && \
mkdir -p metrics/dna/S1/picard_metrics && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx6G -jar $PICARD_HOME/picard.jar CollectGcBiasMetrics \
 VALIDATION_STRINGENCY=SILENT ALSO_IGNORE_DUPLICATES=TRUE \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/S1/S1.sorted.dup.bam \
 OUTPUT=metrics/dna/S1/picard_metrics.qcbias_metrics.txt \
 CHART=metrics/dna/S1/picard_metrics.qcbias_metrics.pdf \
 SUMMARY_OUTPUT=metrics/dna/S1/picard_metrics.qcbias_summary_metrics.txt \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
 MAX_RECORDS_IN_RAM=4000000
picard_collect_gcbias_metrics.S1.57e90f0735b5dcbb5e3d8a54081b54c3.mugqic.done
chmod 755 $COMMAND

metrics_dna_picard_metrics_3_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=8G -N 1 -n 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_dna_picard_metrics_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: metrics_dna_sample_qualimap
#-------------------------------------------------------------------------------
STEP=metrics_dna_sample_qualimap

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: metrics_dna_sample_qualimap_1_JOB_ID: dna_sample_qualimap.S1
#-------------------------------------------------------------------------------
JOB_NAME=dna_sample_qualimap.S1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/metrics_dna_sample_qualimap/dna_sample_qualimap.S1.09d770dd7c066a82059843f806a26d4b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'dna_sample_qualimap.S1.09d770dd7c066a82059843f806a26d4b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/qualimap/2.2.1 && \
mkdir -p metrics/dna/S1/qualimap/S1 && \
qualimap bamqc -nt 11 -gd HUMAN \
  -bam alignment/S1/S1.sorted.dup.bam -outdir metrics/dna/S1/qualimap/S1 \
  --java-mem-size=55G
dna_sample_qualimap.S1.09d770dd7c066a82059843f806a26d4b.mugqic.done
chmod 755 $COMMAND

metrics_dna_sample_qualimap_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=96:00:00 --mem=60G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_dna_sample_qualimap_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: metrics_dna_sambamba_flagstat
#-------------------------------------------------------------------------------
STEP=metrics_dna_sambamba_flagstat

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: metrics_dna_sambamba_flagstat_1_JOB_ID: dna_sambamba_flagstat.S1
#-------------------------------------------------------------------------------
JOB_NAME=dna_sambamba_flagstat.S1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/metrics_dna_sambamba_flagstat/dna_sambamba_flagstat.S1.fbcf47143659e44d51cf305888379fd5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'dna_sambamba_flagstat.S1.fbcf47143659e44d51cf305888379fd5.mugqic.done' > $COMMAND
module purge && \
module load mugqic/sambamba/0.6.6 && \
mkdir -p metrics/dna/S1/flagstat && \
sambamba flagstat -t 5 \
  alignment/S1/S1.sorted.dup.bam \
  > metrics/dna/S1/flagstat/S1.flagstat
dna_sambamba_flagstat.S1.fbcf47143659e44d51cf305888379fd5.mugqic.done
chmod 755 $COMMAND

metrics_dna_sambamba_flagstat_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:00 --mem-per-cpu=4000M -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_dna_sambamba_flagstat_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: metrics_dna_fastqc
#-------------------------------------------------------------------------------
STEP=metrics_dna_fastqc

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: metrics_dna_fastqc_1_JOB_ID: fastqc.S1
#-------------------------------------------------------------------------------
JOB_NAME=fastqc.S1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/metrics_dna_fastqc/fastqc.S1.68a1a3d6c8d9499367e235ede50b5f03.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'fastqc.S1.68a1a3d6c8d9499367e235ede50b5f03.mugqic.done' > $COMMAND
module purge && \
module load mugqic/fastqc/0.11.5 mugqic/java/openjdk-jdk1.8.0_72 && \
mkdir -p metrics/dna/S1/fastqc && \
`cat > metrics/dna/S1/fastqc/adapter.tsv << END
>Adapter1	AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>Adapter2	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
END` && \
fastqc \
  -o metrics/dna/S1/fastqc \
  -t 3 \
  -a metrics/dna/S1/fastqc/adapter.tsv \
  -f bam \
  alignment/S1/S1.sorted.dup.bam
fastqc.S1.68a1a3d6c8d9499367e235ede50b5f03.mugqic.done
chmod 755 $COMMAND

metrics_dna_fastqc_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem-per-cpu=4000M -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_dna_fastqc_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: gatk_callable_loci
#-------------------------------------------------------------------------------
STEP=gatk_callable_loci

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: gatk_callable_loci_1_JOB_ID: gatk_callable_loci.S1
#-------------------------------------------------------------------------------
JOB_NAME=gatk_callable_loci.S1
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_callable_loci/gatk_callable_loci.S1.d0f82fb33c4694cae01e690ec9fbc854.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_callable_loci.S1.d0f82fb33c4694cae01e690ec9fbc854.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx10G -jar $GATK_JAR \
  --analysis_type CallableLoci -dt none --minDepth 10 --maxDepth 500 --minDepthForLowMAPQ 10 --minMappingQuality 10 --minBaseQuality 15 \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --summary alignment/S1/S1.callable.summary.txt \
  --out alignment/S1/S1.callable.bed
gatk_callable_loci.S1.d0f82fb33c4694cae01e690ec9fbc854.mugqic.done
chmod 755 $COMMAND

gatk_callable_loci_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem-per-cpu=4775M -N 1 -n 3 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_callable_loci_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: extract_common_snp_freq
#-------------------------------------------------------------------------------
STEP=extract_common_snp_freq

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: extract_common_snp_freq_1_JOB_ID: extract_common_snp_freq.S1
#-------------------------------------------------------------------------------
JOB_NAME=extract_common_snp_freq.S1
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/extract_common_snp_freq/extract_common_snp_freq.S1.bf6e21de5d21dad8315d809893cf01f0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'extract_common_snp_freq.S1.bf6e21de5d21dad8315d809893cf01f0.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 && \
java -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx8G -jar $BVATOOLS_JAR \
  basefreq \
  --pos /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/common.dbsnp132.q60.tsv \
  --bam alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/S1.commonSNPs.alleleFreq.csv
extract_common_snp_freq.S1.bf6e21de5d21dad8315d809893cf01f0.mugqic.done
chmod 755 $COMMAND

extract_common_snp_freq_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=20G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$extract_common_snp_freq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: baf_plot
#-------------------------------------------------------------------------------
STEP=baf_plot

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: baf_plot_1_JOB_ID: baf_plot.S1
#-------------------------------------------------------------------------------
JOB_NAME=baf_plot.S1
JOB_DEPENDENCIES=$extract_common_snp_freq_1_JOB_ID
JOB_DONE=job_output/baf_plot/baf_plot.S1.02a9cc8aceed2f635f84e7ed34843001.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'baf_plot.S1.02a9cc8aceed2f635f84e7ed34843001.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 && \
java -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx70G -jar $BVATOOLS_JAR \
  ratiobaf --plot --maxDepth 1000  --exclude MT,GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1 \
  --refdict /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.dict \
  --snppos /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/common.dbsnp132.q60.tsv \
  --basefreq alignment/S1/S1.commonSNPs.alleleFreq.csv \
  --prefix alignment/S1/S1.ratioBAF
baf_plot.S1.02a9cc8aceed2f635f84e7ed34843001.mugqic.done
chmod 755 $COMMAND

baf_plot_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$baf_plot_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: gatk_haplotype_caller
#-------------------------------------------------------------------------------
STEP=gatk_haplotype_caller

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_1_JOB_ID: gatk_haplotype_caller.S1.0
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.0
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.0.9c34aa46055e4b83c4f9d59be7a5b4d2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.0.9c34aa46055e4b83c4f9d59be7a5b4d2.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.0.hc.g.vcf.gz \
  --intervals 1
gatk_haplotype_caller.S1.0.9c34aa46055e4b83c4f9d59be7a5b4d2.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_2_JOB_ID: gatk_haplotype_caller.S1.1
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.1
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.1.452206f7f6416c89e5a2391879559a15.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.1.452206f7f6416c89e5a2391879559a15.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.1.hc.g.vcf.gz \
  --intervals 2
gatk_haplotype_caller.S1.1.452206f7f6416c89e5a2391879559a15.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_2_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_3_JOB_ID: gatk_haplotype_caller.S1.2
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.2
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.2.425068b29a4a3dfdfebea9a41cd016bd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.2.425068b29a4a3dfdfebea9a41cd016bd.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.2.hc.g.vcf.gz \
  --intervals 3
gatk_haplotype_caller.S1.2.425068b29a4a3dfdfebea9a41cd016bd.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_3_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_4_JOB_ID: gatk_haplotype_caller.S1.3
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.3
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.3.dd20d3885c7ebe5f014ff3ac6c778f5d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.3.dd20d3885c7ebe5f014ff3ac6c778f5d.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.3.hc.g.vcf.gz \
  --intervals 4
gatk_haplotype_caller.S1.3.dd20d3885c7ebe5f014ff3ac6c778f5d.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_4_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_5_JOB_ID: gatk_haplotype_caller.S1.4
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.4
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.4.a994c6586068001e8f0c39a1da7799f2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.4.a994c6586068001e8f0c39a1da7799f2.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.4.hc.g.vcf.gz \
  --intervals 5
gatk_haplotype_caller.S1.4.a994c6586068001e8f0c39a1da7799f2.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_5_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_6_JOB_ID: gatk_haplotype_caller.S1.5
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.5
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.5.fe990e603ceab7db367389258be8edad.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.5.fe990e603ceab7db367389258be8edad.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.5.hc.g.vcf.gz \
  --intervals 6
gatk_haplotype_caller.S1.5.fe990e603ceab7db367389258be8edad.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_6_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_7_JOB_ID: gatk_haplotype_caller.S1.6
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.6
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.6.fcfe464e6477bcfd7154f457f979d953.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.6.fcfe464e6477bcfd7154f457f979d953.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.6.hc.g.vcf.gz \
  --intervals 7
gatk_haplotype_caller.S1.6.fcfe464e6477bcfd7154f457f979d953.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_7_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_8_JOB_ID: gatk_haplotype_caller.S1.7
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.7
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.7.99b08b4a25eb12261ef3e8fb17247c71.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.7.99b08b4a25eb12261ef3e8fb17247c71.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.7.hc.g.vcf.gz \
  --intervals 8
gatk_haplotype_caller.S1.7.99b08b4a25eb12261ef3e8fb17247c71.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_8_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_9_JOB_ID: gatk_haplotype_caller.S1.8
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.8
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.8.af8c75eac64df11ade613fcd79a46123.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.8.af8c75eac64df11ade613fcd79a46123.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.8.hc.g.vcf.gz \
  --intervals 9
gatk_haplotype_caller.S1.8.af8c75eac64df11ade613fcd79a46123.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_9_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_10_JOB_ID: gatk_haplotype_caller.S1.9
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.9
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.9.199669486f227845f31aa27732709957.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.9.199669486f227845f31aa27732709957.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.9.hc.g.vcf.gz \
  --intervals 10
gatk_haplotype_caller.S1.9.199669486f227845f31aa27732709957.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_10_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_11_JOB_ID: gatk_haplotype_caller.S1.10
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.10
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.10.cfb8787cf5ca9e2b7dcb922b9b658b21.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.10.cfb8787cf5ca9e2b7dcb922b9b658b21.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.10.hc.g.vcf.gz \
  --intervals 11
gatk_haplotype_caller.S1.10.cfb8787cf5ca9e2b7dcb922b9b658b21.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_11_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_12_JOB_ID: gatk_haplotype_caller.S1.11
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.11
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.11.5bf2bbb72c017c478eb3c465bf4e253e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.11.5bf2bbb72c017c478eb3c465bf4e253e.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.11.hc.g.vcf.gz \
  --intervals 12
gatk_haplotype_caller.S1.11.5bf2bbb72c017c478eb3c465bf4e253e.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_12_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_13_JOB_ID: gatk_haplotype_caller.S1.12
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.12
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.12.e3926d633803ff0aa48bcc87b906f1b9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.12.e3926d633803ff0aa48bcc87b906f1b9.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.12.hc.g.vcf.gz \
  --intervals 13
gatk_haplotype_caller.S1.12.e3926d633803ff0aa48bcc87b906f1b9.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_13_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_14_JOB_ID: gatk_haplotype_caller.S1.13
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.13
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.13.1b6816ef4401f358ebf1e57d2a4bdc1d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.13.1b6816ef4401f358ebf1e57d2a4bdc1d.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.13.hc.g.vcf.gz \
  --intervals 14
gatk_haplotype_caller.S1.13.1b6816ef4401f358ebf1e57d2a4bdc1d.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_14_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_15_JOB_ID: gatk_haplotype_caller.S1.14
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.14
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.14.3eeff25f4cad061f6199dfe30fb10471.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.14.3eeff25f4cad061f6199dfe30fb10471.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.14.hc.g.vcf.gz \
  --intervals 15
gatk_haplotype_caller.S1.14.3eeff25f4cad061f6199dfe30fb10471.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_15_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_16_JOB_ID: gatk_haplotype_caller.S1.15
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.15
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.15.d15d01c3941d775c1d26ed1c5e5d59e3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.15.d15d01c3941d775c1d26ed1c5e5d59e3.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.15.hc.g.vcf.gz \
  --intervals 16
gatk_haplotype_caller.S1.15.d15d01c3941d775c1d26ed1c5e5d59e3.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_16_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_17_JOB_ID: gatk_haplotype_caller.S1.16
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.16
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.16.f31cf4a7e414600dad9b679ebc631a11.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.16.f31cf4a7e414600dad9b679ebc631a11.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.16.hc.g.vcf.gz \
  --intervals 17
gatk_haplotype_caller.S1.16.f31cf4a7e414600dad9b679ebc631a11.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_17_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_18_JOB_ID: gatk_haplotype_caller.S1.17
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.17
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.17.855a7a1fc7065348633c59a0c15f35c8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.17.855a7a1fc7065348633c59a0c15f35c8.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.17.hc.g.vcf.gz \
  --intervals 18
gatk_haplotype_caller.S1.17.855a7a1fc7065348633c59a0c15f35c8.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_18_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_19_JOB_ID: gatk_haplotype_caller.S1.18
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.18
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.18.120c384d0d804502205b8d3490ec7fdb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.18.120c384d0d804502205b8d3490ec7fdb.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.18.hc.g.vcf.gz \
  --intervals 19
gatk_haplotype_caller.S1.18.120c384d0d804502205b8d3490ec7fdb.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_19_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_20_JOB_ID: gatk_haplotype_caller.S1.19
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.19
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.19.8070345a6dc0e8796ba0a19d83ea4ccb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.19.8070345a6dc0e8796ba0a19d83ea4ccb.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.19.hc.g.vcf.gz \
  --intervals 20
gatk_haplotype_caller.S1.19.8070345a6dc0e8796ba0a19d83ea4ccb.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_20_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_21_JOB_ID: gatk_haplotype_caller.S1.20
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.20
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.20.59111a76bca611100495fbc094e3ff32.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.20.59111a76bca611100495fbc094e3ff32.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.20.hc.g.vcf.gz \
  --intervals 21
gatk_haplotype_caller.S1.20.59111a76bca611100495fbc094e3ff32.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_21_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_22_JOB_ID: gatk_haplotype_caller.S1.21
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.21
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.21.3a692c8be296a26506e191f699548212.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.21.3a692c8be296a26506e191f699548212.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.21.hc.g.vcf.gz \
  --intervals 22
gatk_haplotype_caller.S1.21.3a692c8be296a26506e191f699548212.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_22_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gatk_haplotype_caller_23_JOB_ID: gatk_haplotype_caller.S1.others
#-------------------------------------------------------------------------------
JOB_NAME=gatk_haplotype_caller.S1.others
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/gatk_haplotype_caller/gatk_haplotype_caller.S1.others.84302517c73e932faf97b64e902d4cf3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_haplotype_caller.S1.others.84302517c73e932faf97b64e902d4cf3.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p alignment/S1/rawHaplotypeCaller && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type HaplotypeCaller --useNewAFCalculator --emitRefConfidence GVCF -dt none -nct 1 -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --input_file alignment/S1/S1.sorted.dup.recal.bam \
  --out alignment/S1/rawHaplotypeCaller/S1.others.hc.g.vcf.gz \
  --excludeIntervals 1 \
  --excludeIntervals 2 \
  --excludeIntervals 3 \
  --excludeIntervals 4 \
  --excludeIntervals 5 \
  --excludeIntervals 6 \
  --excludeIntervals 7 \
  --excludeIntervals 8 \
  --excludeIntervals 9 \
  --excludeIntervals 10 \
  --excludeIntervals 11 \
  --excludeIntervals 12 \
  --excludeIntervals 13 \
  --excludeIntervals 14 \
  --excludeIntervals 15 \
  --excludeIntervals 16 \
  --excludeIntervals 17 \
  --excludeIntervals 18 \
  --excludeIntervals 19 \
  --excludeIntervals 20 \
  --excludeIntervals 21 \
  --excludeIntervals 22 \
  --excludeIntervals GL000207.1 \
  --excludeIntervals GL000226.1 \
  --excludeIntervals GL000229.1 \
  --excludeIntervals GL000231.1 \
  --excludeIntervals GL000210.1 \
  --excludeIntervals GL000239.1 \
  --excludeIntervals GL000235.1 \
  --excludeIntervals GL000201.1 \
  --excludeIntervals GL000247.1 \
  --excludeIntervals GL000245.1 \
  --excludeIntervals GL000197.1 \
  --excludeIntervals GL000203.1 \
  --excludeIntervals GL000246.1 \
  --excludeIntervals GL000249.1 \
  --excludeIntervals GL000196.1 \
  --excludeIntervals GL000248.1 \
  --excludeIntervals GL000244.1 \
  --excludeIntervals GL000238.1 \
  --excludeIntervals GL000202.1 \
  --excludeIntervals GL000234.1 \
  --excludeIntervals GL000232.1 \
  --excludeIntervals GL000206.1 \
  --excludeIntervals GL000240.1 \
  --excludeIntervals GL000236.1 \
  --excludeIntervals GL000241.1 \
  --excludeIntervals GL000243.1 \
  --excludeIntervals GL000242.1 \
  --excludeIntervals GL000230.1 \
  --excludeIntervals GL000237.1 \
  --excludeIntervals GL000233.1 \
  --excludeIntervals GL000204.1 \
  --excludeIntervals GL000198.1 \
  --excludeIntervals GL000208.1 \
  --excludeIntervals GL000191.1 \
  --excludeIntervals GL000227.1 \
  --excludeIntervals GL000228.1 \
  --excludeIntervals GL000214.1 \
  --excludeIntervals GL000221.1 \
  --excludeIntervals GL000209.1 \
  --excludeIntervals GL000218.1 \
  --excludeIntervals GL000220.1 \
  --excludeIntervals GL000213.1 \
  --excludeIntervals GL000211.1 \
  --excludeIntervals GL000199.1 \
  --excludeIntervals GL000217.1 \
  --excludeIntervals GL000216.1 \
  --excludeIntervals GL000215.1 \
  --excludeIntervals GL000205.1 \
  --excludeIntervals GL000219.1 \
  --excludeIntervals GL000224.1 \
  --excludeIntervals GL000223.1 \
  --excludeIntervals GL000195.1 \
  --excludeIntervals GL000212.1 \
  --excludeIntervals GL000222.1 \
  --excludeIntervals GL000200.1 \
  --excludeIntervals GL000193.1 \
  --excludeIntervals GL000194.1 \
  --excludeIntervals GL000225.1 \
  --excludeIntervals GL000192.1
gatk_haplotype_caller.S1.others.84302517c73e932faf97b64e902d4cf3.mugqic.done
chmod 755 $COMMAND

gatk_haplotype_caller_23_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=36G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$gatk_haplotype_caller_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: merge_and_call_individual_gvcf
#-------------------------------------------------------------------------------
STEP=merge_and_call_individual_gvcf

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: merge_and_call_individual_gvcf_1_JOB_ID: merge_and_call_individual_gvcf.merge.S1
#-------------------------------------------------------------------------------
JOB_NAME=merge_and_call_individual_gvcf.merge.S1
JOB_DEPENDENCIES=$gatk_haplotype_caller_1_JOB_ID:$gatk_haplotype_caller_2_JOB_ID:$gatk_haplotype_caller_3_JOB_ID:$gatk_haplotype_caller_4_JOB_ID:$gatk_haplotype_caller_5_JOB_ID:$gatk_haplotype_caller_6_JOB_ID:$gatk_haplotype_caller_7_JOB_ID:$gatk_haplotype_caller_8_JOB_ID:$gatk_haplotype_caller_9_JOB_ID:$gatk_haplotype_caller_10_JOB_ID:$gatk_haplotype_caller_11_JOB_ID:$gatk_haplotype_caller_12_JOB_ID:$gatk_haplotype_caller_13_JOB_ID:$gatk_haplotype_caller_14_JOB_ID:$gatk_haplotype_caller_15_JOB_ID:$gatk_haplotype_caller_16_JOB_ID:$gatk_haplotype_caller_17_JOB_ID:$gatk_haplotype_caller_18_JOB_ID:$gatk_haplotype_caller_19_JOB_ID:$gatk_haplotype_caller_20_JOB_ID:$gatk_haplotype_caller_21_JOB_ID:$gatk_haplotype_caller_22_JOB_ID:$gatk_haplotype_caller_23_JOB_ID
JOB_DONE=job_output/merge_and_call_individual_gvcf/merge_and_call_individual_gvcf.merge.S1.1d39c5064bd1133952aa9fc0103d93f0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'merge_and_call_individual_gvcf.merge.S1.1d39c5064bd1133952aa9fc0103d93f0.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx6G -cp $GATK_JAR \
  org.broadinstitute.gatk.tools.CatVariants  \
  --reference /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --variant alignment/S1/rawHaplotypeCaller/S1.0.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.1.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.2.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.3.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.4.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.5.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.6.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.7.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.8.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.9.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.10.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.11.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.12.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.13.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.14.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.15.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.16.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.17.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.18.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.19.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.20.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.21.hc.g.vcf.gz \
  --variant alignment/S1/rawHaplotypeCaller/S1.others.hc.g.vcf.gz \
  --outputFile alignment/S1/S1.hc.g.vcf.gz
merge_and_call_individual_gvcf.merge.S1.1d39c5064bd1133952aa9fc0103d93f0.mugqic.done
chmod 755 $COMMAND

merge_and_call_individual_gvcf_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=120:00:00 --mem=30G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$merge_and_call_individual_gvcf_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: merge_and_call_individual_gvcf_2_JOB_ID: merge_and_call_individual_gvcf.call.S1
#-------------------------------------------------------------------------------
JOB_NAME=merge_and_call_individual_gvcf.call.S1
JOB_DEPENDENCIES=$merge_and_call_individual_gvcf_1_JOB_ID
JOB_DONE=job_output/merge_and_call_individual_gvcf/merge_and_call_individual_gvcf.call.S1.91d9c21d247b55dae94e594c2b705250.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'merge_and_call_individual_gvcf.call.S1.91d9c21d247b55dae94e594c2b705250.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type GenotypeGVCFs --useNewAFCalculator -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --variant alignment/S1/S1.hc.g.vcf.gz \
  --out alignment/S1/S1.hc.vcf.gz
merge_and_call_individual_gvcf.call.S1.91d9c21d247b55dae94e594c2b705250.mugqic.done
chmod 755 $COMMAND

merge_and_call_individual_gvcf_2_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=120:00:00 --mem=30G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$merge_and_call_individual_gvcf_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: combine_gvcf
#-------------------------------------------------------------------------------
STEP=combine_gvcf

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: combine_gvcf_1_JOB_ID: gatk_combine_gvcf.AllSample.0
#-------------------------------------------------------------------------------
JOB_NAME=gatk_combine_gvcf.AllSample.0
JOB_DEPENDENCIES=$merge_and_call_individual_gvcf_1_JOB_ID
JOB_DONE=job_output/combine_gvcf/gatk_combine_gvcf.AllSample.0.8325bd922c5898d7707a3d9bc9c6fbf0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_combine_gvcf.AllSample.0.8325bd922c5898d7707a3d9bc9c6fbf0.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p variants && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=4194304 -Xmx24G -jar $GATK_JAR \
  --analysis_type CombineGVCFs  \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --variant alignment/S1/S1.hc.g.vcf.gz \
  --out variants/allSamples.0.hc.g.vcf.bgz \
  --intervals 1 \
  --intervals 2 \
  --intervals 3 \
  --intervals 4 \
  --intervals 5
gatk_combine_gvcf.AllSample.0.8325bd922c5898d7707a3d9bc9c6fbf0.mugqic.done
chmod 755 $COMMAND

combine_gvcf_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=60G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$combine_gvcf_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: combine_gvcf_2_JOB_ID: gatk_combine_gvcf.AllSample.1
#-------------------------------------------------------------------------------
JOB_NAME=gatk_combine_gvcf.AllSample.1
JOB_DEPENDENCIES=$merge_and_call_individual_gvcf_1_JOB_ID
JOB_DONE=job_output/combine_gvcf/gatk_combine_gvcf.AllSample.1.f75ea7c98a7cfce455b2b4ea6ee6657b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_combine_gvcf.AllSample.1.f75ea7c98a7cfce455b2b4ea6ee6657b.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p variants && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=4194304 -Xmx24G -jar $GATK_JAR \
  --analysis_type CombineGVCFs  \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --variant alignment/S1/S1.hc.g.vcf.gz \
  --out variants/allSamples.1.hc.g.vcf.bgz \
  --intervals 6 \
  --intervals 7 \
  --intervals 8 \
  --intervals 9 \
  --intervals 10
gatk_combine_gvcf.AllSample.1.f75ea7c98a7cfce455b2b4ea6ee6657b.mugqic.done
chmod 755 $COMMAND

combine_gvcf_2_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=60G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$combine_gvcf_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: combine_gvcf_3_JOB_ID: gatk_combine_gvcf.AllSample.2
#-------------------------------------------------------------------------------
JOB_NAME=gatk_combine_gvcf.AllSample.2
JOB_DEPENDENCIES=$merge_and_call_individual_gvcf_1_JOB_ID
JOB_DONE=job_output/combine_gvcf/gatk_combine_gvcf.AllSample.2.170e705819895bb11161187896cf463a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_combine_gvcf.AllSample.2.170e705819895bb11161187896cf463a.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
mkdir -p variants && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=4194304 -Xmx24G -jar $GATK_JAR \
  --analysis_type CombineGVCFs  \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --variant alignment/S1/S1.hc.g.vcf.gz \
  --out variants/allSamples.2.hc.g.vcf.bgz \
  --intervals 11 \
  --intervals 12 \
  --intervals 13 \
  --intervals 14
gatk_combine_gvcf.AllSample.2.170e705819895bb11161187896cf463a.mugqic.done
chmod 755 $COMMAND

combine_gvcf_3_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=60G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$combine_gvcf_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: combine_gvcf_4_JOB_ID: gatk_combine_gvcf.AllSample.others
#-------------------------------------------------------------------------------
JOB_NAME=gatk_combine_gvcf.AllSample.others
JOB_DEPENDENCIES=$merge_and_call_individual_gvcf_1_JOB_ID
JOB_DONE=job_output/combine_gvcf/gatk_combine_gvcf.AllSample.others.5865ad3fa049efa6e08318610f2c9c7c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gatk_combine_gvcf.AllSample.others.5865ad3fa049efa6e08318610f2c9c7c.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=4194304 -Xmx24G -jar $GATK_JAR \
  --analysis_type CombineGVCFs  \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
   \
  --variant alignment/S1/S1.hc.g.vcf.gz \
  --out variants/allSamples.others.hc.g.vcf.bgz \
  --excludeIntervals 1 \
  --excludeIntervals 2 \
  --excludeIntervals 3 \
  --excludeIntervals 4 \
  --excludeIntervals 5 \
  --excludeIntervals 6 \
  --excludeIntervals 7 \
  --excludeIntervals 8 \
  --excludeIntervals 9 \
  --excludeIntervals 10 \
  --excludeIntervals 11 \
  --excludeIntervals 12 \
  --excludeIntervals 13 \
  --excludeIntervals 14 \
  --excludeIntervals GL000207.1 \
  --excludeIntervals GL000226.1 \
  --excludeIntervals GL000229.1 \
  --excludeIntervals GL000231.1 \
  --excludeIntervals GL000210.1 \
  --excludeIntervals GL000239.1 \
  --excludeIntervals GL000235.1 \
  --excludeIntervals GL000201.1 \
  --excludeIntervals GL000247.1 \
  --excludeIntervals GL000245.1 \
  --excludeIntervals GL000197.1 \
  --excludeIntervals GL000203.1 \
  --excludeIntervals GL000246.1 \
  --excludeIntervals GL000249.1 \
  --excludeIntervals GL000196.1 \
  --excludeIntervals GL000248.1 \
  --excludeIntervals GL000244.1 \
  --excludeIntervals GL000238.1 \
  --excludeIntervals GL000202.1 \
  --excludeIntervals GL000234.1 \
  --excludeIntervals GL000232.1 \
  --excludeIntervals GL000206.1 \
  --excludeIntervals GL000240.1 \
  --excludeIntervals GL000236.1 \
  --excludeIntervals GL000241.1 \
  --excludeIntervals GL000243.1 \
  --excludeIntervals GL000242.1 \
  --excludeIntervals GL000230.1 \
  --excludeIntervals GL000237.1 \
  --excludeIntervals GL000233.1 \
  --excludeIntervals GL000204.1 \
  --excludeIntervals GL000198.1 \
  --excludeIntervals GL000208.1 \
  --excludeIntervals GL000191.1 \
  --excludeIntervals GL000227.1 \
  --excludeIntervals GL000228.1 \
  --excludeIntervals GL000214.1 \
  --excludeIntervals GL000221.1 \
  --excludeIntervals GL000209.1 \
  --excludeIntervals GL000218.1 \
  --excludeIntervals GL000220.1 \
  --excludeIntervals GL000213.1 \
  --excludeIntervals GL000211.1 \
  --excludeIntervals GL000199.1 \
  --excludeIntervals GL000217.1 \
  --excludeIntervals GL000216.1 \
  --excludeIntervals GL000215.1 \
  --excludeIntervals GL000205.1 \
  --excludeIntervals GL000219.1 \
  --excludeIntervals GL000224.1 \
  --excludeIntervals GL000223.1 \
  --excludeIntervals GL000195.1 \
  --excludeIntervals GL000212.1 \
  --excludeIntervals GL000222.1 \
  --excludeIntervals GL000200.1 \
  --excludeIntervals GL000193.1 \
  --excludeIntervals GL000194.1 \
  --excludeIntervals GL000225.1 \
  --excludeIntervals GL000192.1
gatk_combine_gvcf.AllSample.others.5865ad3fa049efa6e08318610f2c9c7c.mugqic.done
chmod 755 $COMMAND

combine_gvcf_4_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=60G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$combine_gvcf_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: merge_and_call_combined_gvcf
#-------------------------------------------------------------------------------
STEP=merge_and_call_combined_gvcf

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: merge_and_call_combined_gvcf_1_JOB_ID: merge_and_call_combined_gvcf.merge.AllSample
#-------------------------------------------------------------------------------
JOB_NAME=merge_and_call_combined_gvcf.merge.AllSample
JOB_DEPENDENCIES=$combine_gvcf_1_JOB_ID:$combine_gvcf_2_JOB_ID:$combine_gvcf_3_JOB_ID:$combine_gvcf_4_JOB_ID
JOB_DONE=job_output/merge_and_call_combined_gvcf/merge_and_call_combined_gvcf.merge.AllSample.87ad03d0f1613e03ce8793a82abc2cc0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'merge_and_call_combined_gvcf.merge.AllSample.87ad03d0f1613e03ce8793a82abc2cc0.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx6G -cp $GATK_JAR \
  org.broadinstitute.gatk.tools.CatVariants  \
  --reference /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --variant variants/allSamples.0.hc.g.vcf.bgz \
  --variant variants/allSamples.1.hc.g.vcf.bgz \
  --variant variants/allSamples.2.hc.g.vcf.bgz \
  --variant variants/allSamples.others.hc.g.vcf.bgz \
  --outputFile variants/allSamples.hc.g.vcf.gz
merge_and_call_combined_gvcf.merge.AllSample.87ad03d0f1613e03ce8793a82abc2cc0.mugqic.done
chmod 755 $COMMAND

merge_and_call_combined_gvcf_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=120:00:00 --mem=60G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$merge_and_call_combined_gvcf_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: merge_and_call_combined_gvcf_2_JOB_ID: merge_and_call_combined_gvcf.call.AllSample
#-------------------------------------------------------------------------------
JOB_NAME=merge_and_call_combined_gvcf.call.AllSample
JOB_DEPENDENCIES=$merge_and_call_combined_gvcf_1_JOB_ID
JOB_DONE=job_output/merge_and_call_combined_gvcf/merge_and_call_combined_gvcf.call.AllSample.f265c60109fb147abe805494ca040afd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'merge_and_call_combined_gvcf.call.AllSample.f265c60109fb147abe805494ca040afd.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/GenomeAnalysisTK/3.8 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=1048576 -Xmx30G -jar $GATK_JAR \
  --analysis_type GenotypeGVCFs --useNewAFCalculator -G StandardAnnotation -G StandardHCAnnotation \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --reference_sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --variant variants/allSamples.hc.g.vcf.gz \
  --out variants/allSamples.hc.vcf.gz
merge_and_call_combined_gvcf.call.AllSample.f265c60109fb147abe805494ca040afd.mugqic.done
chmod 755 $COMMAND

merge_and_call_combined_gvcf_2_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=120:00:00 --mem=60G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$merge_and_call_combined_gvcf_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: variant_recalibrator
#-------------------------------------------------------------------------------
STEP=variant_recalibrator

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: variant_recalibrator_1_JOB_ID: variant_recalibrator.tranch.allSamples
#-------------------------------------------------------------------------------
JOB_NAME=variant_recalibrator.tranch.allSamples
JOB_DEPENDENCIES=$merge_and_call_combined_gvcf_2_JOB_ID
JOB_DONE=job_output/variant_recalibrator/variant_recalibrator.tranch.allSamples.7389658b08ae71d38ce03be974fba31d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'variant_recalibrator.tranch.allSamples.7389658b08ae71d38ce03be974fba31d.mugqic.done' > $COMMAND
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
variant_recalibrator.tranch.allSamples.7389658b08ae71d38ce03be974fba31d.mugqic.done
chmod 755 $COMMAND

variant_recalibrator_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=30G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$variant_recalibrator_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: variant_recalibrator_2_JOB_ID: variant_recalibrator.apply.allSamples
#-------------------------------------------------------------------------------
JOB_NAME=variant_recalibrator.apply.allSamples
JOB_DEPENDENCIES=$merge_and_call_combined_gvcf_2_JOB_ID:$variant_recalibrator_1_JOB_ID
JOB_DONE=job_output/variant_recalibrator/variant_recalibrator.apply.allSamples.1782dbbf24bd911bd923e5d8b7e03bdb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'variant_recalibrator.apply.allSamples.1782dbbf24bd911bd923e5d8b7e03bdb.mugqic.done' > $COMMAND
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
variant_recalibrator.apply.allSamples.1782dbbf24bd911bd923e5d8b7e03bdb.mugqic.done
chmod 755 $COMMAND

variant_recalibrator_2_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=30G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$variant_recalibrator_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: haplotype_caller_decompose_and_normalize
#-------------------------------------------------------------------------------
STEP=haplotype_caller_decompose_and_normalize

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: haplotype_caller_decompose_and_normalize_1_JOB_ID: decompose_and_normalize
#-------------------------------------------------------------------------------
JOB_NAME=decompose_and_normalize
JOB_DEPENDENCIES=$variant_recalibrator_2_JOB_ID
JOB_DONE=job_output/haplotype_caller_decompose_and_normalize/decompose_and_normalize.e76f00146b826b41f6b086fa627e7a5e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'decompose_and_normalize.e76f00146b826b41f6b086fa627e7a5e.mugqic.done' > $COMMAND
module purge && \
module load mugqic/htslib/1.8 mugqic/vt/0.57 && \
zless variants/allSamples.hc.vqsr.vcf | sed 's/ID=AD,Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa -  \
          | \
bgzip -cf \
 > \
variants/allSamples.hc.vqsr.vt.vcf.gz && tabix -pvcf variants/allSamples.hc.vqsr.vt.vcf.gz \
        
decompose_and_normalize.e76f00146b826b41f6b086fa627e7a5e.mugqic.done
chmod 755 $COMMAND

haplotype_caller_decompose_and_normalize_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu=4000M -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$haplotype_caller_decompose_and_normalize_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: haplotype_caller_flag_mappability
#-------------------------------------------------------------------------------
STEP=haplotype_caller_flag_mappability

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: haplotype_caller_flag_mappability_1_JOB_ID: haplotype_caller_flag_mappability
#-------------------------------------------------------------------------------
JOB_NAME=haplotype_caller_flag_mappability
JOB_DEPENDENCIES=$haplotype_caller_decompose_and_normalize_1_JOB_ID
JOB_DONE=job_output/haplotype_caller_flag_mappability/haplotype_caller_flag_mappability.e7c9e7b4e3b15848394797bf799e4319.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'haplotype_caller_flag_mappability.e7c9e7b4e3b15848394797bf799e4319.mugqic.done' > $COMMAND
module purge && \
module load mugqic/vcftools/0.1.14 mugqic/tabix/0.2.6 mugqic/htslib/1.8 && \
vcf-annotate \
  -d key=INFO,ID=MIL,Number=1,Type=String,Description='Mappability annotation. 300IS 40SD 1SHI. HC = too high coverage (>400), LC = too low coverage (<50), MQ = too low mean mapQ (<20), ND = no data at the position' \
  -c CHROM,FROM,TO,INFO/MIL \
  -a /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/mappabilityGC/Illu_PE.exclusion.bed.gz \
  variants/allSamples.hc.vqsr.vt.vcf.gz | \
bgzip -cf \
 > \
variants/allSamples.hc.vqsr.vt.mil.vcf.gz && tabix -pvcf variants/allSamples.hc.vqsr.vt.mil.vcf.gz \
        
haplotype_caller_flag_mappability.e7c9e7b4e3b15848394797bf799e4319.mugqic.done
chmod 755 $COMMAND

haplotype_caller_flag_mappability_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu=4000M -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$haplotype_caller_flag_mappability_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: haplotype_caller_snp_id_annotation
#-------------------------------------------------------------------------------
STEP=haplotype_caller_snp_id_annotation

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: haplotype_caller_snp_id_annotation_1_JOB_ID: haplotype_caller_snp_id_annotation
#-------------------------------------------------------------------------------
JOB_NAME=haplotype_caller_snp_id_annotation
JOB_DEPENDENCIES=$haplotype_caller_flag_mappability_1_JOB_ID
JOB_DONE=job_output/haplotype_caller_snp_id_annotation/haplotype_caller_snp_id_annotation.0d887f52b31f371bb02274783472a7ef.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'haplotype_caller_snp_id_annotation.0d887f52b31f371bb02274783472a7ef.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/snpEff/4.3 mugqic/htslib/1.8 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=2 -Xmx8G -jar $SNPEFF_HOME/SnpSift.jar annotate \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.dbSNP142.vcf.gz \
  variants/allSamples.hc.vqsr.vt.mil.vcf.gz | \
bgzip -cf \
 > \
variants/allSamples.hc.vqsr.vt.mil.snpId.vcf.gz && tabix -pvcf variants/allSamples.hc.vqsr.vt.mil.snpId.vcf.gz \
        
haplotype_caller_snp_id_annotation.0d887f52b31f371bb02274783472a7ef.mugqic.done
chmod 755 $COMMAND

haplotype_caller_snp_id_annotation_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$haplotype_caller_snp_id_annotation_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: haplotype_caller_snp_effect
#-------------------------------------------------------------------------------
STEP=haplotype_caller_snp_effect

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: haplotype_caller_snp_effect_1_JOB_ID: haplotype_caller_snp_effect
#-------------------------------------------------------------------------------
JOB_NAME=haplotype_caller_snp_effect
JOB_DEPENDENCIES=$haplotype_caller_snp_id_annotation_1_JOB_ID
JOB_DONE=job_output/haplotype_caller_snp_effect/haplotype_caller_snp_effect.931bd170ae5ce7ea1952853c85b705ae.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'haplotype_caller_snp_effect.931bd170ae5ce7ea1952853c85b705ae.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/snpEff/4.3 mugqic/htslib/1.8 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Xmx8G -jar $SNPEFF_HOME/snpEff.jar eff -lof \
   \
  -c $SNPEFF_HOME/snpEff.config \
  -i vcf \
  -o vcf \
  -csvStats variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.stats.csv \
  -stats variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.stats.html \
  hg19 \
  variants/allSamples.hc.vqsr.vt.mil.snpId.vcf.gz > variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf && \
bgzip -cf \
 \
 variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf > \
variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.gz && tabix -pvcf variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.gz \
        
haplotype_caller_snp_effect.931bd170ae5ce7ea1952853c85b705ae.mugqic.done
chmod 755 $COMMAND

haplotype_caller_snp_effect_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=16G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$haplotype_caller_snp_effect_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: haplotype_caller_dbnsfp_annotation
#-------------------------------------------------------------------------------
STEP=haplotype_caller_dbnsfp_annotation

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: haplotype_caller_dbnsfp_annotation_1_JOB_ID: dbnsfp_annotation
#-------------------------------------------------------------------------------
JOB_NAME=dbnsfp_annotation
JOB_DEPENDENCIES=$haplotype_caller_snp_effect_1_JOB_ID
JOB_DONE=job_output/haplotype_caller_dbnsfp_annotation/dbnsfp_annotation.62954a15ae6d49eb7bd66a177385e920.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'dbnsfp_annotation.62954a15ae6d49eb7bd66a177385e920.mugqic.done' > $COMMAND
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/snpEff/4.3 mugqic/htslib/1.8 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=2 -Xmx24G -jar $SNPEFF_HOME/SnpSift.jar dbnsfp \
  -v -db /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/annotations/dbNSFPv3.5c/dbNSFP3.5c.txt.gz \
  variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.gz \
  > variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf && \
bgzip -cf \
 \
 variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf > \
variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf.gz && tabix -pvcf variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf.gz \
        
dbnsfp_annotation.62954a15ae6d49eb7bd66a177385e920.mugqic.done
chmod 755 $COMMAND

haplotype_caller_dbnsfp_annotation_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem=40G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$haplotype_caller_dbnsfp_annotation_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: haplotype_caller_gemini_annotations
#-------------------------------------------------------------------------------
STEP=haplotype_caller_gemini_annotations

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: haplotype_caller_gemini_annotations_1_JOB_ID: gemini_annotations
#-------------------------------------------------------------------------------
JOB_NAME=gemini_annotations
JOB_DEPENDENCIES=$haplotype_caller_dbnsfp_annotation_1_JOB_ID
JOB_DONE=job_output/haplotype_caller_gemini_annotations/gemini_annotations.587eb6e4332dacf3c29d2022460bf135.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'gemini_annotations.587eb6e4332dacf3c29d2022460bf135.mugqic.done' > $COMMAND
module purge && \
module load mugqic/gemini/0.20.1 mugqic/htslib/1.8 && \
gemini load -v variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf.gz \
  -t snpEff --cores ${SLURM_CPUS_ON_NODE} --save-info-string \
  --tempdir ${SLURM_TMPDIR} \
  variants/allSamples.gemini.db
gemini_annotations.587eb6e4332dacf3c29d2022460bf135.mugqic.done
chmod 755 $COMMAND

haplotype_caller_gemini_annotations_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:00 --mem-per-cpu=4000M -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$haplotype_caller_gemini_annotations_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: haplotype_caller_metrics_vcf_stats
#-------------------------------------------------------------------------------
STEP=haplotype_caller_metrics_vcf_stats

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: haplotype_caller_metrics_vcf_stats_1_JOB_ID: haplotype_caller_metrics_change_rate
#-------------------------------------------------------------------------------
JOB_NAME=haplotype_caller_metrics_change_rate
JOB_DEPENDENCIES=$haplotype_caller_snp_effect_1_JOB_ID:$haplotype_caller_dbnsfp_annotation_1_JOB_ID
JOB_DONE=job_output/haplotype_caller_metrics_vcf_stats/haplotype_caller_metrics_change_rate.02e98e415a434a917d935066aaa7840f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'haplotype_caller_metrics_change_rate.02e98e415a434a917d935066aaa7840f.mugqic.done' > $COMMAND
module purge && \
module load mugqic/python/2.7.14 mugqic/mugqic_tools/2.2.2 && \
python $PYTHON_TOOLS/vcfStats.py \
  -v variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf.gz \
  -d /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.dict \
  -o variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf.part_changeRate.tsv \
  -f variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.stats.csv
haplotype_caller_metrics_change_rate.02e98e415a434a917d935066aaa7840f.mugqic.done
chmod 755 $COMMAND

haplotype_caller_metrics_vcf_stats_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu=4000M -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$haplotype_caller_metrics_vcf_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: run_multiqc
#-------------------------------------------------------------------------------
STEP=run_multiqc

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: run_multiqc_1_JOB_ID: multiqc_all_samples
#-------------------------------------------------------------------------------
JOB_NAME=multiqc_all_samples
JOB_DEPENDENCIES=$metrics_dna_picard_metrics_1_JOB_ID:$metrics_dna_picard_metrics_2_JOB_ID:$metrics_dna_picard_metrics_3_JOB_ID:$metrics_dna_sample_qualimap_1_JOB_ID:$metrics_dna_sambamba_flagstat_1_JOB_ID:$metrics_dna_fastqc_1_JOB_ID
JOB_DONE=job_output/run_multiqc/multiqc_all_samples.64fbd56bf58a8dc3deb3ee9d0605cb24.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'multiqc_all_samples.64fbd56bf58a8dc3deb3ee9d0605cb24.mugqic.done' > $COMMAND
module purge && \
module load mugqic/MultiQC/v1.6 && \
multiqc -f  \
 \
  metrics/dna/S1 \
-n metrics/dna/multiqc_report
multiqc_all_samples.64fbd56bf58a8dc3deb3ee9d0605cb24.mugqic.done
chmod 755 $COMMAND

run_multiqc_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu=4000M -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$run_multiqc_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: cram_output
#-------------------------------------------------------------------------------
STEP=cram_output

mkdir -p $JOB_OUTPUT_DIR/$STEP

#-------------------------------------------------------------------------------
# JOB: cram_output_1_JOB_ID: cram_output.S1
#-------------------------------------------------------------------------------
JOB_NAME=cram_output.S1
JOB_DEPENDENCIES=$recalibration_2_JOB_ID
JOB_DONE=job_output/cram_output/cram_output.S1.3d6dda598b8e9bde177b0f975a7ad843.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'cram_output.S1.3d6dda598b8e9bde177b0f975a7ad843.mugqic.done' > $COMMAND
module purge && \
module load mugqic/samtools/1.4.1 && \
samtools view -h -T $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa -C \
  alignment/S1/S1.sorted.dup.recal.bam \
  > alignment/S1/S1.sorted.dup.recal.cram
cram_output.S1.3d6dda598b8e9bde177b0f975a7ad843.mugqic.done
chmod 755 $COMMAND

cram_output_1_JOB_ID=$(echo "rm -f $JOB_DONE &&   $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem-per-cpu=4000M -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cram_output_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
LOG_MD5=$(echo $USER-'10.70.78.25-DnaSeq-S1.S1_R1' | md5sum | awk '{ print $1 }')
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=blg8425.int.ets1.calculquebec.ca&ip=10.70.78.25&pipeline=DnaSeq&steps=picard_sam_to_fastq,sym_link_fastq,trimmomatic,merge_trimmomatic_stats,skewer_trimming,bwa_mem_picard_sort_sam,sambamba_merge_sam_files,gatk_indel_realigner,sambamba_merge_realigned,fix_mate_by_coordinate,picard_mark_duplicates,recalibration,sym_link_final_bam,metrics_dna_picard_metrics,metrics_dna_sample_qualimap,metrics_dna_sambamba_flagstat,metrics_dna_fastqc,picard_calculate_hs_metrics,gatk_callable_loci,extract_common_snp_freq,baf_plot,gatk_haplotype_caller,merge_and_call_individual_gvcf,combine_gvcf,merge_and_call_combined_gvcf,variant_recalibrator,haplotype_caller_decompose_and_normalize,haplotype_caller_flag_mappability,haplotype_caller_snp_id_annotation,haplotype_caller_snp_effect,haplotype_caller_dbnsfp_annotation,haplotype_caller_gemini_annotations,haplotype_caller_metrics_vcf_stats,run_multiqc,cram_output&samples=1&md5=$LOG_MD5" --quiet --output-document=/dev/null

