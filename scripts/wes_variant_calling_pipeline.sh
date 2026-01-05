#!/usr/bin/env bash

###############################################################################
# WES chr22 Variant Calling Pipeline (GATK Best Practices)
# Sample: NA12878
# Environment: WSL Linux
###############################################################################

# ------------------------------- #
# Error Handling
# ------------------------------- #
set -euo pipefail

# ------------------------------- #
# Project directories
# ------------------------------- #
PROJECT_DIR="$HOME/wes_chr22_project"
FASTQ_DIR="${PROJECT_DIR}/fastq"
REF_DIR="${PROJECT_DIR}/references"
SCRIPT_DIR="${PROJECT_DIR}/scripts"
RESULTS_DIR="${PROJECT_DIR}/results"
LOG_DIR="${PROJECT_DIR}/logs"
TMP_DIR="${PROJECT_DIR}/tmp"

# Create directories
mkdir -p "${RESULTS_DIR}" "${LOG_DIR}" "${TMP_DIR}"

# ------------------------------- #
# Sample and region
# ------------------------------- #
SAMPLE="NA12878"
TARGET_CHR="chr22"

# ------------------------------- #
# Input FASTQ files
# ------------------------------- #
R1="${FASTQ_DIR}/SRR1518011_1.fastq.gz"
R2="${FASTQ_DIR}/SRR1518011_2.fastq.gz"

# ------------------------------- #
# Reference and known sites
# ------------------------------- #
REF="${REF_DIR}/chr22.fasta"
DBSNP="${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
MILLS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# ------------------------------- #
# Output files
# ------------------------------- #
ALIGNED_SAM="${RESULTS_DIR}/${SAMPLE}.chr22.sam"
ALIGNED_BAM="${RESULTS_DIR}/${SAMPLE}.chr22.bam"
BAM_SORTED="${RESULTS_DIR}/${SAMPLE}.chr22.sorted.bam"
RG_BAM="${RESULTS_DIR}/${SAMPLE}.chr22.sorted.rg.bam"
BAM_DEDUP="${RESULTS_DIR}/${SAMPLE}.chr22.dedup.bam"
RECAL_TABLE="${RESULTS_DIR}/${SAMPLE}.chr22.recal.table"
BAM_RECAL="${RESULTS_DIR}/${SAMPLE}.chr22.recal.bam"
VCF_RAW="${RESULTS_DIR}/${SAMPLE}.chr22.raw.vcf.gz"
VCF_FILTERED="${RESULTS_DIR}/${SAMPLE}.chr22.filtered.vcf.gz"
VCF_PASS="${RESULTS_DIR}/${SAMPLE}.chr22.filtered.pass.vcf.gz"

THREADS=2

###############################################################################
# Step 1: FastQC
###############################################################################
echo "Running FastQC..."
fastqc "${R1}" "${R2}" \
  --outdir "${RESULTS_DIR}" \
  > "${LOG_DIR}/fastqc.log" 2>&1

echo "Done!"


###############################################################################
# Step 2: BWA-MEM Alignment
###############################################################################
echo "Running BWA-MEM alignment..."

bwa mem -t "${THREADS}" "${REF}" "${R1}" "${R2}" > "${ALIGNED_SAM}"\
  2> "${LOG_DIR}/bwa.log" | \
samtools sort -@ "${THREADS}" -o "${ALIGNED_BAM}" \
  2> "${LOG_DIR}/samtools_sort.log"

samtools index "${BAM_SORTED}"

echo "Done!"


###############################################################################
# Step 3: Mark Duplicates
###############################################################################
echo "Adding Read Groups..."
gatk --java-options "-Xmx2g -Djava.io.tmpdir=${TMP_DIR}" AddOrReplaceReadGroups \
  -I "${BAM_SORTED}" \
  -O "${RG_BAM}" \
  -RGID NA12878 \
  -RGLB WES \
  -RGPL ILLUMINA \
  -RGPU unit1 \
  -RGSM NA12878 \
  > logs/add_rg.log 2>&1

samtools index "${RG_BAM}" \
  2> logs/add_rg_index.log

echo "Done!"

echo "Marking duplicates..."
gatk --java-options "-Xmx2g -Djava.io.tmpdir=${TMP_DIR}" MarkDuplicates \
  -I "${BAM_SORTED}" \
  -O "${BAM_DEDUP}" \
  -M "${RESULTS_DIR}/${SAMPLE}.dup.metrics.txt" \
  > "${LOG_DIR}/markduplicates.log" 2>&1

samtools index "${BAM_DEDUP}" \
  2> logs/dedup_index.log

echo "Done!"


###############################################################################
# Step 4: Base Quality Score Recalibration (BQSR)
###############################################################################
echo "Running BaseRecalibrator..."
gatk --java-options "-Xmx2g -Djava.io.tmpdir=${TMP_DIR}" BaseRecalibrator \
  -R "${REF}" \
  -I "${BAM_DEDUP}" \
  --known-sites "${DBSNP}" \
  --known-sites "${MILLS}" \
  -O "${RECAL_TABLE}" \
  > "${LOG_DIR}/baserecalibrator.log" 2>&1

echo "Applying BQSR..."
gatk --java-options "-Xmx2g -Djava.io.tmpdir=${TMP_DIR}" ApplyBQSR \
  -R "${REF}" \
  -I "${BAM_DEDUP}" \
  --bqsr-recal-file "${RECAL_TABLE}" \
  -O "${BAM_RECAL}" \
  > "${LOG_DIR}/applybqsr.log" 2>&1

samtools index "${BAM_RECAL}"

echo "Done!"


###############################################################################
# Step 5: Variant Calling (HaplotypeCaller)
###############################################################################
echo "Calling variants on ${TARGET_CHR}..."
gatk --java-options "-Xmx2g -Djava.io.tmpdir=${TMP_DIR}" HaplotypeCaller \
  -R "${REF}" \
  -I "${BAM_RECAL}" \
  -L "${TARGET_CHR}" \
  -O "${VCF_RAW}" \
  > "${LOG_DIR}/haplotypecaller.log" 2>&1

echo "Done!"


###############################################################################
# Step 6: Hard Variant Filtering
###############################################################################
echo "Filtering variants..."
gatk VariantFiltration \
  -R "${REF}" \
  -V "${VCF_RAW}" \
  -O "${VCF_FILTERED}" \
  --filter-name "QD_lt_2" --filter-expression "QD < 2.0" \
  --filter-name "FS_gt_60" --filter-expression "FS > 60.0" \
  --filter-name "MQ_lt_40" --filter-expression "MQ < 40.0" \
  > "${LOG_DIR}/variantfiltration.log" 2>&1

bcftools view -f PASS "${VCF_FILTERED}" \
  -Oz -o "${VCF_PASS}"
2> logs/bcftools_passvariants.log
bcftools index "${VCF_PASS}"

echo "Done!"


###############################################################################
echo "Pipeline completed successfully!"
echo "Final VCF: ${VCF_PASS}"
###############################################################################
