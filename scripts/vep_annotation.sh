#!/bin/bash

set -euo pipefail


# -------------------------
#  Project directories
# -------------------------

PROJECT_DIR="$HOME/wes_chr22_project"
REF_DIR="${PROJECT_DIR}/references"
RESULTS_DIR="${PROJECT_DIR}/results"
LOG_DIR="${PROJECT_DIR}/logs"
TMP_DIR="${PROJECT_DIR}/tmp"


# -------------------------
#  Input File
# -------------------------
VCF_PASS="${RESULTS_DIR}/NA12878.chr22.filtered.pass.vcf.gz"
REF="${REF_DIR}/chr22.fasta"


# -------------------------
#  Output File
# -------------------------
ANN="${RESULTS_DIR}/NA12878.chr22.vep.vcf"
STATS="${RESULTS_DIR}/NA12878.chr22.vep.summary.html"



# -------------------------
#  STEP: Annotation
# -------------------------
echo "Starting VEP annotation...."
vep \
  --input_file "${PASS_VCF}" \
  --output_file "${ANN}" \
  --vcf \
  --cache \
  --offline \
  --assembly GRCh38 \
  --species homo_sapiens \
  --fasta "${REF}" \
  --symbol \
  --canonical \
  --biotype \
  --hgvs \
  --protein \
  --numbers \
  --variant_class \
  --af_gnomad \
  --sift b \
  --polyphen b \
  --fork 2 \
  --stats_file ${STATS} \
  > ${LOG_DIR}/vep.log 2>&1

echo "VEP annotation completed."
