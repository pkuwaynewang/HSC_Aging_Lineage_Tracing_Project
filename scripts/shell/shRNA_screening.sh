#!/bin/bash

# shRNA Screening Library Analysis Pipeline
# Usage: bash script.sh R1.fastq.gz R2.fastq.gz sample_name reference_index

set -e  # Exit on error

# Check arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <R1.fastq.gz> <R2.fastq.gz> <sample_name> <reference_index>"
    echo "Example: $0 sample_54_S19_R1_001.fastq.gz sample_54_S19_R2_001.fastq.gz 54 ../reference/shRNA_screening_1st"
    exit 1
fi

# Input parameters
R1_FASTQ_GZ=$1
R2_FASTQ_GZ=$2
SAMPLE_NAME=$3
REFERENCE_INDEX=$4

# Create output directory structure
BASE_DIR="analysis_${SAMPLE_NAME}"
mkdir -p "${BASE_DIR}"/{unzipped,UMI,sam,bam,merged,final_result,without_UMI}

echo "=== Starting shRNA Screening Analysis Pipeline ==="
echo "Sample: ${SAMPLE_NAME}"
echo "R1: ${R1_FASTQ_GZ}"
echo "R2: ${R2_FASTQ_GZ}"
echo "Reference: ${REFERENCE_INDEX}"
echo ""

# ============================================
# STEP 1: Unzip FASTQ.GZ files
# ============================================
echo "[Step 1] Unzipping FASTQ.GZ files..."
R1_FASTQ="${BASE_DIR}/unzipped/$(basename ${R1_FASTQ_GZ} .gz)"
R2_FASTQ="${BASE_DIR}/unzipped/$(basename ${R2_FASTQ_GZ} .gz)"

gunzip -c "${R1_FASTQ_GZ}" > "${R1_FASTQ}"
gunzip -c "${R2_FASTQ_GZ}" > "${R2_FASTQ}"
echo "✓ Unzipping complete"
echo ""

# ============================================
# STEP 2: Extract UMI sequences
# ============================================
echo "[Step 2] Extracting UMI sequences from R1..."
UMI_FILE="${BASE_DIR}/UMI/${SAMPLE_NAME}.UMI.txt"

paste <(awk 'NR % 4 == 1 {sub(/ .*/, ""); print}' "${R1_FASTQ}") \
      <(awk 'NR % 4 == 2 {print substr($0, 1, 10)}' "${R1_FASTQ}") \
      > "${UMI_FILE}"
echo "✓ UMI extraction complete: ${UMI_FILE}"
echo ""

# ============================================
# STEP 3: Clean UMI file
# ============================================
echo "[Step 3] Cleaning UMI file..."
CLEAN_UMI="${BASE_DIR}/UMI/${SAMPLE_NAME}.clean.UMI.txt"

sed 's/^@//' "${UMI_FILE}" | awk -F " " '{print $1"\t"$2}' > "${CLEAN_UMI}"
echo "✓ UMI cleaning complete: ${CLEAN_UMI}"
echo ""

# ============================================
# STEP 4: Bowtie alignment
# ============================================
echo "[Step 4] Running Bowtie alignment..."
SAM_FILE="${BASE_DIR}/sam/${SAMPLE_NAME}.sam"

bowtie -v 0 --trim5 25 -m 1 "${REFERENCE_INDEX}" \
       -1 "${R1_FASTQ}" \
       -2 "${R2_FASTQ}" \
       -S "${SAM_FILE}"
echo "✓ Bowtie alignment complete: ${SAM_FILE}"
echo ""

# ============================================
# STEP 5: Convert SAM to BAM and extract fields
# ============================================
echo "[Step 5] Converting SAM to BAM and extracting fields..."
BAM_FILE="${BASE_DIR}/bam/${SAMPLE_NAME}.bam"
READY_TO_MERGE="${BASE_DIR}/bam/${SAMPLE_NAME}.ready_to_merge.txt"

samtools view -bS "${SAM_FILE}" > "${BAM_FILE}"
samtools view "${BAM_FILE}" | awk '$2 == 99 || $2 == 147 || $2 == 83 || $2 == 163 {print $1"\t"$2"\t"$3}'  > "${READY_TO_MERGE}"
echo "✓ BAM conversion and extraction complete"
echo ""

# ============================================
# STEP 6: Merge BAM data with UMI
# ============================================
echo "[Step 6] Merging alignment data with UMI..."
MERGED_FILE="${BASE_DIR}/merged/${SAMPLE_NAME}.merged.txt"

join -t $'\t' \
     <(sort -k1,1 "${READY_TO_MERGE}") \
     <(sort -k1,1 "${CLEAN_UMI}") \
     > "${MERGED_FILE}"
echo "✓ Merging complete: ${MERGED_FILE}"
echo ""

# ============================================
# STEP 7: Generate final results WITH UMI
# ============================================
echo "[Step 7] Generating final results with UMI deduplication..."
FINAL_WITH_UMI="${BASE_DIR}/final_result/${SAMPLE_NAME}.txt"

awk -F'\t' '{a[$3][$4]=1} END {for (group in a) printf "%s\t%d\n", group, length(a[group])}' \
    "${MERGED_FILE}" | sort > "${FINAL_WITH_UMI}"
echo "✓ Final results (with UMI): ${FINAL_WITH_UMI}"
echo ""

# ============================================
# STEP 8: Generate final results WITHOUT UMI
# ============================================
echo "[Step 8] Generating final results without UMI deduplication..."
FINAL_WITHOUT_UMI="${BASE_DIR}/without_UMI/${SAMPLE_NAME}.txt"

cut -f 3 "${READY_TO_MERGE}" | sort | uniq -c > "${FINAL_WITHOUT_UMI}"
echo "✓ Final results (without UMI): ${FINAL_WITHOUT_UMI}"
echo ""

# ============================================
# Summary
# ============================================
echo "=== Analysis Complete ==="
echo "Output directory: ${BASE_DIR}/"
echo ""
echo "Key output files:"
echo "  - UMI file: ${CLEAN_UMI}"
echo "  - SAM file: ${SAM_FILE}"
echo "  - BAM file: ${BAM_FILE}"
echo "  - Merged file: ${MERGED_FILE}"
echo "  - Final (with UMI): ${FINAL_WITH_UMI}"
echo "  - Final (without UMI): ${FINAL_WITHOUT_UMI}"
echo ""
echo "Summary statistics:"
echo "  Total reads in R1: $(wc -l < "${R1_FASTQ}" | awk '{print $1/4}')"
echo "  UMIs extracted: $(wc -l < "${CLEAN_UMI}")"
echo "  Aligned reads: $(samtools view -c "${BAM_FILE}")"
echo "  Unique shRNAs (with UMI): $(wc -l < "${FINAL_WITH_UMI}")"
echo "  Unique shRNAs (without UMI): $(wc -l < "${FINAL_WITHOUT_UMI}")"
