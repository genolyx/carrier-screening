#!/bin/bash
#
# Quick Analysis Script - Interactive Mode
# 대화형으로 샘플 분석 실행
#

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "======================================"
echo "Dark Gene Pipeline - Quick Analysis"
echo "======================================"
echo ""

# Work directory 입력
read -p "Work directory (YYMM, e.g., 2601) [default: $(date +%y%m)]: " WORK_DIR
WORK_DIR=${WORK_DIR:-$(date +%y%m)}

# 사용 가능한 샘플 목록 표시
DATA_DIR="${SCRIPT_DIR}/data"
FASTQ_BASE="${DATA_DIR}/fastq/${WORK_DIR}"

if [ -d "$FASTQ_BASE" ]; then
    echo ""
    echo "Available samples in ${WORK_DIR}:"
    ls -1 "$FASTQ_BASE" 2>/dev/null | while read sample; do
        if [ -d "${FASTQ_BASE}/${sample}" ]; then
            if [ -f "${FASTQ_BASE}/${sample}/analysis.completed" ]; then
                echo "  ✓ $sample (completed)"
            else
                echo "  • $sample"
            fi
        fi
    done
    echo ""
fi

# Sample name 입력
read -p "Sample name: " SAMPLE_NAME

if [ -z "$SAMPLE_NAME" ]; then
    echo "Error: Sample name is required"
    exit 1
fi

# Cleanup 옵션
read -p "Clean up work directory after completion? (y/N): " CLEANUP_CHOICE
if [[ "$CLEANUP_CHOICE" =~ ^[Yy]$ ]]; then
    CLEANUP_OPT="--cleanup"
else
    CLEANUP_OPT=""
fi

echo ""
echo "Starting analysis with:"
echo "  Work Directory: ${WORK_DIR}"
echo "  Sample: ${SAMPLE_NAME}"
echo "  Cleanup: ${CLEANUP_OPT:-disabled}"
echo ""

read -p "Continue? (Y/n): " CONFIRM
if [[ "$CONFIRM" =~ ^[Nn]$ ]]; then
    echo "Cancelled."
    exit 0
fi

# 분석 실행
exec "${SCRIPT_DIR}/run_analysis.sh" \
    --work-dir "$WORK_DIR" \
    --sample "$SAMPLE_NAME" \
    $CLEANUP_OPT
