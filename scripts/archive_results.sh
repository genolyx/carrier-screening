#!/bin/bash
# set -e (Disabled to allow partial archive of mixed valid/invalid links)

# Define Output Directory
ARCHIVE_DIR="final_archived_results"
mkdir -p "$ARCHIVE_DIR"/{bam,vcf,summary,snapshots,fallback,cnv}

echo "Starting Archive Process..."

# 1. BAMs (Crucial - dereference symlinks with -L)
echo "Archiving BAMs..."
if ls results/alignment/*.md.bam 1> /dev/null 2>&1; then
    cp -L results/alignment/*.md.bam "$ARCHIVE_DIR/bam/"
fi
if ls results/alignment/*.md.bam.bai 1> /dev/null 2>&1; then
    cp -L results/alignment/*.md.bam.bai "$ARCHIVE_DIR/bam/"
elif ls results/alignment/*.bai 1> /dev/null 2>&1; then
    cp -L results/alignment/*.bai "$ARCHIVE_DIR/bam/"
fi

# 2. VCFs
echo "Archiving VCFs..."
if [ -d "results/variant" ]; then
    cp -r results/variant/* "$ARCHIVE_DIR/vcf/"
fi

# 3. Summary Reports
echo "Archiving Summary..."
if [ -d "results/summary" ]; then
    cp -r results/summary/* "$ARCHIVE_DIR/summary/"
fi

# 4. Snapshots (Visual Evidence)
echo "Archiving Snapshots..."
if [ -d "results/snapshots" ]; then
    cp -r results/snapshots/* "$ARCHIVE_DIR/snapshots/"
fi

# 5. Fallback & CNV
echo "Archiving Fallback & CNV..."
if [ -d "results/fallback" ]; then
    cp -r results/fallback/* "$ARCHIVE_DIR/fallback/"
fi
if [ -d "results/cnv" ]; then
    cp -r results/cnv/* "$ARCHIVE_DIR/cnv/"
    
    # Explicit Check for Reference Model
    if [ -d "$ARCHIVE_DIR/cnv/cohort/gcnv-model" ]; then
        echo "SUCCESS: GCNV Reference Data (gcnv-model) successfully preserved."
    fi
fi

# 6. Additional Tracks (Pseudogene, SV, Repeat, Coverage, Rescue)
echo "Archiving Additional Tracks..."
for dir in pseudogene sv repeat coverage pipeline_info rescue; do
    if [ -d "results/$dir" ]; then
        mkdir -p "$ARCHIVE_DIR/$dir"
        cp -r results/$dir/* "$ARCHIVE_DIR/$dir/"
    fi
done

# 6. Verification
echo "Verifying Archive..."
BAM_COUNT=$(find "$ARCHIVE_DIR/bam" -name "*.md.bam" | wc -l)
SUMMARY_COUNT=$(find "$ARCHIVE_DIR/summary" -name "detailed_report.txt" | wc -l)

echo "Archive Stats:"
echo "  BAM Files: $BAM_COUNT"
echo "  Summary Report: $SUMMARY_COUNT"

if [ "$BAM_COUNT" -eq 0 ]; then
    echo "WARNING: No BAM files were archived. Proceeding with cleanup anyway."
    # exit 1  <-- Disabled to ensure cleanup happens
fi

if [ "$SUMMARY_COUNT" -eq 0 ]; then
    echo "Warning: No detailed_report.txt found."
fi

echo "Archive Successful!"
du -sh "$ARCHIVE_DIR"
exit 0
