#!/bin/bash

# Migrate fastq_old files to new structure
# Extract sample name from filename pattern

echo "Migrating fastq_old to new directory structure..."

WORK_DIR="2601"  # Current month
SOURCE_DIR="fastq_old"
TARGET_BASE="fastq/$WORK_DIR"

mkdir -p "$TARGET_BASE"

# Group files by sample name
for f in "$SOURCE_DIR"/*.fq.gz "$SOURCE_DIR"/*.fastq.gz; do
    [ -e "$f" ] || continue
    
    filename=$(basename "$f")
    
    # Extract sample name (remove _R1/_R2 and extensions)
    sample_name=$(echo "$filename" | sed -E 's/_R[12]_.*//; s/_[12]\.fq.*//; s/\.fq.*//; s/\.fastq.*//')
    
    # Create sample directory
    sample_dir="$TARGET_BASE/$sample_name"
    mkdir -p "$sample_dir"
    
    # Copy file
    cp "$f" "$sample_dir/"
    
    echo "  Copied: $filename -> $sample_dir/"
done

# Mark all migrated samples as completed
for sample_dir in "$TARGET_BASE"/*; do
    [ -d "$sample_dir" ] || continue
    
    # Check if R1/R2 pair exists
    r1_count=$(ls "$sample_dir"/*R1* 2>/dev/null | wc -l)
    r2_count=$(ls "$sample_dir"/*R2* 2>/dev/null | wc -l)
    
    if [ "$r1_count" -gt 0 ] && [ "$r2_count" -gt 0 ]; then
        # Mark as completed since they came from fastq_old
        touch "$sample_dir/analysis.completed"
        echo "  Marked $(basename $sample_dir) as completed"
    fi
done

echo "Migration complete!"
