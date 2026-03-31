process DEPTH_ANALYSIS {
    tag "$sample_id"
    label 'paraphase' // Reuse existing env with samtools+mosdepth
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path bed
    path dark_genes_bed

    output:
    path "${sample_id}_target_coverage.txt", emit: cov
    tuple val(sample_id), path("${sample_id}_target_coverage.txt"), emit: cov_tuple // For dependency joining
    path "${sample_id}_qc_metrics.txt", emit: qc
    path "${sample_id}_intron_depth.txt", emit: intron_report

    script:
    """
    # --- Micromamba (paraphase 컨테이너 + NXF_DOCKER_TASK_USER → HOME 필수)
    export TMPDIR=\$PWD
    export HOME=\$PWD
    export PIP_CACHE_DIR=\$PWD/.cache/pip
    export XDG_CACHE_HOME=\$PWD/.cache/xdg
    export PATH=\$PWD:\$PATH
    export CONDA_PKGS_DIRS=\$PWD/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    
    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX \$PIP_CACHE_DIR \$XDG_CACHE_HOME

    wget -q --no-check-certificate https://curl.se/ca/cacert.pem || true
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 \
            && chmod +x micromamba_bin || true
    fi
    
    if [ -f micromamba_bin ] && [ ! -x ./env/bin/samtools ]; then
        ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge samtools=1.16.1 mosdepth=0.3.3 -y
    fi
    [ -x ./env/bin/samtools ] && [ -x ./env/bin/mosdepth ] || { echo "ERROR: samtools/mosdepth env missing"; exit 1; }

    ST=\$PWD/env/bin/samtools
    MD=\$PWD/env/bin/mosdepth
    
    # --- Analysis ---

    if [ ! -f "${bam}.bai" ]; then
        ln -s ${bai} ${bam}.bai
    fi

    \$ST bedcov ${bed} ${bam} > ${sample_id}_target_coverage.txt

    \$MD --by ${bed} --thresholds 20,50,100 ${sample_id} ${bam}

    # 3. Parse Summary for User
    echo "Sample: ${sample_id}" > ${sample_id}_qc_metrics.txt
    echo "QC Metrics (on targets)" >> ${sample_id}_qc_metrics.txt
    
    MEAN_DEPTH=\$(grep "total_region" ${sample_id}.mosdepth.summary.txt | awk '{print \$4}')
    echo "Mean Depth: \$MEAN_DEPTH" >> ${sample_id}_qc_metrics.txt

    PCT_20X=\$(grep "^total.*	20	" ${sample_id}.mosdepth.region.dist.txt | awk '{printf "%.2f%%", \$3 * 100}')
    PCT_50X=\$(grep "^total.*	50	" ${sample_id}.mosdepth.region.dist.txt | awk '{printf "%.2f%%", \$3 * 100}')
    PCT_100X=\$(grep "^total.*	100	" ${sample_id}.mosdepth.region.dist.txt | awk '{printf "%.2f%%", \$3 * 100}')

    echo "% Bases > 20x: \$PCT_20X" >> ${sample_id}_qc_metrics.txt
    echo "% Bases > 50x: \$PCT_50X" >> ${sample_id}_qc_metrics.txt
    echo "% Bases > 100x: \$PCT_100X" >> ${sample_id}_qc_metrics.txt
    
    
    # --- 4. Intron Depth Verification (Merged) ---
    echo "Region\\tMean_Depth\\tMax_Depth\\tPercent_Covered" > ${sample_id}_intron_depth.txt
    
    while read -r chrom start end name; do
        \$ST depth -r "\$chrom:\$start-\$end" -a ${bam} | \\
        awk -v name="\$name" '{sum+=\$3; if(\$3>0) covered++; if(\$3>max) max=\$3; count++} END {if (count>0) print name "\t" sum/count "\t" max "\t" (covered/count)*100; else print name "\t0\t0\t0"}' >> ${sample_id}_intron_depth.txt
    done < ${dark_genes_bed}

    # Same windows as HBA1_HBA2_Region + SMN1_SMN2_Region rows above (dark_genes_plus.bed) — mean depth, not a separate hardcoded locus/median.
    echo "\\n[Dark-gene SMA / HBA depth QC]" >> ${sample_id}_intron_depth.txt
    ALPHA_MEAN=\$(awk -F'\\t' '\$1=="HBA1_HBA2_Region" {print \$2+0}' ${sample_id}_intron_depth.txt)
    SMA_MEAN=\$(awk -F'\\t' '\$1=="SMN1_SMN2_Region" {print \$2+0}' ${sample_id}_intron_depth.txt)
    echo "HBA1_HBA2_Region mean depth (alpha): \$ALPHA_MEAN" >> ${sample_id}_intron_depth.txt
    echo "SMN1_SMN2_Region mean depth (SMA/SMN): \$SMA_MEAN" >> ${sample_id}_intron_depth.txt

    WARN_MIN=${params.depth_warning_min_dark}
    INC_SMN=${params.depth_warning_include_smn}
    if [ "\$INC_SMN" != "true" ]; then
        echo "Note: depth QC threshold uses HBA only (depth_warning_include_smn=false); SMN/SMA rely on Paraphase/SMAca." >> ${sample_id}_intron_depth.txt
    fi
    if [ "\$WARN_MIN" = "0" ] || [ "\$WARN_MIN" = "0.0" ]; then
        echo "Note: depth QC warning disabled (depth_warning_min_dark=0)." >> ${sample_id}_intron_depth.txt
    else
        if [ "\$INC_SMN" = "true" ]; then
            TRIGGER=\$(awk -v a="\$ALPHA_MEAN" -v s="\$SMA_MEAN" -v t="\$WARN_MIN" 'BEGIN{print ((a+0 < t+0) || (s+0 < t+0)) ? 1 : 0}')
        else
            TRIGGER=\$(awk -v a="\$ALPHA_MEAN" -v t="\$WARN_MIN" 'BEGIN{print (a+0 < t+0) ? 1 : 0}')
        fi
        if [ "\$TRIGGER" = "1" ]; then
            if [ "\$INC_SMN" = "true" ]; then
                FAIL_WHY=\$(awk -v a="\$ALPHA_MEAN" -v s="\$SMA_MEAN" -v t="\$WARN_MIN" 'BEGIN{
                    fa=(a+0<t); fs=(s+0<t);
                    if (fa && fs) print "Both regions fail.";
                    else if (fa) print "HBA1_HBA2_Region fails only (SMN1_SMN2_Region is OK).";
                    else if (fs) print "SMN1_SMN2_Region fails only (HBA1_HBA2_Region is OK).";
                    else print "";
                }')
                echo "WARNING: Dark-gene mean depth below \$WARN_MIN x — \$FAIL_WHY Values: HBA1_HBA2_Region=\$ALPHA_MEAN, SMN1_SMN2_Region=\$SMA_MEAN (same rows as Region table; not mosdepth panel-wide mean)." >> ${sample_id}_intron_depth.txt
                echo "WARNING: SMN window often looks much lower than exome-wide depth (SMN1/SMN2 homology); use Paraphase/SMAca for SMA context." >> ${sample_id}_intron_depth.txt
            else
                echo "WARNING: HBA1_HBA2_Region mean depth below \$WARN_MIN x (value=\$ALPHA_MEAN). SMN depth not used for this gate (depth_warning_include_smn=false)." >> ${sample_id}_intron_depth.txt
            fi
        fi
    fi
    """
}
