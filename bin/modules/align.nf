process INDEX_BWA {
    tag "bwa_index"
    label 'bwa'
    storeDir "refs/bwa_index"
    input:
    path ref_fasta
    output:
    tuple path(ref_fasta), path("${ref_fasta}.*"), emit: indices
    script:
    """
    # --- Dependency Setup (Micromamba) ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export HOME=\$PWD
    export XDG_CACHE_HOME=\$PWD/.cache
    export CONDA_PKGS_DIRS=\$HOME/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    mkdir -p \$XDG_CACHE_HOME
    mkdir -p \$CONDA_PKGS_DIRS
    mkdir -p \$MAMBA_ROOT_PREFIX

    # SSL Fix
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    # Install Micromamba & BWA
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge bwa=0.7.17 samtools=1.16.1 -y
    export PATH=\$PWD/env/bin:\$PATH

    bwa index ${ref_fasta}
    """
}

process ALIGN_AND_SORT {
    tag "$sample_id"
    label 'bwa'
    publishDir "${params.outdir}/alignment", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path ref_fasta
    path ref_fai
    path ref_bwa_indices
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam
    
    script:
    """
    # --- Dependency Setup (Micromamba) ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export CONDA_PKGS_DIRS=/home/sam/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba
    
    mkdir -p \$CONDA_PKGS_DIRS \$MAMBA_ROOT_PREFIX

    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    # Install BWA and Samtools
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge bwa=0.7.17 samtools=1.16.1 -y
    export PATH=\$PWD/env/bin:\$PATH

    # Pipe BWA directly to Samtools Sort
    # Optimization: Reduces I/O by avoiding large intermediate SAM file
    bwa mem -t ${task.cpus} -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' ${ref_fasta} ${reads[0]} ${reads[1]} | \\
    samtools sort -@ ${task.cpus} -m 768M -o ${sample_id}.bam -
    
    samtools index ${sample_id}.bam
    """
}

process MARK_DUPLICATES {
    tag "$sample_id"
    label 'gatk'
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.md.bam"), path("${sample_id}.md.bai"), emit: bam
    path "${sample_id}.duplicate_metrics.txt", emit: metrics

    script:
    """
    # --- Dependency Setup (Micromamba) ---
    export TMPDIR=\$PWD
    export PATH=\$PWD:\$PATH
    export HOME=\$PWD
    export XDG_CACHE_HOME=\$PWD/.cache
    export CONDA_PKGS_DIRS=\$HOME/.cache/micromamba_pkgs
    export MAMBA_ROOT_PREFIX=\$PWD/micromamba

    mkdir -p \$XDG_CACHE_HOME
    mkdir -p \$CONDA_PKGS_DIRS
    mkdir -p \$MAMBA_ROOT_PREFIX

    # SSL Fix
    wget -q --no-check-certificate https://curl.se/ca/cacert.pem
    export SSL_CERT_FILE=\$PWD/cacert.pem
    export MAMBA_SSL_VERIFY=false

    # Install Micromamba & GATK
    if [ ! -f "micromamba_bin" ]; then
        wget -qO micromamba_bin https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
        chmod +x micromamba_bin
    fi
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge gatk4=4.4.0.0 -y
    export PATH=\$PWD/env/bin:\$PATH

    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_id}.md.bam \\
        -M ${sample_id}.duplicate_metrics.txt \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY SILENT
    """
}
