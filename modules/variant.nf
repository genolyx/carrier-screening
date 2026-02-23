process CALL_VARIANTS {
    tag "$sample_id"
    label 'gatk'
    publishDir "${params.outdir}/variant", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    path ref_dict
    path backbone_bed

    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), path("${sample_id}_filtered.vcf.gz.tbi"), emit: vcf

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
    
    ./micromamba_bin create -r \$MAMBA_ROOT_PREFIX -p ./env -c bioconda -c conda-forge gatk4=4.4.0.0 -y
    export PATH=\$PWD/env/bin:\$PATH

    gatk HaplotypeCaller \\
        -R ${ref_fasta} \\
        -I ${bam} \\
        -L ${backbone_bed} \\
        --interval-padding 100 \\
        -O ${sample_id}_raw.vcf.gz \\
        -ERC NONE \\
        --create-output-variant-index true
        
    gatk VariantFiltration \\
        -R ${ref_fasta} \\
        -V ${sample_id}_raw.vcf.gz \\
        -O ${sample_id}_filtered.vcf.gz \\
        --filter-expression "QD < 1.5" --filter-name "QD1.5" \\
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \\
        --filter-expression "SOR > 4.0" --filter-name "SOR4" \\
        --filter-expression "FS > 80.0" --filter-name "FS80" \\
        --filter-expression "MQ < 30.0" --filter-name "MQ30" \\
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
    """
}
