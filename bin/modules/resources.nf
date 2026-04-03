process PREPARE_VIZ_RESOURCES {
    tag "setup_resources"
    label 'local'
    publishDir "${params.outdir}/resources", mode: 'copy',
        saveAs: { filename -> filename == "viz_env" ? null : filename }
    errorStrategy 'retry'
    maxRetries 2

    input:
    path ref_fasta
    path ref_fai

    output:
    path "viz_env",              emit: env
    path "refGene.sorted.gtf.gz",     emit: gtf
    path "refGene.sorted.gtf.gz.tbi", emit: gtf_index
    path "smn_ref_wide.fa",      emit: smn_ref
    path "smn_ref_wide.fa.*",    emit: smn_index

    script:
    // viz 툴이 Docker 이미지에 pre-install 되어 있으므로 viz_env → /opt/conda 심볼릭 링크
    def gtf_provided = params.refgene_gtf ? true : false
    """
    export TMPDIR=\$PWD

    # --- 1. Shared Viz Environment ---
    # igv-reports, pysam, bwa, samtools, tabix are pre-installed in the Docker image (/opt/conda).
    # viz_env is a symlink so GENERATE_VISUAL_EVIDENCE can use: export PATH=\$PWD/viz_env/bin:\$PATH
    ln -sf /opt/conda viz_env
    export PATH=/opt/conda/bin:\$PATH

    # --- 2. RefSeq GTF ---
    if [ "${gtf_provided}" = "true" ] && [ -f "${params.refgene_gtf}" ]; then
        echo "Using pre-downloaded RefSeq GTF: ${params.refgene_gtf}"
        cp "${params.refgene_gtf}" refGene.gtf.gz
    else
        echo "Downloading RefSeq GTF from UCSC..."
        wget -q --no-check-certificate -O refGene.gtf.gz \\
            https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
    fi

    echo "Sorting and indexing GTF..."
    gzip -d -c refGene.gtf.gz | grep -v "^#" | sort -k1,1 -k4,4n | bgzip > refGene.sorted.gtf.gz
    tabix -p gff refGene.sorted.gtf.gz

    # --- 3. SMN Mini-Reference ---
    echo "Generating SMN Mini-Reference..."
    samtools faidx ${ref_fasta} chr5:70920000-70960000 > smn_ref_wide.fa
    sed -i "s/>.*/>chr5/" smn_ref_wide.fa
    bwa index smn_ref_wide.fa
    """
}
