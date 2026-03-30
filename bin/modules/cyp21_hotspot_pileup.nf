// BAM pileup at the seven NM_000500.9 sites in cyp21a2_hotspots.tsv (cross-check vs SNV VCF).
// Imports cyp21_paralog_pileup.py from the same task work dir.

process CYP21_HOTSPOT_PILEUP {
    tag "$sample_id"
    publishDir "${params.outdir}/fallback", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(hotspots_tsv), path(hotspot_py), path(paralog_pileup_py)

    output:
    path "${sample_id}_cyp21_hotspot_pileup.tsv", emit: tsv

    script:
    """
    export TMPDIR=\$PWD
    export HOME=\$PWD
    export PATH=\$PWD:\$PATH

    # Stage helper next to this script so "import cyp21_paralog_pileup" resolves
    cp ${paralog_pileup_py} cyp21_paralog_pileup.py

    python3 ${hotspot_py} \\
        --bam ${bam} \\
        --sites ${hotspots_tsv} \\
        -o ${sample_id}_cyp21_hotspot_pileup.tsv
    """
}
