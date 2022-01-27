// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
options        = initOptions(params.options)

process CIRCLEMAP_REALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::circle-map=1.1.4 conda-forge::biopython=1.77" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/circle-map:1.1.4--pyh864c0ab_1"
    } else {
        container "quay.io/biocontainers/circle-map:1.1.4--pyh864c0ab_1"
    }

    input:
    tuple val(meta), path(re_bam), path(re_bai), path(qname), path(sbam), path(sbai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    circle_map.py \\
        Realign \\
        $options.args \\
        -i $re_bam \\
        -qbam $qname \\
        -sbam $sbam \\
        -fasta $fasta \\
        --threads $task.cpus \\
        --cmapq $params.circdna_filter_mapq \\
        --bases $params.coverage_bases \\
        --extension $params.coverage_extension \\
        --split $params.circdna_filter_nSplit \\
        -o ${prefix}_circularDNA_coordinates.bed

    echo \$(Circle-Map --help 2<&1) | grep -o "version=[0-9].[0-9].[0-9]" > ${software}.version.txt
    """
}
