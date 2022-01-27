// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
options        = initOptions(params.options)

process CIRCLEMAP_READEXTRACTOR {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::circle-map=1.1.4 biopython=1.77" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/circle-map:1.1.4--pyh864c0ab_1"
    } else {
        container "quay.io/biocontainers/circle-map:1.1.4--pyh864c0ab_1"
    }

    input:
    tuple val(meta), path(qname_bam)

    output:
    tuple val(meta), path("*.circular_read_candidates.bam"), emit: circ_bam
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    circle_map.py ReadExtractor -i $qname_bam \\
        -o ${prefix}.circular_read_candidates.bam
    echo \$(Circle-Map --help 2<&1) | grep -o "version=[0-9].[0-9].[0-9]" > ${software}.version.txt
    """
}
