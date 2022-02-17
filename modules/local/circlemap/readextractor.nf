process CIRCLEMAP_READEXTRACTOR {
    tag "$meta.id"
    label 'process_high'

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
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    circle_map.py ReadExtractor -i $qname_bam \\
        -o ${prefix}.circular_read_candidates.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Circle-Map: \$(echo \$(Circle-Map --help 2<&1) | grep -o "version=[0-9].[0-9].[0-9]")
    END_VERSIONS
    """
}
