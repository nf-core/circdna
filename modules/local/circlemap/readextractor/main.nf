process CIRCLEMAP_READEXTRACTOR {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::circle-map=1.1.4=pyh5e36f6f_2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circle-map:1.1.4--pyh5e36f6f_2':
        'quay.io/biocontainers/circle-map:1.1.4--pyh5e36f6f_2' }"

    input:
    tuple val(meta), path(qname_bam)

    output:
    tuple val(meta), path("*.circular_read_candidates.bam"), emit: bam
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    circle_map.py \\
        ReadExtractor -i $qname_bam \\
        -o ${prefix}.circular_read_candidates.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Circle-Map: \$(echo \$(circle_map.py --help 2<&1 | grep -o "version=[0-9].[0-9].[0-9]" | sed 's/version=//g'))
    END_VERSIONS
    """
}
