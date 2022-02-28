process CIRCLEMAP_REPEATS {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::circle-map=1.1.4 biopython=1.77" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circle-map:1.1.4--pyh864c0ab_1':
        'quay.io/biocontainers/circle-map:1.1.4--pyh864c0ab_1' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    circle_map.py \\
        Repeats \\
        $args \\
        -i $bam \\
        -o ${prefix}_circularDNA_repeats_coordinates.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Circle-Map: \$(echo \$(Circle-Map --help 2<&1) | grep -o "version=[0-9].[0-9].[0-9]")
    END_VERSIONS
    """
}
