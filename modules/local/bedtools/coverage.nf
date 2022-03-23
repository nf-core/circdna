process BEDTOOLS_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_2"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2"
    }

    input:
    tuple val(meta), path(bed), path(bam)

    output:
    tuple val(meta), path("*.filterd.bed"), emit: bed
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools \\
        coverage \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        -T $prefix \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version)
    END_VERSIONS
    """
}
