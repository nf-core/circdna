// Import generic module functions
process SAMBLASTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samblaster=0.1.26 bioconda::samtools=1.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samblaster:0.1.26--h7d875b9_1"
    } else {
        container "quay.io/biocontainers/samblaster:0.1.26--h7d875b9_1"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.split.bam"), emit: split_bam
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        view \\
        -h \\
        -@ $task.cpus \\
        $bam |
    samblaster \\
        $args \\
        -s ${prefix}.split.sam \\
        > /dev/null

    samtools view \\
        -@ $task.cpus \\
        -o ${prefix}.split.bam \\
        -bS ${prefix}.split.sam

    rm ${prefix}.split.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samblaster: \$(echo \$(samblaster --version 2>&1))
    END_VERSIONS
    """
}
