process SAMBLASTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samblaster=0.1.26 bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:ba4a02b56f3e524a6e006bcd99fe8cc1d7fe09eb-0' :
        'quay.io/biocontainers/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:ba4a02b56f3e524a6e006bcd99fe8cc1d7fe09eb-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.split.bam"), emit: split_bam
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

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
        samblaster: \$(echo \$(samblaster --version 2>&1) | sed 's/samblaster: Version //g')
    END_VERSIONS
    """
}
