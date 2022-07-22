process SEQTK_SEQ {
    tag "$meta.id"
    label 'process_low'
    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--hed695b0_2':
        'quay.io/biocontainers/seqtk:1.3--hed695b0_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*unicycler.fastq.gz"), emit: fastq
    path "versions.yml"                         , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqtk \\
        seq \\
        $args \\
        -F "#" \\
        $fasta | \\
        gzip -c > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
