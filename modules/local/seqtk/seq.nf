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
    tuple val(meta), path("*unicycler.circular.fastq.gz"), emit: fastq_circular, optional: true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqtk \\
        seq \\
        $args \\
        -F "#" \\
        $fasta > \\
        ${prefix}.unicycler.fastq

    grep -A3 "circular=true" ${prefix}.unicycler.fastq | \\
        grep -v "^--" | \\
        gzip --no-name > \\
        ${prefix}.unicycler.circular.fastq.gz

    gzip -n ${prefix}.unicycler.fastq
    """
}
