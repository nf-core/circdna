process GETCIRCULARREADS {
    tag "$meta.id"
    label 'process_low'
    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*unicycler.circular.fastq.gz"), emit: fastq

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat $fastq | grep -A3 "circular=true" | \\
        grep -v "^--" | \\
        gzip --no-name > \\
        ${prefix}.unicycler.circular.fastq.gz
    """
}
