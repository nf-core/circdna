process GETCIRCULARREADS {
    tag "$meta.id"
    label 'process_low'
    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*unicycler.circular.fastq.gz"), optional: true, emit: fastq

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    zcat $fastq > temp.fastq
    if grep -q "circular=true" temp.fastq; then
        cat temp.fastq | grep -A3 "circular=true" | \\
            grep -v "^--" | \\
            gzip --no-name > \\
            ${prefix}.unicycler.circular.fastq.gz
    fi
    rm temp.fastq
    """
}
