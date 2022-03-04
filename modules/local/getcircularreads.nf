process GETCIRCULARREADS {
    tag "$meta.id"
    label 'process_low'
    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*unicycler.circular.fastq.gz")   , emit: fastq
    path "versions.yml"                                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    zcat $fastq | grep -A3 "circular=true" | \\
        grep -v "^--" | \\
        gzip --no-name > \\
        ${prefix}.unicycler.circular.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zcat: \$(zcat --version | grep zcat | sed 's/zcat (gzip) //g')
    END_VERSIONS
    """
}
