process GETCIRCULARREADS {
    tag "$meta.id"
    label 'process_low'
    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*unicycler.circular.fastq.gz"), optional: true, emit: fastq
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    zcat $fastq > temp.fastq
    if grep -q "circular=true" temp.fastq; then
        cat temp.fastq | grep -A3 "circular=true" | \\
            grep -v "^--" | \\
            gzip --no-name > \\
            ${prefix}.fastq.gz
    fi
    rm temp.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
