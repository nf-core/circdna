process GETCIRCULARREADS {
    tag "$meta.id"
    label 'process_low'
    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.unicycler.circular.fastq.gz"), optional: true, emit: fastq
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    if gzip -t ${fastq} 2>/dev/null; then
        ZCAT_CMD="gzip -cd"
    else
        ZCAT_CMD="cat"
    fi

    \$ZCAT_CMD ${fastq} > temp.fastq
    if grep -q "circular=true" temp.fastq; then
        cat temp.fastq | grep -A3 "circular=true" | \\
            grep -v "^--" | \\
            gzip --no-name > \\
            ${prefix}.unicycler.circular.fastq.gz
    else
        touch ${prefix}.unicycler.circular.fastq.gz
    fi
    rm temp.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > '${prefix}.unicycler.circular.fastq.gz'
    touch versions.yml
    """
}
