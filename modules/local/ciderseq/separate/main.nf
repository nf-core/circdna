process CIDERSEQ_SEPARATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/bioinfortools/ciderseq2:2.0.0'

    input:
    tuple val(meta), path(reads)
    path(config)
    path(blastdb)

    output:
    tuple val(meta), path("${prefix}.*.fa"), emit: separated
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ciderseq_separate.py \\
        --config ${config} \\
        --blastdb ${blastdb} \\
        --outdir . \\
        --prefix ${prefix} \\
        $args \\
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciderseq2: 2.0.0
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nohit.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciderseq2: 2.0.0
    END_VERSIONS
    """
}
