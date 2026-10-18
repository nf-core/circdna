process CIDERSEQ_DECONCAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/bioinfortools/ciderseq2:2.0.0'

    input:
    tuple val(meta), path(reads)
    path(config)

    output:
    tuple val(meta), path("${prefix}.fa"),   emit: deconcat
    tuple val(meta), path("${prefix}.stat"), emit: stat
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: (meta.genome ? "${meta.id}.${meta.genome}" : "${meta.id}")
    """
    ciderseq_deconcat.py \\
        --config ${config} \\
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
    prefix = task.ext.prefix ?: (meta.genome ? "${meta.id}.${meta.genome}" : "${meta.id}")
    """
    touch ${prefix}.fa
    touch ${prefix}.stat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciderseq2: 2.0.0
    END_VERSIONS
    """
}
