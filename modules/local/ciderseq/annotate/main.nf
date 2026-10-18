process CIDERSEQ_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/bioinfortools/ciderseq2:2.0.0'

    input:
    tuple val(meta), path(reads)
    path(config)
    path(protein_db)

    output:
    tuple val(meta), path("${prefix}.json"), emit: annotation
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: (meta.genome ? "${meta.id}.${meta.genome}" : "${meta.id}")
    """
    ciderseq_annotate.py \\
        --config ${config} \\
        --protein_db ${protein_db} \\
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
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciderseq2: 2.0.0
    END_VERSIONS
    """
}
