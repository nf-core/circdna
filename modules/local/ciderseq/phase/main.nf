process CIDERSEQ_PHASE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/bioinfortools/ciderseq2:2.0.0'

    input:
    tuple val(meta), path(deconcat)
    path(annotation)
    path(config)
    val(genome)

    output:
    tuple val(meta), path("${prefix}.{fasta,gb,genbank,fa}"), emit: phased
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${genome}"
    """
    ciderseq_phase.py \\
        --config ${config} \\
        --genome ${genome} \\
        --outdir . \\
        --prefix ${prefix} \\
        $args \\
        ${deconcat} ${annotation}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciderseq2: 2.0.0
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.${genome}"
    """
    touch ${prefix}.fasta
    touch ${prefix}.gb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciderseq2: 2.0.0
    END_VERSIONS
    """
}
