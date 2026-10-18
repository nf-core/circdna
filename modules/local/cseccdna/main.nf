process CS_ECCDNA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/bioinfortools/ciderseq2:2.0.0'

    input:
    tuple val(meta), path(input_dir)
    path(genome)

    output:
    tuple val(meta), path("*.fasta"), emit: candidates
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ciderseq_eccdna.py \\
        --genome ${genome} \\
        --input_dir ${input_dir} \\
        --outdir . \\
        --blast_threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciderseq2: 2.0.0
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.assess1.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciderseq2: 2.0.0
    END_VERSIONS
    """
}
