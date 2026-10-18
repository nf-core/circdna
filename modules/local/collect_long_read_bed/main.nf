process COLLECT_LONG_READ_BED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    // One entry per engine for this sample, grouped by meta.id upstream.
    tuple val(meta), path(beds)

    output:
    tuple val(meta), path("${prefix}.eccdna_candidates.bed"), emit: candidates
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    collect_long_read_bed.py \\
        --out ${prefix}.eccdna_candidates.bed \\
        ${beds}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collect_long_read_bed: 1.0.0
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.eccdna_candidates.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        collect_long_read_bed: 1.0.0
    END_VERSIONS
    """
}
