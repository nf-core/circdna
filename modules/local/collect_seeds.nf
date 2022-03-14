process COLLECT_SEEDS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/python:3.9--1' :
            'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(cns)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cngain = params.aa_cngain
    """
    collect_seeds.py \\
        --sample $prefix \\
        --cns $cns \\
        --cngain $cngain

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
