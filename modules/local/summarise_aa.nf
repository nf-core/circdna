process SUMMARISE_AA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/python:3.9--1' :
            'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(summary_file), path(class_file)

    output:
    tuple val(meta), path(tsv), emit: tsv
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<-END_VERSIONS > versions.yml
    summarise_aa.py \\
        --summary $summary_file \\
        --class_file $class_file \\
        --id ${meta.id} \\
        --output ${prefix}.aa_results_summary.tsv
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
