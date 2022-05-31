process SUMMARISE_AA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "pandas=1.1.5" : null)
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
            'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(summary_file), path(class_file)

    output:
    tuple val(meta), path("*aa_results_summary.tsv"), emit: txt
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    summarise_aa.py \\
        --summary $summary_file \\
        --class_file $class_file \\
        --id ${meta.id} \\
        --output ${prefix}.aa_results_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.aa_results_summary.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
