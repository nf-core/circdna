process AMPLICONARCHITECT_AMPLICONCLASSIFIER {
    tag "$meta.id"
    label 'process_medium'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "conda-forge::python=3.8 conda-forge::matplotlib=2.2.5 conda-forge::intervaltree" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(cycles), path(graph)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*"), emit: all
    path "*.classifier_stdout.log", emit: log
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    make_AmpliconClassifier_input.sh ./ ${meta.id}.AmpliconClassifier
    REF=${params.reference_build}
    AA_DATA_REPO=${params.aa_data_repo}

    amplicon_classifier.py \\
        --ref \$REF \\
        $args \\
        --input ${meta.id}.AmpliconClassifier.input \\
        > ${meta.id}.classifier_stdout.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampliconarchitect:
    END_VERSIONS
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python /' )
    END_VERSIONS
    """
}
