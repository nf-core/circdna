process AMPLICONARCHITECT_AMPLICONSIMILARITY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::intervaltree=3.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(cycles), path(graph)

    output:
    tuple val(meta), path("*_scores.tsv")   , emit: scores
    path("*")
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    REF=${params.reference_build}
    export AA_DATA_REPO=${params.aa_data_repo}
    export AA_SRC=${projectDir}/bin

    make_AmpliconClassifier_input.sh ./ ${meta.id}.AmpliconSimilarity

    amplicon_similarity.py \\
        --ref \$REF \\
        $args \\
        --input ${meta.id}.AmpliconSimilarity.input \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    REF=${params.reference_build}
    export AA_DATA_REPO=${params.aa_data_repo}
    export AA_SRC=${projectDir}/bin

    touch "${prefix}.similarity_scores.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
