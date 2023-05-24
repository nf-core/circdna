process AMPLICONCLASSIFIER_AMPLICONSIMILARITY {
    tag "AA Amplicons"
    label 'process_low'

    conda "bioconda::ampliconclassifier=0.4.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampliconclassifier:0.4.14--hdfd78af_0':
        'quay.io/biocontainers/ampliconclassifier:0.4.14--hdfd78af_0' }"

    input:
    path(input)

    output:
    path("*_scores.tsv")   , emit: scores
    path("*")
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    REF=${params.reference_build}
    export AA_DATA_REPO=${params.aa_data_repo}
    export AA_SRC=${projectDir}/bin

    amplicon_similarity.py \\
        --ref \$REF \\
        $args \\
        --input $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    REF=${params.reference_build}
    export AA_DATA_REPO=${params.aa_data_repo}
    export AA_SRC=${projectDir}/bin

    amplicon_similarity.py --help
    touch "ampliconclassifier_similarity_scores.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
    END_VERSIONS
    """
}
