process AMPLICONCLASSIFIER_MAKEINPUT {
    tag 'AA Amplicons'
    label 'process_low'

    conda "bioconda::ampliconclassifier=0.4.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampliconclassifier:0.4.14--hdfd78af_0':
        'quay.io/biocontainers/ampliconclassifier:0.4.14--hdfd78af_0' }"

    input:
    path(graph)
    path(cycles)

    output:
    path "*.input"      , emit: input
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    make_input.sh ./ ampliconclassifier

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    make_input.py --help
    touch "ampliconclassifier.input"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
    END_VERSIONS
    """
}
