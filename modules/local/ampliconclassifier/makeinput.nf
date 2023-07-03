process AMPLICONCLASSIFIER_MAKEINPUT {
    tag 'AA Amplicons'
    label 'process_low'

    conda "bioconda::ampliconclassifier=0.4.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampliconclassifier:0.4.14--hdfd78af_0':
        'quay.io/biocontainers/ampliconclassifier:0.4.14--hdfd78af_0' }"

    input:
    val(id)
    path(summary)

    output:
    path "*.input"      , emit: input
    // path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # Take vectors as input
    vector1=(\$(echo $id | sed 's/\\[//g' | sed 's/, / /g' | sed 's/\\]//g' ))
    vector2=(\$(echo $summary | sed 's/\\[//g' | sed 's/, /,/g' ))

    echo \$vector1
    echo \$vector2

    # Check that vectors are of equal length
    if [ \${#vector1[@]} -ne \${#vector2[@]} ]; then
        echo "Vectors are not of equal length."
        exit 1
    fi

    # Sort the vectors
    vector1_sorted=(\$(printf '%s\n' "\${vector1[@]}"|sort))
    vector2_sorted=(\$(printf '%s\n' "\${vector2[@]}"|sort))

    # Write to file
    for index in \${!vector1_sorted[@]}; do
        echo \${vector1_sorted[\$index]}\t\${vector2_sorted[\$index]}
    done > run_metadata_list.input

#    "${task.process}":
#        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
#    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch "ampliconclassifier.input"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
    END_VERSIONS
    """
}
