process AMPLICONARCHITECT_AMPLICONCLASSIFIER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ampliconclassifier=0.4.5=hdfd78af_1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampliconclassifier:0.4.5--hdfd78af_1':
        'quay.io/biocontainers/ampliconclassifier:0.4.5--hdfd78af_1' }"

    input:
    tuple val(meta), path(cycles), path(graph)

    output:
    tuple val(meta), path("*_profiles.tsv")     , emit: class_tsv   , optional: true
    tuple val(meta), path("*gene_list.tsv")     , emit: gene_list   , optional: true
    tuple val(meta), path("*ecDNA_counts.tsv")  , emit: ecDNA_counts, optional: true
    path "*.classifier_stdout.log"              , emit: log         , optional: true
    path "versions.yml"                         , emit: versions    , optional: true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    REF=${params.reference_build}
    export AA_DATA_REPO=${params.aa_data_repo}
    export AA_SRC=${projectDir}/bin

    # Make AmpliconClassifier Input from graph and cycles files
    make_AmpliconClassifier_input.sh ./ ${meta.id}.AmpliconClassifier

    amplicon_classifier.py \\
        --ref \$REF \\
        $args \\
        --input ${meta.id}.AmpliconClassifier.input \\
        > ${meta.id}.classifier_stdout.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export AA_DATA_REPO=${params.aa_data_repo}
    export MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}
    export AA_SRC=${projectDir}/bin
    REF=${params.reference_build}

    touch "${prefix}.amplicon_classification_profiles.tsv"
    touch "${prefix}.classifier_stdout.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
    END_VERSIONS
    """
}
