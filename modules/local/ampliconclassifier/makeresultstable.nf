process AMPLICONCLASSIFIER_MAKERESULTSTABLE {
    tag 'AA Amplicons'
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ampliconclassifier=0.4.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampliconclassifier:0.4.14--hdfd78af_0':
        'quay.io/biocontainers/ampliconclassifier:0.4.14--hdfd78af_0' }"

    input:
    path (input_file)
    path (class_file)
    path (gene_list)
    path (feature_entropy)
    path (basic_properties)
    path (bed_files)

    output:
    path "*result_data.json"    , emit: json
    path "*result_table.tsv"    , emit: tsv
    path "index.html"           , emit: html
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Create subdirectories in working directory
    mkdir ampliconclassifier_classification_bed_files
    mv $bed_files ampliconclassifier_classification_bed_files/

    make_results_table.py \\
        $args \\
        --input $input_file \\
        --classification_file $class_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch ampliconclasifier_result_data.json
    touch ampliconclasifier_result_table.tsv
    touch index.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
    END_VERSIONS
    """
}
