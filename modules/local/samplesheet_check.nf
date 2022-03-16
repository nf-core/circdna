<<<<<<< HEAD
// Import generic module functions
process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/python:3.9--1' :
            'quay.io/biocontainers/python:3.9--1' }"
=======
process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"
>>>>>>> TEMPLATE

    input:
    path samplesheet

    output:
<<<<<<< HEAD
    path '*.csv'        , emit: csv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when
=======
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions
>>>>>>> TEMPLATE

    script: // This script is bundled with the pipeline, in nf-core/circdna/bin/
    """
    check_samplesheet.py \\
        $samplesheet \\
<<<<<<< HEAD
        samplesheet.valid.csv \\
        $params.input_format
=======
        samplesheet.valid.csv
>>>>>>> TEMPLATE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
