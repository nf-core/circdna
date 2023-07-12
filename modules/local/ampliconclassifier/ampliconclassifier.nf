process AMPLICONCLASSIFIER_AMPLICONCLASSIFIER {
    tag "AA Amplicons"
    label 'process_low'

    conda "bioconda::ampliconsuite=0.1555.2 mosek::mosek=10.1b1"
    container 'quay.io/nf-core/prepareaa:1.0.0'

    input:
    path (graphs)
    path (cycles)
    path (cnseg)

    output:
    path ("*"               )   , emit: all         , optional: true
    path ("versions.yml"    )   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "ampliconarchitect"

    """
    export AA_DATA_REPO=${params.aa_data_repo}
    if ! command -v AmpliconArchitect.py &> /dev/null; then
        export AA_SRC=\$(dirname \$(python -c "import ampliconarchitectlib; print(ampliconarchitectlib.__file__)"))
    else
        export AA_SRC=\$(dirname \$(readlink -f \$(which AmpliconArchitect.py)))
    fi

    if ! command -v amplicon_classifier.py &> /dev/null; then
        export AC_SRC=\$(dirname \$(python -c "import ampliconclassifierlib; print(ampliconclassifierlib.__file__)"))
    else
        export AC_SRC=\$(dirname \$(readlink -f \$(which amplicon_classifier.py)))
    fi

    AmpliconSuite-pipeline.py \\
        -s $prefix \\
        --completed_AA_runs ./ \\
        -t $task.cpus \\
        --ref $params.reference_build

    mv ampliconarchitect_classification/* ./
    rmdir ampliconarchitect_classification

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconSuite-pipeline.py: \$(AmpliconSuite-pipeline.py --version | sed 's/AmpliconSuite-pipeline version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    export AA_DATA_REPO=${params.aa_data_repo}
    export MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}
    export AA_SRC=${projectDir}/bin
    REF=${params.reference_build}

    touch "ampliconclassifier_amplicon_classification_profiles.tsv"
    touch "ampliconclassifier_classifier_stdout.log"

    amplicon_classifier.py --help

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconSuite-pipeline.py: \$(AmpliconSuite-pipeline.py --version | sed 's/AmpliconSuite-pipeline version //')
    END_VERSIONS
    """
}
