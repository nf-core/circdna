process PREPAREAA {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ampliconsuite=0.1555.2 mosek::mosek=10.1b1"
    container 'quay.io/nf-core/prepareaa:1.0.0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*CNV_SEEDS.bed") , emit: bed
    path "*.log"                            , emit: log
    path "*run_metadata.json"               , emit: run_metadata_json
    path "*sample_metadata.json"            , emit: sample_metadata_json
    path "*timing_log.txt"                  , emit: timing_log
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cngain = params.aa_cngain
    def ref = params.reference_build
    """
    export AA_DATA_REPO=${params.aa_data_repo}
    export MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}
    # Define Variables AA_SRC and AC_SRC
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
    REF=${params.reference_build}

    AmpliconSuite-pipeline.py \\
        $args \\
        -s $prefix \\
        -t $task.cpus \\
        --bam $bam \\
        --ref $ref \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconSuite-pipeline.py: \$(AmpliconSuite-pipeline.py --version | sed 's/AmpliconSuite-pipeline version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cngain = params.aa_cngain
    def ref = params.reference_build
    """
    export AA_DATA_REPO=${params.aa_data_repo}
    export MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}
    REF=${params.reference_build}

    touch "${prefix}_CNV_SEEDS.bed"
    touch "${prefix}.log"
    touch "${prefix}.run_metadata.json"
    touch "${prefix}.sample_metadata.json"
    touch "${prefix}.timing_log.txt"
    touch "${prefix}_summary.txt"

    PrepareAA.py --help

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconSuite-pipeline.py: \$(AmpliconSuite-pipeline.py --version | sed 's/AmpliconSuite-pipeline version //')
    END_VERSIONS
    """
}
