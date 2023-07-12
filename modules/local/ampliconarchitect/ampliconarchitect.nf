process AMPLICONARCHITECT_AMPLICONARCHITECT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ampliconsuite=0.1555.2 mosek::mosek=10.1b1"
    container 'quay.io/nf-core/prepareaa:1.0.0'

    input:
    tuple val(meta), path(bam), path(bai), path(bed)

    output:
    tuple val(meta), path("*cycles.txt")    , optional: true, emit: cycles
    tuple val(meta), path("*graph.txt")     , optional: true, emit: graph
    tuple val(meta), path("*cnseg.txt")     , optional: true, emit: cnseg
    tuple val(meta), path("*.out")          , optional: true, emit: out
    tuple val(meta), path("*.{pdf,png}")    , optional: true, emit: svview
    tuple val(meta), path("*_summary.txt")  , optional: true, emit: summary
    tuple val(meta), path("*{log.txt,flag.txt}")            , emit: log
    tuple val(meta), path("*sample_metadata.json")          , emit: s_json
    tuple val(meta), path("*run_metadata.json")             , emit: r_json
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export AA_DATA_REPO=${params.aa_data_repo}
    export MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}

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
        -t $task.cpus \\
        --bam $bam \\
        --bed $bed \\
        --ref \$REF \\
        -s "${prefix}" \\
        --run_AA
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconSuite-pipeline.py: \$(echo \$(AmpliconSuite-pipeline.py --version) | sed 's/^.*PrepareAA version //')
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

    touch "${prefix}.logs.txt"
    touch "${prefix}.cycles.txt"
    touch "${prefix}.graph.txt"
    touch "${prefix}.out"
    touch "${prefix}_cnseg.txt"
    touch "${prefix}.pdf"
    touch "${prefix}.png"
    touch "${prefix}_summary.txt"

    AmpliconArchitect.py --help

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconSuite-pipeline.py: \$(echo \$(AmpliconSuite-pipeline.py --version) | sed 's/^.*PrepareAA version //')
    END_VERSIONS
    """
}
