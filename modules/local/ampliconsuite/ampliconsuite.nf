process AMPLICONSUITE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'nf-core/prepareaa:1.0.4'

    input:
    tuple val(meta), path(bam)
    path(mosek_license_dir)
    path(aa_data_repo)

    output:
    path "*.bed"                    , emit: bed
    path "*.cns"                    , emit: cns, optional: true
    path "*.cnr.gz"                 , emit: cnr, optional: true
    path "*.log"                    , emit: log
    path "*run_metadata.json"       , emit: run_metadata_json
    path "*sample_metadata.json"    , emit: sample_metadata_json
    path "*timing_log.txt"          , emit: timing_log
    path "*.input"                  , emit: ac_input, optional: true
    path "*logs.txt"                , emit: logs, optional: true
    path "*cycles.txt"              , emit: cycles, optional: true
    path "*graph.txt"               , emit: graph, optional: true
    path "*"
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cngain = params.aa_cngain
    def ref = params.reference_build
    """
    export AA_DATA_REPO=\$(echo $aa_data_repo)
    export MOSEKLM_LICENSE_FILE=\$(echo $mosek_license_dir)
    # Define Variables AA_SRC and AC_SRC
    export AA_SRC=\$(dirname \$(python -c "import ampliconarchitectlib; print(ampliconarchitectlib.__file__)"))
    export AC_SRC=\$(dirname \$(which amplicon_classifier.py))
    REF=${params.reference_build}

    AmpliconSuite-pipeline.py \\
        $args \\
        -s $prefix \\
        -t $task.cpus \\
        --bam $bam \\
        --ref $ref \\
        --run_AA --run_AC \\
        $args

    # Move Files to base work directory
    find ${prefix}_cnvkit_output/ -type f -print0 | xargs -0 mv -t ./
    find ${prefix}_AA_results/ -type f -print0 | xargs -0 mv -t ./
    find ${prefix}_classification/ -type f -print0 | xargs -0 mv -t ./

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
    export AA_DATA_REPO=\$(echo $aa_data_repo)
    export MOSEKLM_LICENSE_FILE=\$(echo $mosek_license_dir)
    # Define Variables AA_SRC and AC_SRC
    export AA_SRC=\$(dirname \$(python -c "import ampliconarchitectlib; print(ampliconarchitectlib.__file__)"))
    export AC_SRC=\$(dirname \$(which amplicon_classifier.py))
    REF=${params.reference_build}

    touch "${prefix}_CNV_SEEDS.bed"
    touch "${prefix}.log"
    touch "${prefix}.run_metadata.json"
    touch "${prefix}.sample_metadata.json"
    touch "${prefix}.timing_log.txt"
    touch "${prefix}_summary.txt"

    AmpliconSuite-pipeline.py --help

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconSuite-pipeline.py: \$(AmpliconSuite-pipeline.py --version | sed 's/AmpliconSuite-pipeline version //')
    END_VERSIONS
    """
}
