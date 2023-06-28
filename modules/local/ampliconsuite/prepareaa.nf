process PREPAREAA {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ampliconsuite=0.1555.2 mosek::mosek=10.1b1"
    container '/home/local/BICR/dschreye/src/AmpliconSuite-pipeline/docker/test/ampliconsuite.img'

    input:
    tuple val(meta), path(bam), path(cns)

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
    export AA_SRC=\$(dirname \$(readlink -f \$(which AmpliconArchitect.py)))
    export AC_SRC=\$(dirname \$(readlink -f \$(which amplicon_classifier.py)))
    REF=${params.reference_build}

    AmpliconSuite-pipeline.py \\
        $args \\
        -s $prefix \\
        -t $task.cpus \\
        --cnv_bed $cns \\
        --bam $bam \\
        --ref $ref \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconSuite-pipeline.py: \$(echo \$(AmpliconSuite-pipeline.py --version) | sed 's/^.*PrepareAA version //')
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
        prepareaa: \$(echo \$(PrepareAA.py --version) | sed 's/^.*PrepareAA version //')
    END_VERSIONS
    """
}
