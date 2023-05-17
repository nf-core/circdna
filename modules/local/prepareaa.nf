process PREPAREAA {
    tag "$meta.id"
    label 'process_low'

    //conda "conda-forge::python=3.6.15 anaconda::intervaltree=3.0.2 conda-forge::numpy=1.19.5 bioconda::pysam=0.15.3"
    //conda "conda-forge::python=2.7 bioconda::pysam=0.15.2 anaconda::flask=1.1.2 anaconda::cython=0.29.14 anaconda::numpy=1.16.6 anaconda::scipy=1.2.1 conda-forge::matplotlib=2.2.5 mosek::mosek=8.0.60 anaconda::future=0.18.2 anaconda::intervaltree=3.0.2"


    // Python3
    // Python2
    //conda "conda-forge::python=2.7 bioconda::pysam=0.20.0 anaconda::flask=1.1.2 anaconda::numpy=1.16.5 conda-forge::matplotlib=2.2.4 mosek::mosek=8.0.60 anaconda::intervaltree=3.0.2 anaconda::scipy=1.2.0 anaconda::future=0.18.2"
    conda "conda-forge::python=3.7 bioconda::pysam=0.16.0 anaconda::flask=2.2.2 conda-forge::numpy=1.21.6 conda-forge::matplotlib=3.2.2 anaconda::scipy=1.7.3 conda-forge::intervaltree=3.0.2 anaconda::future=0.18.2 mosek::mosek=9.0.88"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.gala3hhhhhhhhhorg/singularity/mulled-v2-6eefa51f13933d65b4f8155ca2f8cd81dea474ba:baa777f7c4e89a2ec4d1eab7d424a1f46503bac7-0':
        'quay.io/biocontainers/mulled-v2-6eefa51f13933d65b4f8155ca2f8cd81dea474ba:baa777f7c4e89a2ec4d1eab7d424a1f46503bac7-0' }"


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
    export AA_SRC=${projectDir}/bin
    REF=${params.reference_build}

    PrepareAA.py \\
        $args \\
        -s $prefix \\
        -t $task.cpus \\
        --cnv_bed $cns \\
        --sorted_bam $bam \\
        --cngain $cngain \\
        --ref $ref

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prepareaa: \$(echo \$(PrepareAA.py --version) | sed 's/^.*PrepareAA version //')
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
