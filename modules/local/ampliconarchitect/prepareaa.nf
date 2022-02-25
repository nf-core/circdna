// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process AMPLICONARCHITECT_PREPAREAA {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::python=2.7 anaconda::numpy=1.15.4 conda-forge::matplotlib=2.2.5 conda-forge:intervaltree=3.0.2 bioconda::pysam=0.17.0 mosek::mosek=8.0.60 anaconda::scipy=1.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(bam), path(bai), path(cns)

    output:
    tuple val(meta), path("*CNV_SEEDS.bed"), emit: bed
    path "versions.yml"          , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def cn_gain_threshold = params.aa_cn_gain_threshold ? "--cngain $params.aa_cn_gain_threshold" : ""
    def no_filter   = params.aa_no_filter ? "--no_filter" : ""

    """
    AA_DATA_REPO=${params.aa_data_repo}
    MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}
    REF=${params.reference_build}

    PrepareAA.py \\
        -s ${meta.id} \\
        -t ${task.cpus} \\
        $options.args \\
        --sorted_bam $bam \\
        --ref \$REF \\
        --cnv_bed $cns \\
        $cn_gain_threshold \\
        $no_filter
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python /' )
    END_VERSIONS
    """
}
