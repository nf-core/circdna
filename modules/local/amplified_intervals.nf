process AMPLIFIED_INTERVALS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=2.7 conda-forge::matplotlib=2.2.5 anaconda::numpy=1.15.4 bioconda::pysam=0.17.0 mosek::mosek=8.0.60 anaconda::scipy=1.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(bed), path(bam), path(bai)

    output:
    tuple val(meta), path("*CNV_SEEDS.bed"), emit: bed
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cngain = params.aa_cngain
    def ref = params.reference_build
    """
    export AA_DATA_REPO=${params.aa_data_repo}
    export MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}
    REF=${params.reference_build}

    amplified_intervals.py \\
        $args \\
        --bed $bed \\
        --out ${prefix}_AA_CNV_SEEDS \\
        --bam $bam \\
        --gain $cngain \\
        --ref $ref

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: echo \$(python --version 2<&1 | sed 's/Python //g')
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

    touch ${prefix}_AA_CNV_SEEDS.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: echo \$(python --version 2<&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
