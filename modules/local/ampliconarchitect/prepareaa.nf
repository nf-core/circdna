process AMPLICONARCHITECT_PREPAREAA {
    tag "$meta.id"
    label 'process_medium'

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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export AA_DATA_REPO=${params.aa_data_repo}
    export MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}
    export AA_SRC=${projectDir}/bin

    PrepareAA.py \\
        -s ${prefix} \\
        -t ${task.cpus} \\
        $args \\
        --sorted_bam $bam \\
        --ref $params.reference_build \\
        --cnv_bed $cns
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python /' )
    END_VERSIONS
    """
}
