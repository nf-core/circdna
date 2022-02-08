// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process AMPLICONARCHITECT_PREPAREAA {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::python=2.7 anaconda::numpy=1.15.4 conda-forge::matplotlib=2.2.5 conda-forge:intervaltree=3.1.0 bioconda::pysam=0.16.0.1 mosek::mosek=9.2.38 anaconda::scipy=1.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "quay.io/biocontainers/YOUR-TOOL-HERE"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta), path(cns)

    output:
    tuple val(meta), path("*CNV_SEEDS.bed"), emit: bed

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
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
        --cnv_bed $cns
    """
}
