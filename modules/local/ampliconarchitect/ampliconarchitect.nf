// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process AMPLICONARCHITECT_AMPLICONARCHITECT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    //conda (params.enable_conda ? "conda-forge::python=2.7 bioconda::pysam=0.16.0.1 flask=1.1.2 cython=0.29.22 numpy=1.15.4 scipy=1.1.3 conda-forge::matplotlib=2.2.5" : null)
    conda (params.enable_conda ? "conda-forge::python=2.7 bioconda::pysam=0.17.0 flask=1.1.2 cython=0.29.15 numpy=1.16.5 scipy=1.2.1 conda-forge::matplotlib=2.2.5 mosek::mosek=8.0.60" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "quay.io/biocontainers/YOUR-TOOL-HERE"
    }

    input:
    tuple val(meta), path(bam), path(bai), path(bed)

    output:
    tuple val(meta), path("*"), emit: bam
    path "*.version.txt"          , emit: version
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*cycles.txt"), optional: true, emit: cycles
    tuple val(meta), path("*graph.txt"), optional: true, emit: graph

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    AA_DATA_REPO=${params.aa_data_repo}
    MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}
    # output=${params.outdir}/ampliconarchitect
    AmpliconArchitect.py $options.args \\
        --bam $bam --bed $bed --ref "GRCh38" --out "${prefix}"

    AmpliconArchitect.py --version > ${software}.version.txt
    """
}
