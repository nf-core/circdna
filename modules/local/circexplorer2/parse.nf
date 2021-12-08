// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process CIRCEXPLORER2_PARSE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::circexplorer2=2.3.8" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/circexplorer2:2.3.8--pyh864c0ab_1"
    } else {
        container "quay.io/biocontainers/circexplorer2:2.3.8--pyh864c0ab_1"
    }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bed")  , emit: bed
    path "*.log"                    , emit: log
    path "*.version.txt"            , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    CIRCexplorer2 parse $options.args $bam -b ${prefix}.temp.bed > ${prefix}_CIRCexplorer2_parse.log
    cat ${prefix}.temp.bed | tr "/" "\t" > ${prefix}.circexplorer_circdna.bed
    rm ${prefix}.temp.bed
    
    echo \$(CIRCexplorer2 --version) > ${software}.version.txt
    """
}
