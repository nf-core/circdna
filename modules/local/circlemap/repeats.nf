// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process CIRCLEMAP_REPEATS {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::circle-map=1.1.4 conda-forge::biopython=1.77" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-f7a0df9f53b95b416b41af4bb3c1e214667f88c0"
    } else {
        container "quay.io/biocontainers/circle-map:1.1.4--pyh864c0ab_1"
    }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    circle_map.py \\
        Repeats \\
        $options.args \\
        -i $bam \\
        -o ${prefix}_circularDNA_repeats_coordinates.bed

    echo \$(Circle-Map --help 2<&1) | grep -o "version=[0-9].[0-9].[0-9]" > ${software}.version.txt
    """
}
