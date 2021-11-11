// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SAMBLASTER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::samblaster=0.1.26 bioconda::samtools=1.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samblaster:0.1.26--h7d875b9_1"
    } else {
        container "quay.io/biocontainers/samblaster:0.1.26--h7d875b9_1"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.split.bam"), emit: split_bam

    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    samtools \\
        view \\
        -h \\
        -@ $task.cpus \\
        $bam |
    samblaster \\
        $options.args \\
        -s ${prefix}.split.sam \\
        > /dev/null

    samtools view \\
        -@ $task.cpus \\
        -o ${prefix}.split.bam \\
        -bS ${prefix}.split.sam 

    rm ${prefix}.split.sam 

    echo \$(samblaster --version 2>&1) > ${software}.version.txt
    """
}
