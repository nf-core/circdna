// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)


process SEQTK_SEQ {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--hed695b0_2':
        'quay.io/biocontainers/seqtk:1.3--hed695b0_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*unicycler.fastq.gz"), emit: fastq
    tuple val(meta), path("*unicycler.circular.fastq.gz"), emit: fastq_circular
    // tuple val(meta), path("*.unicycler.circular.fastq.gz"), emit: fastq, optional: true
    // TODO nf-core: List additional required output channels/values here

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqtk \\
        seq \\
        $args \\
        -F "#" \\
        $fasta > \\
        ${prefix}.unicycler.fastq

    grep -A3 "circular=true" ${prefix}.unicycler.fastq | \\
        grep -v "^--" | \\
        gzip --no-name > \\
        ${prefix}.unicycler.circular.fastq.gz

    gzip -n ${prefix}.unicycler.fastq
    """
}
