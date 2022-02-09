// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)
process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::minimap2=2.21' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:e1ea28074233d7265a5dc2111d6e55130dff5653-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:e1ea28074233d7265a5dc2111d6e55130dff5653-0' }"

    input:
    path reference
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.paf"), emit: paf
    path "versions.yml" , emit: versions
    path "*"

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        $reference \\
        $reads \\
        > ${prefix}.unicycler.minimap2.sam


    samtools sort -@ $task.cpus -o ${prefix}.unicycler.minimap2.sorted.sam
    paftools.js sam2paf ${prefix}.unicycler.minimap2.sorted.sam > ${prefix}.unicycler.minimap2.paf
    samtools view -S -b ${prefix}.unicycler.minimap2.sorted.sam > ${prefix}.unicycler.minimap2.sorted.bam
    # rm ${prefix}.unicycler.minimap2.sam ${prefix}.unicycler.minimap2.sorted.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
