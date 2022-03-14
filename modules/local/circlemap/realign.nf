process CIRCLEMAP_REALIGN {
    tag "$meta.id"
    label 'process_high'
    conda (params.enable_conda ? "bioconda::circle-map=1.1.4 biopython=1.77" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-3423532b7f92d89bbf04791d2e48d40cb0a4111e:2a8abad654fa8c270bfcc01eef47487e2573badb-0':
        'quay.io/biocontainers/mulled-v2-3423532b7f92d89bbf04791d2e48d40cb0a4111e:2a8abad654fa8c270bfcc01eef47487e2573badb-0' }"

    input:
    tuple val(meta), path(re_bam), path(re_bai), path(qname), path(sbam), path(sbai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bed"), emit: bed, optional: true
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    circle_map.py \\
        Realign \\
        $args \\
        -i $re_bam \\
        -qbam $qname \\
        -sbam $sbam \\
        -fasta $fasta \\
        --threads $task.cpus \\
        -o ${prefix}_circularDNA_coordinates.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Circle-Map: \$(echo \$(Circle-Map --help 2<&1 | grep -o "version=[0-9].[0-9].[0-9]" | sed 's/version=//g'))
    END_VERSIONS
    """
}
