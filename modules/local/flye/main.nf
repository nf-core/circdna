process FLYE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye:2.9.6--py313h7fbb527_1' :
        'quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("assembly.fasta"), emit: assembly
    tuple val(meta), path("assembly_graph.gfa"), emit: gfa
    tuple val(meta), path("assembly_info.txt"), emit: info
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    flye \\
        $args \\
        --out-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$(flye --version 2>&1 | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    cat <<-EOF > assembly.fasta
    >contig_1
    ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
    GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
    EOF

    cat <<-EOF > assembly_graph.gfa
    H	VN:Z:1.0
    S	1	ATCGATCGATCGATCG	RC:i:0
    EOF

    cat <<-EOF > assembly_info.txt
    #seq_name	length	cov	circ	repeat	mult
    contig_1	1000	50	0	0	1
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: "2.9.6"
    END_VERSIONS
    """
}
