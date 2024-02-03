process CNVKIT_BATCH {
    tag "$meta.id"
    label 'process_low'

    conda 'bioconda::cnvkit=0.9.9'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.9--pyhdfd78af_0' :
        'quay.io/biocontainers/cnvkit:0.9.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  cnn

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.cnn"), emit: cnn, optional: true
    tuple val(meta), path("*.cnr"), emit: cnr, optional: true
    tuple val(meta), path("*.cns"), emit: cns, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def reference_args = cnn ? "--reference $cnn" : ""
    def fasta_args = cnn ? "" : "--fasta $fasta"
""

    """
    cnvkit.py \\
        batch \\
        $bam \\
        $fasta_args \\
        $reference_args \\
        --processes $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e "s/cnvkit v//g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_args = cnn ? "" : "--fasta $fasta"
    def reference_args = cnn ? "--reference $cnn" : ""

    """
    touch ${prefix}.bed
    touch ${prefix}.cnn
    touch ${prefix}.cnr
    touch ${prefix}.cns

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e "s/cnvkit v//g")
    END_VERSIONS
    """
}
