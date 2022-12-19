process CNVKIT_BATCH {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::cnvkit=0.9.9' : null)
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
    def fasta_args = reference ? "" : "--fasta $fasta"
    def reference_args = reference ? "--reference $cnn" : ""

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
    def reference = params.aa_data_repo + "/" + params.reference_build + "/" + params.reference_build + "_cnvkit_filtered_ref.cnn"
    def fasta_args = reference ? "" : "--fasta $fasta"
    def reference_args = reference ? "--reference $cnn" : ""

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
