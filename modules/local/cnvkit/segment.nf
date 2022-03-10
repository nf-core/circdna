process CNVKIT_SEGMENT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::cnvkit=0.9.9' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.9--pyhdfd78af_0' :
        'quay.io/biocontainers/cnvkit:0.9.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(cnr)

    output:
    tuple val(meta), path("*.cns"), emit: cns
    tuple val(meta), path("*.cnr"), emit: cnr
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cnvkit.py \\
        segment \\
        $cnr \\
        -p $task.cpus \\
        -m "cbs" \\
        -o ${prefix}.cnvkit.segment.cns
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e "s/cnvkit v//g")
    END_VERSIONS
    """
}
