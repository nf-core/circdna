process CNVKIT_SEGMENT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::cnvkit=0.9.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cnvkit:0.9.9--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/cnvkit:0.9.9--pyhdfd78af_0"
    }

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/software/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(cnr)

    output:
    tuple val(meta), path("*.cns"), emit: cns
    // path "*.version.txt"          , emit: version

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
