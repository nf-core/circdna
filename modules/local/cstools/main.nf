process CSTOOLS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/bioinfortools/ciderseq2:2.0.0'

    input:
    tuple val(meta), path(reads)
    path(config)

    output:
    tuple val(meta), path("${prefix}*"), emit: output
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args carries the sub-command + options, e.g.:
    //   "split --format fastq --numjobs 4 --mode full"
    //   "join --clean"
    //   "plot"
    def args = task.ext.args ?: 'plot'
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # cs-tools.py sub-commands ask for interactive confirmation (click.prompt /
    # click.confirm).  Answer with empty lines so click falls back to the
    # option defaults and the confirmation is declined (split only creates the
    # job files; the actual ciderseq jobs are run by Nextflow itself).
    printf '\\n\\n' | cs-tools.py \\
        $args \\
        ${config} \\
        ${reads}

    # collect everything the tool produced (split dir, plots, joined results)
    if [ -d "${reads}.dir" ]; then
        mv "${reads}.dir" ${prefix}.dir
    fi
    if [ -d plots ]; then
        mkdir -p ${prefix}_plots
        mv plots/* ${prefix}_plots/ 2>/dev/null || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciderseq2: 2.0.0
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_cstools.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciderseq2: 2.0.0
    END_VERSIONS
    """
}
