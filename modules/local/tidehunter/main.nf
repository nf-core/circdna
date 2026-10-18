process TIDEHUNTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tidehunter:1.5.6--h7f5d12c_0' :
        'quay.io/biocontainers/tidehunter:1.5.6--h7f5d12c_0' }"

    input:
    tuple val(meta), path(reads)
    path adapter5     // optional: 5' adapter FASTA (with adapter3 enables full-length mode)
    path adapter3     // optional: 3' adapter FASTA

    output:
    tuple val(meta), path("${prefix}.fasta"),      emit: cons_fa,   optional: true
    tuple val(meta), path("${prefix}.cons.out"),   emit: cons_tab,  optional: true
    tuple val(meta), path("${prefix}.cons.fq"),    emit: cons_fq,   optional: true
    tuple val(meta), path("${prefix}.unit.fasta"), emit: unit_fa,   optional: true
    tuple val(meta), path("${prefix}.unit.out"),   emit: unit_tab,  optional: true
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_unit = args.contains('-u') || args.contains('--unit-seq')
    def is_tab  = args.contains('-f 2')
    def is_fq   = args.contains('-f 3')
    def out_name = is_unit ? (is_tab ? "${prefix}.unit.out" : "${prefix}.unit.fasta") :
                    (is_tab ? "${prefix}.cons.out" : is_fq ? "${prefix}.cons.fq" : "${prefix}.fasta")
    def adapter_args = (adapter5 && adapter3) ? "-5 ${adapter5} -3 ${adapter3}" : ''
    """
    # TideHunter: de novo tandem-repeat detection from long reads (no reference)
    TideHunter ${adapter_args} ${reads} $args > ${out_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidehunter: 1.5.6
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    touch ${prefix}.cons.out
    touch ${prefix}.cons.fq
    touch ${prefix}.unit.fasta
    touch ${prefix}.unit.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidehunter: 1.5.6
    END_VERSIONS
    """
}
