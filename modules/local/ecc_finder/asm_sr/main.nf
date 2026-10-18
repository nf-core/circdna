process ECC_FINDER_ASM_SR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/bioinfortools/ecc_finder:1.0.0' :
        'quay.io/bioinfortools/ecc_finder:1.0.0' }"

    input:
    tuple val(meta), path(query1)
    tuple val(meta2), path(query2)

    output:
    tuple val(meta), path("${prefix}.fasta"), emit: fasta
    tuple val("${task.process}"), val('ecc_finder'), eval("ecc_finder.py --version"), topic: versions, emit: versions_ecc_finder

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "ecc.asm.sr"

    """
    # ecc_finder.py 子命令用相对路径调用子脚本，需在 work 目录软链容器脚本
    ln -sf /app/*.py . 2>/dev/null || true
    ecc_finder.py asm-sr \\
        ${query1} \\
        ${query2} \\
        -t ${task.cpus} \\
        -o . \\
        -x ${prefix} \\
        $args

    mv eccFinder_asm_output/${prefix}.fasta ${prefix}.fasta
    """

    stub:
    prefix = task.ext.prefix ?: "ecc.asm.sr"
    """
    touch ${prefix}.fasta
    """
}
