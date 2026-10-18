process ECC_FINDER_ASM_ONT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/bioinfortools/ecc_finder:1.0.0' :
        'quay.io/bioinfortools/ecc_finder:1.0.0' }"

    input:
    tuple val(meta), path(query)

    output:
    tuple val(meta), path("${prefix}.fasta"), emit: fasta
    tuple val("${task.process}"), val('ecc_finder'), eval("ecc_finder.py --version"), topic: versions, emit: versions_ecc_finder

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "ecc.asm.ont"

    """
    # ecc_finder.py 子命令用相对路径调用子脚本，需在 work 目录软链容器脚本
    ln -sf /app/*.py . 2>/dev/null || true
    ecc_finder.py asm-ont \\
        ${query} \\
        -t ${task.cpus} \\
        -o . \\
        -x ${prefix} \\
        $args

    # -o . 时 ecc_finder 直接写 work 根目录；部分版本写入 eccFinder_asm_output/ 子目录
    if [ -d eccFinder_asm_output ]; then
        mv eccFinder_asm_output/${prefix}.fasta ${prefix}.fasta
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "ecc.asm.ont"
    """
    touch ${prefix}.fasta
    """
}
