process ECC_FINDER_MAP_SR {
    tag "$meta2.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/bioinfortools/ecc_finder:1.0.0' :
        'quay.io/bioinfortools/ecc_finder:1.0.0' }"

    input:
    tuple val(meta), path(idx)
    tuple val(meta2), path(query1)
    tuple val(meta3), path(query2)
    tuple val(meta4), path(ref)

    output:
    tuple val(meta), path("${prefix}.csv"), emit: csv
    tuple val(meta), path("${prefix}.fasta"), emit: fasta
    tuple val("${task.process}"), val('ecc_finder'), eval("ecc_finder.py --version"), topic: versions, emit: versions_ecc_finder

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "ecc.sr"

    """
    # ecc_finder.py 子命令用相对路径调用子脚本，需在 work 目录软链容器脚本
    ln -sf /app/*.py . 2>/dev/null || true
    # bwa 索引目录 → bwa 索引前缀文件（map-sr.py 校验 idx 为文件；bwa 用相邻 .bwt/.pac 索引）
    IDX_BWT=\$(ls ${idx}/*.bwt 2>/dev/null | head -1)
    IDX_PREFIX=\${IDX_BWT%.bwt}
    [ -n "\$IDX_PREFIX" ] && touch "\$IDX_PREFIX"
    ecc_finder.py map-sr \\
        \$IDX_PREFIX \\
        ${query1} \\
        ${query2} \\
        -r ${ref} \\
        -t ${task.cpus} \\
        -o . \\
        -x ${prefix} \\
        $args

    mv eccFinder_output/${prefix}.csv ${prefix}.csv 2>/dev/null || true
    mv eccFinder_output/${prefix}.fasta ${prefix}.fasta 2>/dev/null || true
    """

    stub:
    prefix = task.ext.prefix ?: "ecc.sr"
    """
    touch ${prefix}.csv
    touch ${prefix}.fasta
    """
}
