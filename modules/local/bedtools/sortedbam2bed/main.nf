process BEDTOOLS_SORTEDBAM2BED {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1':
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
    tuple val(meta), path(sorted_bam), path(sorted_bai)

    output:
    tuple val(meta), path("*.concordant.txt"), emit: conc_txt
    path  "versions.yml"          , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # 输出 8 个空白分隔字段：chr start end readID readNo mapq strand cigar
    # CIRCLEFINDER step10 会按 tab 第 8 列取 CIGAR 做 M/S/H 展开并分类
    # （first/second/confusing），所以 CIGAR 必须落在第 8 列。
    # 旧实现把 bedtools bamtobed -cigar 的第 7 列（CIGAR）写进第 7 列、
    # 第 8 列留空，导致 step10 的分类 awk 恒不匹配 → 终集恒为空。
    # 另外把 readNo 显式成列，避免依赖 read name 的 _1/1、_2/2 后缀。
    bedtools bamtobed $args -i $sorted_bam | \\
        awk 'BEGIN{FS=OFS="\t"} {
            name=\$4; rn="0";
            if (substr(name, length(name)-3) == "_1/1") { rn="1"; name=substr(name,1,length(name)-4) }
            else if (substr(name, length(name)-3) == "_2/2") { rn="2"; name=substr(name,1,length(name)-4) }
            else if (substr(name, length(name)-1) == "/1") { rn="1"; name=substr(name,1,length(name)-2) }
            else if (substr(name, length(name)-1) == "/2") { rn="2"; name=substr(name,1,length(name)-2) }
            print \$1, \$2, \$3, name, rn, \$5, \$6, \$7
        }' > '${prefix}.concordant.txt'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch '${prefix}.concordant.txt'
    touch versions.yml
    """
}
