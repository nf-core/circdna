process BEDTOOLS_SPLITBAM2BED {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1' :
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2'}"

    input:
    tuple val(meta), path(split_bam)

    output:
    tuple val(meta), path("*.split.txt"), emit: split_txt
    path  "versions.yml"          , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # 输出 9 个空白分隔字段：chr start end readID readNo mapq strand cigar tag
    # CIRCLEFINDER 的后续 awk 以 FS=" " 解析，字段位置是硬约束：
    #   \$1=chr \$2=start \$3=end \$4=readID \$6=mapq \$7=strand
    # 且两行合并后 \$10 必须是第二行的 chr —— 所以每行必须恰好 9 个字段。
    # 旧实现的三个问题：
    #   1) CIGAR 在 bedtools bamtobed -cigar 的第 7 列，却去展开/匹配第 8 列 → 恒不匹配；
    #   2) 依赖 read name 的 _1/1、_2/2 后缀，而 SRA 数据名为 ERR10889838.1234 → 永不命中；
    #   3) 字段数因此不足 9，step8 的 \$10/\$16 全部错位。
    # 结果 split.txt 恒为空 → CIRCLEFINDER 以 "No split reads found" 提前退出（全样本 0 检出）。
    bedtools bamtobed $args -i $split_bam | \\
        awk 'BEGIN{FS=OFS="\t"} {
            name=\$4; rn="0";
            if (substr(name, length(name)-3) == "_1/1") { rn="1"; name=substr(name,1,length(name)-4) }
            else if (substr(name, length(name)-3) == "_2/2") { rn="2"; name=substr(name,1,length(name)-4) }
            else if (substr(name, length(name)-1) == "/1") { rn="1"; name=substr(name,1,length(name)-2) }
            else if (substr(name, length(name)-1) == "/2") { rn="2"; name=substr(name,1,length(name)-2) }
            cig=\$7; last=substr(cig, length(cig), 1); tag="";
            if (cig ~ /^[0-9]+M/ && (last == "S" || last == "H")) tag="first";
            else if (cig ~ /^[0-9]+[SH]/ && last == "M") tag="second";
            if (tag == "") next;
            print \$1, \$2, \$3, name, rn, \$5, \$6, cig, tag
        }' > '${prefix}.txt'

    # Software Version
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch '${prefix}.txt'
    touch versions.yml
    """
}
