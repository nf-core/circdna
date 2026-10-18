process FLED {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/bioinfortools/fled:1.7.0'

    input:
    tuple val(meta), path(fastq)
    path genome_fasta

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.fled_junctions.txt"), emit: junctions
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.DiGraph.*Junction.fa"), emit: sequences
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    # FLED requires an indexed reference genome (.fai) and plain fasta.
    # Gzipped inputs are decompressed since FLED detects format by extension.
    if [[ "${genome_fasta}" == *.gz ]]; then
        zcat "${genome_fasta}" > reference.fa
        GENOME="reference.fa"
    else
        GENOME="${genome_fasta}"
    fi
    if [ ! -f \${GENOME}.fai ]; then
        samtools faidx \${GENOME}
    fi

    if [[ "${fastq}" == *.gz ]]; then
        zcat "${fastq}" > reads_input.fastq
        FQ="reads_input.fastq"
    else
        FQ="${fastq}"
    fi
    # PseudoReference dedupe / CSI / --split-prefix fixes are baked into fled:1.7.0.
    FLED Detection \\
        -ref \${GENOME} \\
        -fq \${FQ} \\
        -o ${prefix} \\
        -dir . \\
        -t ${task.cpus} \\
        $args

    # Combine both junction outputs into a single candidates file.
    # MulsegFullJunction.out is only produced when multi-segment eccDNAs exist.
    # `set -e` 下 cat 缺失文件会中止脚本，需 || true 容错。
    cat ${prefix}.DiGraph.OnesegJunction.out 2>/dev/null > ${prefix}.fled_junctions.txt || true
    cat ${prefix}.DiGraph.MulsegFullJunction.out 2>/dev/null >> ${prefix}.fled_junctions.txt || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fled: \$(FLED 2>&1 | grep -oP 'version=\\K[\\d.]+' || echo 'unknown')
        samtools: \$(samtools --version 2>&1 | grep -oP 'samtools \\K[\\d.]+')
    END_VERSIONS
    """

    stub:
    def stub_prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${stub_prefix}.fled_junctions.txt
    touch ${stub_prefix}.DiGraph.OnesegJunction.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fled: "1.7.0"
        samtools: "1.20"
    END_VERSIONS
    """
}
