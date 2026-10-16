process CIRCLESEEKER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circleseeker:1.6.1--pyhdfd78af_0' :
        'quay.io/biocontainers/circleseeker:1.6.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}_eccDNA_summary.csv"), emit: merged
    tuple val(meta), path("${prefix}.circleseeker.bed"),   emit: bed
    tuple val(meta), path("${prefix}_summary.txt"),        emit: summary
    tuple val(meta), path("${prefix}_report.html"),        emit: report
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # CircleSeeker requires FASTA input (pysam.FastaFile) and detects format
    # by file extension, so decompress gz and convert FASTQ to FASTA first.
    if [[ "${reads}" == *.gz ]]; then
        zcat "${reads}" > reads_input.fastq
        READS_IN="reads_input.fastq"
    else
        READS_IN="${reads}"
    fi
    if [[ "\${READS_IN}" == *.fastq || "\${READS_IN}" == *.fq ]]; then
        awk 'NR%4==1 {print ">" substr(\$0,2)} NR%4==2 {print}' "\${READS_IN}" > reads_input.fasta
        READS_IN="reads_input.fasta"
    fi

    # Reference genome: decompress if gz (CircleSeeker auto-creates .mmi/.fai
    # indexes next to the file, so run it in the writable workdir).
    if [[ "${fasta}" == *.gz ]]; then
        zcat "${fasta}" > reference.fa
        FASTAIN="reference.fa"
    else
        FASTAIN="${fasta}"
    fi

    circleseeker \\
        -i \${READS_IN} \\
        -r \${FASTAIN} \\
        -o . \\
        -p ${prefix} \\
        -t ${task.cpus} \\
        $args

    # 无候选时 CircleSeeker 不生成 summary/report 文件；先补空文件再转 BED（§8）
    if [ ! -f ${prefix}_eccDNA_summary.csv ]; then
        echo "eccDNA_id,type,state,chr,start,end,strand,length" > ${prefix}_eccDNA_summary.csv
    fi
    touch ${prefix}_summary.txt
    touch ${prefix}_report.html

    # Convert the eccDNA summary CSV into a BED6+read_count table so the
    # standard LONG_READ_FILTERING machinery can be reused downstream.
    circleseeker_to_bed.py \\
        ${prefix}_eccDNA_summary.csv \\
        ${prefix}.circleseeker.bed

    # External aligners spawned by CircleSeeker can outlive the main process
    # and keep the container PID namespace alive; wait for any remaining children.
    wait || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circleseeker: \$(circleseeker --version 2>&1 | sed 's/CircleSeeker //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_eccDNA_summary.csv
    touch ${prefix}.circleseeker.bed
    touch ${prefix}_summary.txt
    touch ${prefix}_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circleseeker: "stub"
    END_VERSIONS
    """
}
