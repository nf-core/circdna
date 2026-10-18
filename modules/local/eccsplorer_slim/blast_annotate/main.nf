process ECCSPLORER_BLAST_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.17.0--d4fb881691596759' :
        'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759' }"

    input:
    tuple val(meta), path(cand_fa)
    tuple val(meta2), path(db_dir)
    val(db_name)

    output:
    tuple val(meta), path("${prefix}_blast.m6"), emit: blast_m6
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # per-candidate blastn annotation (replicates eccMapper analyze_candidate_region)
    blastn \\
        -query ${cand_fa} \\
        -db ${db_name} \\
        -out ${prefix}_blast.m6 \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \\
        -max_hsps 25 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | grep -oP 'blastn: \\K[0-9.]+' | head -1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_blast.m6
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: 2.17.0
    END_VERSIONS
    """
}
