process ECCSPLORER_BLAST_COMBINEDDB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.17.0--d4fb881691596759' :
        'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759' }"

    input:
    tuple val(meta), path(db_fa)

    output:
    tuple val(meta), path("${prefix}_db*"), emit: combined_db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # makeblastdb (replicates ECCsplorer basic_setup; single-db equivalent of
    # makeblastdb + blastdb_aliastool when -d is one database)
    makeblastdb -in ${db_fa} -dbtype nucl -out ${prefix}_db -title ${prefix}_db $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | grep -oP 'blastn: \\K[0-9.]+' | head -1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_db.nhr
    touch ${prefix}_db.nsq
    touch ${prefix}_db.nin
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: 2.17.0
    END_VERSIONS
    """
}
