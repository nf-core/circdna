process ECCSPLORER_COMPARATIVE_BLAST {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.17.0--d4fb881691596759' :
        'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759' }"

    input:
    tuple val(meta), path(map_cand_fa)     // mapper candidate sequences (ECC-SEQUENCES.fasta)
    tuple val(meta2), path(cluster_contigs) // cluster contig fastas (per-cluster)

    output:
    tuple val(meta), path("${prefix}_comparative-blast.m6"), emit: blast_m6
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # makeblastdb from mapper candidates (replicates eccComparer make_mapper_fasta_db)
    makeblastdb -in ${map_cand_fa} -dbtype nucl -out ${prefix}_mapperdb -title mapperdb $args

    # blast each cluster contig vs mapper candidates (replicates eccComparer blast_cand)
    : > ${prefix}_comparative-blast.m6
    for f in ${cluster_contigs}; do
        cl=\$(basename \$f | sed 's/\\..*//')
        blastn -query \$f -db ${prefix}_mapperdb \\
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \\
            -max_hsps 25 >> ${prefix}_comparative-blast.m6 2>/dev/null || true
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | grep -oP 'blastn: \\K[0-9.]+' | head -1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_comparative-blast.m6
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: 2.17.0
    END_VERSIONS
    """
}
