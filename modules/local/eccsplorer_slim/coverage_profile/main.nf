process ECCSPLORER_COVERAGE_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eccsplorer_slim:1.0.0--0' :
        'quay.io/bioinfortools/eccsplorer_slim:1.0.0' }"

    input:
    tuple val(meta), path(candidates_bed)
    tuple val(meta2), path(beds)
    val(names)
    tuple val(meta3), path(ref_fa)

    output:
    tuple val(meta), path("*_coverage.tsv"), emit: coverage
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # v2: multi-alignment input (TR/CO × all/SR/DR)
    alignments_list=\$(ls ${beds} | sed 's/^/./; s/\$/,/g' | tr -d '\n' | sed 's/,\$//')
    names_list="${names}"
    if [ \$(echo "\$names_list" | tr ',' '\\n' | wc -l) -ne \$(echo "\$alignments_list" | tr ',' '\\n' | wc -l) ]; then
        # fallback: derive names from file basenames
        names_list=\$(ls ${beds} | sed 's/\\.bed\$//' | paste -sd, -)
    fi

    coverage_profile.py \\
        --candidates ${candidates_bed} \\
        --alignments "\${alignments_list}" \\
        --ref ${ref_fa} \\
        --outdir cov_${prefix} \\
        --names "\${names_list}" \\
        $args

    # coverage_profile.py writes <candidate_id>_coverage.tsv into per-candidate subdirs;
    # move them to the process root so the module glob matches
    mv cov_${prefix}/*/*_coverage.tsv . 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eccsplorer_slim: 1.0.0
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}_coverage"
    touch "${prefix}_coverage.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eccsplorer_slim: 1.0.0
    END_VERSIONS
    """
}
