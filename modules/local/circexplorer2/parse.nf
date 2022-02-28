// Import generic module functions
process CIRCEXPLORER2_PARSE {
    tag "$meta.id"
    label 'process_medium'
    conda (params.enable_conda ? "bioconda::circexplorer2=2.3.8" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/circexplorer2:2.3.8--pyh864c0ab_1"
    } else {
        container "quay.io/biocontainers/circexplorer2:2.3.8--pyh864c0ab_1"
    }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bed")  , emit: bed
    path "*.log"                    , emit: log
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    CIRCexplorer2 parse $args $bam -b ${prefix}.temp.bed > ${prefix}_CIRCexplorer2_parse.log
    cat ${prefix}.temp.bed | tr "/" "\t" > ${prefix}.circexplorer_circdna.bed
    rm ${prefix}.temp.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CIRCexplorer2: \$(echo \$(CIRCexplorer2 --version))
    END_VERSIONS
    """
}
