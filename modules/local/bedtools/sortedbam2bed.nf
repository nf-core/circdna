process BEDTOOLS_SORTEDBAM2BED {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.30.0"
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
    bedtools bamtobed $args -i $sorted_bam | \
        sed -e 's/\\// /g' | \
        awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\\n",\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8)}' > '${prefix}.concordant.txt'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
