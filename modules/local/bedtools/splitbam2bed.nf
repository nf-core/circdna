process BEDTOOLS_SPLITBAM2BED {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2"
    }

    input:
    tuple val(meta), path(split_bam)

    output:
    tuple val(meta), path("*.split.txt"), emit: split_txt
    path  "versions.yml"          , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools bamtobed $args -i $split_bam | \
        sed -e 's/_2\\/2/ 2/g' | \
        sed -e 's/_1\\/1/ 1/g' |
        awk '{printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\\n", \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8}' |
        awk 'BEGIN{FS=OFS="\t"} {gsub("M", " M ", \$8)} 1' | \
        awk 'BEGIN{FS=OFS="\t"} {gsub("S", " S ", \$8)} 1' | \
        awk 'BEGIN{FS=OFS="\t"} {gsub("H", " H ", \$8)} 1' | \
        awk 'BEGIN{FS=OFS=" "} {if ((\$9=="M" && \$NF=="H") || \
        (\$9=="M" && \$NF=="S"))  {printf ("%s\tfirst\\n",\$0)} else if ((\$9=="S" && \$NF=="M") || \
        (\$9=="H" && \$NF=="M")) {printf ("%s\tsecond\\n",\$0)} }' | \
        awk 'BEGIN{FS=OFS="\t"} {gsub(" ", "", \$8)} 1' > '${prefix}.split.txt'

    # Software Version
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version)
    END_VERSIONS

    """
}
