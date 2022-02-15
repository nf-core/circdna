// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/software
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join

// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

params.options = [:]
options        = initOptions(params.options)

process BEDTOOLS_SORTEDBAM2BED {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2"
    }

    input:
    tuple val(meta),
            path(sorted_bam),
            path(sorted_bai)

    output:
    tuple val(meta),
        path("*.concordant.txt"),
        emit: conc_txt

    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bedtools bamtobed $options.args -i $sorted_bam | \
        sed -e 's/\\// /g' | \
        awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\\n",\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8)}' > '${prefix}.concordant.txt'

    # Software Version
    bedtools --version > ${software}.version.txt

    """
}
