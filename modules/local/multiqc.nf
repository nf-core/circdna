process MULTIQC {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::multiqc=1.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0' }"

    input:
    path multiqc_config
    path multiqc_custom_config
    path software_versions
    path workflow_summary
    path ('fastqc/*')
    path ('trimgalore/fastqc/*')
    path ('trimgalore/*')
    path ('samtools/stats/*')
    path ('samtools/flagstat/*')
    path ('samtools/idxstats/*')
    path ('picard/markduplicates/stats/*')
    path ('picard/markduplicates/flagstat/*')
    path ('picard/markduplicates/idxstats/*')
    path ('picard/markduplicates/metrics/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    script:
    def args = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc \\
        -f \\
        $args \\
        $custom_config \\
        .
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
