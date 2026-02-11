process MULTIQC {
    label 'process_medium'

    conda 'bioconda::multiqc=1.19'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.18--pyhdfd78af_0' :
        'biocontainers/multiqc:1.19--pyhdfd78af_0' }"

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

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_multiqc_config = task.ext.extra_multiqc_config ?: null
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def multiqc_logo = task.ext.multiqc_logo ?: null
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace_names = task.ext.replace_names ?: null
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def sample_names = task.ext.sample_names ?: null
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    """
    multiqc \\
        -f \\
        $args \\
        $config \\
        $prefix \\
        $extra_config \\
        $logo \\
        $replace \\
        $samples \\
        .
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    mkdir multiqc_plots
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
