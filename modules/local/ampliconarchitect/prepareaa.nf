// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process AMPLICONARCHITECT_PREPAREAA {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::python=3.8 bioconda::cnvkit=0.9.9 bioconda::canvas=1.35.1.1316 bioconda::freebayes=1.3.5 numpy matplotlib intervaltree" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "quay.io/biocontainers/YOUR-TOOL-HERE"
    }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bed"), emit: bed

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    cnvkit=\$(which cnvkit.py)
    python3=\$(which python3)
    AA_DATA_REPO=${params.aa_data_repo}
    MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}
    REF=${params.reference_build}

    PrepareAA.py \\
        -s ${meta.id} \\
        -t ${task.cpus} \\
        $options.args \\
        --sorted_bam $bam \\
        --ref \$REF \\
        --cnvkit_dir \$cnvkit \\
        --python3_path \$python3
    """
}
