process AMPLICONARCHITECT_AMPLICONARCHITECT {
    tag "$meta.id"
    label 'process_low'

   // conda "conda-forge::python=2.7 bioconda::pysam=0.15.2 anaconda::flask=1.1.2 anaconda::cython=0.29.14 anaconda::numpy=1.16.6 anaconda::scipy=1.2.1 conda-forge::matplotlib=2.2.5 mosek::mosek=8.0.60 anaconda::future=0.18.2"
   // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
   //     'https://depot.galaxyproject.org/singularity/mulled-v2-6eefa51f13933d65b4f8155ca2f8cd81dea474ba:baa777f7c4e89a2ec4d1eab7d424a1f46503bac7-0':
   //     'quay.io/biocontainers/mulled-v2-6eefa51f13933d65b4f8155ca2f8cd81dea474ba:baa777f7c4e89a2ec4d1eab7d424a1f46503bac7-0' }"
    conda "bioconda::ampliconsuite=0.1555.2 mosek::mosek=10.1b1"
    container '/home/local/BICR/dschreye/ampliconsuite.sif'

    input:
    tuple val(meta), path(bam), path(bai), path(bed)

    output:
    path "versions.yml"                     , emit: versions
    tuple val(meta), path("*cycles.txt")    , optional: true, emit: cycles
    tuple val(meta), path("*graph.txt")     , optional: true, emit: graph
    tuple val(meta), path("*cnseg.txt")    , optional: true, emit: cnseg
    tuple val(meta), path("*.out")          , optional: true, emit: out
    tuple val(meta), path("*.{pdf,png}")    , optional: true, emit: svview
    tuple val(meta), path("*_summary.txt")  , optional: true, emit: summary
    tuple val(meta), path("*")              , optional: true, emit: all

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export AA_DATA_REPO=${params.aa_data_repo}
    export MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}
    # export AA_SRC=${projectDir}/bin
    export AA_SRC=\$(dirname \$(readlink -f \$(which AmpliconArchitect.py)))
    export AC_SRC=\$(dirname \$(readlink -f \$(which amplicon_classifier.py)))
    REF=${params.reference_build}

    AmpliconSuite-pipeline.py \\
    -t $task.cpus \\
    --bam $bam \\
    --bed $bed \\
    --ref \$REF \\
    -s "${prefix}" \\
    --run_AA
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconArchitect: \$(echo \$(AmpliconArchitect.py --version 2>&1) | sed 's/AmpliconArchitect version //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export AA_DATA_REPO=${params.aa_data_repo}
    export MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}
    export AA_SRC=${projectDir}/bin
    REF=${params.reference_build}

    touch "${prefix}.logs.txt"
    touch "${prefix}.cycles.txt"
    touch "${prefix}.graph.txt"
    touch "${prefix}.out"
    touch "${prefix}_cnseg.txt"
    touch "${prefix}.pdf"
    touch "${prefix}.png"
    touch "${prefix}_summary.txt"

    AmpliconArchitect.py --help

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconArchitect: \$(echo \$(AmpliconArchitect.py --version 2>&1) | sed 's/AmpliconArchitect version //g')
    END_VERSIONS
    """
}
