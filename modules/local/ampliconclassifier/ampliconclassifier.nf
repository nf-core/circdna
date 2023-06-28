process AMPLICONCLASSIFIER_AMPLICONCLASSIFIER {
    tag "AA Amplicons"
    label 'process_low'

    conda "conda-forge::python=3.7 bioconda::pysam=0.16.0 anaconda::flask=2.2.2 conda-forge::numpy=1.21.6 conda-forge::matplotlib=3.2.2 anaconda::scipy=1.7.3 conda-forge::intervaltree=3.0.2 anaconda::future=0.18.2 mosek::mosek=9.0.88"
    container '/home/local/BICR/dschreye/ampliconsuite.sif'

    input:
    path (graphs)
    path (cycles)
    path (cnseg)

    output:
    path ("*amplicon_classification_profiles.tsv"   ), emit: class_tsv       , optional: true
    path ("*edge_classification_profiles.tsv"       ), emit: edge_tsv        , optional: true
    path ("*gene_list.tsv"                  )        , emit: gene_list       , optional: true
    path ("*ecDNA_counts.tsv"               )        , emit: ecDNA_counts    , optional: true
    path ("*.bed"                           )        , emit: bed             , optional: true
    path ("*annotated_cycles.txt"           )        , emit: annotated_cycles, optional: true
    path ("*class_radar.{png,pdf}"          )        , emit: radar_plot      , optional: true
    path ("*feature_entropy.tsv"            )        , emit: entropy         , optional: true
    path ("*features_to_graph.txt"          )        , emit: features_to_graph, optional: true
    path ("*feature_basic_properties.tsv"   )        , emit: basic_properties, optional: true
    path ("*classification_bed_files/*"     )        , emit: bed_files       , optional: true
    path ("*annotated_cycles_files/"        )        , emit: cycles_files    , optional: true
    path ("*.classifier_stdout.log"         )        , emit: log             , optional: true
    path ("*"                               )        , emit: all             , optional: true
    path ("versions.yml"                    )        , emit: versions        , optional: true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "ampliconarchitect"

    """
    REF=${params.reference_build}
    export AA_DATA_REPO=${params.aa_data_repo}
    export AA_SRC=\$(dirname \$(readlink -f \$(which AmpliconArchitect.py)))
    export AC_SRC=\$(dirname \$(readlink -f \$(which amplicon_classifier.py)))

    AmpliconSuite-pipeline.py \\
        -s $prefix \\
        --completed_AA_runs ./ \\
        -t $task.cpus \\
        --ref "GRCh38"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    export AA_DATA_REPO=${params.aa_data_repo}
    export MOSEKLM_LICENSE_FILE=${params.mosek_license_dir}
    export AA_SRC=${projectDir}/bin
    REF=${params.reference_build}

    touch "ampliconclassifier_amplicon_classification_profiles.tsv"
    touch "ampliconclassifier_classifier_stdout.log"

    amplicon_classifier.py --help

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
    END_VERSIONS
    """
}
