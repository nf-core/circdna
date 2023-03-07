process AMPLICONCLASSIFIER_AMPLICONCLASSIFIER {
    tag "AA Amplicons"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ampliconclassifier=0.4.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampliconclassifier:0.4.14--hdfd78af_0':
        'quay.io/biocontainers/ampliconclassifier:0.4.14--hdfd78af_0' }"

    input:
    path (input_file)

    output:
    path ("*amplicon_classification_profiles.tsv"   ), emit: class_tsv       , optional: true
    path ("*edge_classification_profiles.tsv"       ), emit: edge_tsv        , optional: true
    path ("*gene_list.tsv"                  )        , emit: gene_list       , optional: true
    path ("*ecDNA_counts.tsv"               )        , emit: ecDNA_counts    , optional: true
    path ("*.bed"                           )        , emit: bed             , optional: true
    path ("*annotated_cycles.txt"           )        , emit: annotated_cycles, optional: true
    path ("*class_radar.{png,pdf}"          )        , emit: radar_plot      , optional: true
    path ("*feature_entropy.tsv"            )        , emit: entropy         , optional: true
    path ("*feature_basic_properties.tsv"   )        , emit: basic_properties, optional: true
    path ("*classification_bed_files/*"     )        , emit: bed_files       , optional: true
    path ("*annotated_cycles_files/"        )        , emit: cycles_files    , optional: true
    path ("*.classifier_stdout.log"         )        , emit: log             , optional: true
    path ("*"                               )        , emit: all             , optional: true
    path ("versions.yml"                    )        , emit: versions        , optional: true

    script:
    def args = task.ext.args ?: ''

    """
    REF=${params.reference_build}
    export AA_DATA_REPO=${params.aa_data_repo}
    export AA_SRC=${projectDir}/bin

    amplicon_classifier.py \\
        --ref \$REF \\
        $args \\
        --input $input_file \\
        > ampliconclassifier.classifier_stdout.log

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

    touch "${prefix}.amplicon_classification_profiles.tsv"
    touch "${prefix}.classifier_stdout.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliconClassifier: \$(echo \$(amplicon_classifier.py --version | sed 's/amplicon_classifier //g' | sed 's/ .*//g'))
    END_VERSIONS
    """
}
