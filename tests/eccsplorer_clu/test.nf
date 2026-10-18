#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ECCSPLORER_CLU_PREPARE   } from '../../modules/local/eccsplorer_slim/clu_prepare/main'
include { ECCSPLORER_CLU_CANDIDATES } from '../../modules/local/eccsplorer_slim/clu_candidates/main'
include { REPEATEXPLORER2          } from '../../modules/local/repeatexplorer2/main'

workflow {
    def data_dir = "${projectDir}/../../testdatasets/ngs"

    ch_treatment = channel.of(
        [ [id: 'circdna_1', pair: 'p1'],
            file("${data_dir}/circdna_1_R1.fastq.gz"),
            file("${data_dir}/circdna_1_R2.fastq.gz") ]
    )
    ch_control = channel.of(
        [ [id: 'gdna_1', pair: 'p1'],
            file("${data_dir}/gdna_1_R1.fastq.gz"),
            file("${data_dir}/gdna_1_R2.fastq.gz") ]
    )
    taxon = 'vir'

    ECCSPLORER_CLU_PREPARE ( ch_treatment, ch_control, taxon )
    REPEATEXPLORER2 (
        ECCSPLORER_CLU_PREPARE.out.ready_fa,
        taxon
    )
    ECCSPLORER_CLU_CANDIDATES (
        REPEATEXPLORER2.out.clu_dir,
        'TR',
        'CO'
    )

    ECCSPLORER_CLU_CANDIDATES.out.cluster_candidates
        .view { _meta, csv -> "CLU CANDIDATES OK: ${csv}" }
    ECCSPLORER_CLU_CANDIDATES.out.comparative_table
        .view { _meta, csv -> "CLU TABLE OK: ${csv}" }
}
