/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCircdna.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input =  Channel.fromPath(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta =  Channel.fromPath(params.fasta) } else { exit 1, 'Fasta reference genome not specified!' }

// Modify fasta channel to include meta data
ch_fasta_meta = ch_fasta.map{ it -> [[id:it[0].baseName], it] }

branch = params.circle_identifier.split(",")
run_circexplorer2 = ("circexplorer2" in branch)
run_circle_map_realign = ("circle_map_realign" in branch)
run_circle_map_repeats = ("circle_map_repeats" in branch)
run_circle_finder = ("circle_finder" in branch)
run_ampliconarchitect = ("ampliconarchitect" in branch)
run_unicycler = ("unicycler" in branch)

if (!(run_unicycler | run_circle_map_realign | run_circle_map_repeats | run_circle_finder | run_ampliconarchitect | run_circexplorer2)) {
    exit 1, 'circle_identifier param not valid. Please check!'
}

if (run_unicycler && !params.input_format == "FASTQ") {
        exit 1, 'Unicycler needs FastQ input. Please specify input_format == "FASTQ", if possible, or don`t run unicycler.'
}

if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check if BWA Index is given
if (params.bwa_index) {
    ch_bwa_index = Channel.fromPath(params.bwa_index).collect()
    bwa_index_exists = true
    } else {
        ch_bwa_index = Channel.empty()
        bwa_index_exists = false
    }

// AMPLICON ARCHITECT INPUT
if (run_ampliconarchitect) {
    mosek_license_dir = file(params.mosek_license_dir)
    if (!mosek_license_dir.exists()) {
        exit 1, "Mosek License Directory is missing! Please specifiy directory containing mosek license using --mosek_license_dir and rename license to 'mosek.lic'."
    }
    if (!params.aa_data_repo) { exit 1, "AmpliconArchitect Data Repository Missing! Please see https://github.com/jluebeck/AmpliconArchitect for more information and specify --aa_data_repo." }
    if (params.reference_build != "hg19" & params.reference_build != "GRCh38" & params.reference_build != "GRCh37" & params.reference_build != "mm10"){
        exit 1, "Reference Build not given! Please specify --reference_build 'mm10', 'hg19', 'GRCh38', or 'GRCh37'."
    }

    if (!params.cnvkit_cnn) {
        ch_cnvkit_reference = file(params.aa_data_repo + "/" + params.reference_build + "/" + params.reference_build + "_cnvkit_filtered_ref.cnn", checkIfExists: true)
    } else {
        ch_cnvkit_reference = file(params.cnvkit_cnn)
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK           } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS & LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


// CONCATENATE FASTQ
include { CAT_FASTQ     }     from '../modules/nf-core/cat/fastq/main'

// QUALITY CONTROL
include { FASTQC        }     from '../modules/nf-core/fastqc/main'

// TRIMMING
include { TRIMGALORE    }    from '../modules/nf-core/trimgalore/main'

// Genome Preparation
include { BWA_INDEX     }   from '../modules/nf-core/bwa/index/main'

// Alignment
include { BWA_MEM                                   }   from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_BAM        }   from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BAM      }   from '../modules/nf-core/samtools/index/main'

// PICARD
include { SAMTOOLS_FAIDX                            }   from '../modules/nf-core/samtools/faidx/main'
include { BAM_MARKDUPLICATES_PICARD                 }   from '../subworkflows/local/bam_markduplicates_picard/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FILTER     }   from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FILTERED   }   from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED }   from '../modules/nf-core/samtools/index/main'

// SAMTOOLS STATISTICS
include { SAMTOOLS_STATS                            }   from '../modules/nf-core/samtools/stats/main'

// BAM STATS
include { BAM_STATS_SAMTOOLS as BAM_STATS_SAMTOOLS_RAW }   from '../subworkflows/local/bam_stats_samtools/main'

// CIRCLE-MAP
include { CIRCLEMAP_READEXTRACTOR                   }   from '../modules/local/circlemap/readextractor.nf'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_RE         }   from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RE       }   from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QNAME_CM   }   from '../modules/nf-core/samtools/sort/main'
include { CIRCLEMAP_REALIGN                         }   from '../modules/local/circlemap/realign.nf'
include { CIRCLEMAP_REPEATS                         }   from '../modules/local/circlemap/repeats.nf'

// CIRCLE_FINDER
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QNAME_CF   }   from '../modules/nf-core/samtools/sort/main'
include { SAMBLASTER                                }     from '../modules/local/samblaster.nf'
include { BEDTOOLS_SORTEDBAM2BED                    }     from '../modules/local/bedtools/sortedbam2bed.nf'
include { BEDTOOLS_SPLITBAM2BED                     }     from '../modules/local/bedtools/splitbam2bed.nf'
include { CIRCLEFINDER                              }     from '../modules/local/circlefinder.nf'

// CIRCexplorer2
include { CIRCEXPLORER2_PARSE       }     from '../modules/local/circexplorer2/parse.nf'

// AmpliconArchitect
include { CNVKIT_BATCH                              }     from '../modules/local/cnvkit/batch/main.nf'
include { CNVKIT_SEGMENT                            }     from '../modules/local/cnvkit/segment.nf'
include { COLLECT_SEEDS                             }     from '../modules/local/collect_seeds.nf'
include { AMPLIFIED_INTERVALS                       }     from '../modules/local/amplified_intervals.nf'
include { AMPLICONARCHITECT_AMPLICONARCHITECT       }     from '../modules/local/ampliconarchitect/ampliconarchitect.nf'
include { AMPLICONCLASSIFIER_AMPLICONCLASSIFIER      }     from '../modules/local/ampliconclassifier/ampliconclassifier.nf'
include { AMPLICONCLASSIFIER_AMPLICONSIMILARITY      }     from '../modules/local/ampliconclassifier/ampliconsimilarity.nf'
include { AMPLICONCLASSIFIER_MAKEINPUT               }     from '../modules/local/ampliconclassifier/makeinput.nf'
include { AMPLICONCLASSIFIER_MAKERESULTSTABLE        }     from '../modules/local/ampliconclassifier/makeresultstable.nf'

// Unicycler
include { UNICYCLER           }     from '../modules/local/unicycler/main.nf'
include { SEQTK_SEQ           }     from '../modules/local/seqtk/seq.nf'
include { GETCIRCULARREADS    }     from '../modules/local/getcircularreads.nf'
include { MINIMAP2_ALIGN      }     from '../modules/nf-core/minimap2/align/main.nf'


// MULTIQC
include { MULTIQC }     from '../modules/local/multiqc.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CIRCDNA {
    ch_versions = Channel.empty()

    // Define Empty Channels for MultiQC
    ch_samtools_stats           = Channel.empty()
    ch_samtools_flagstat        = Channel.empty()
    ch_samtools_idxstats        = Channel.empty()
    ch_markduplicates_stats     = Channel.empty()
    ch_markduplicates_flagstat  = Channel.empty()
    ch_markduplicates_idxstats  = Channel.empty()
    ch_markduplicates_multiqc   = Channel.empty()

    // Check file format
    if (params.input_format == "FASTQ") {
        //
        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        //
        INPUT_CHECK (
            ch_input
        )
        .reads
        .map {
            meta, fastq ->
                meta.id = meta.id.split('_')[0..-2].join('_')
                [ meta, fastq ] }
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

        //
        // MODULE: Concatenate FASTQs from the same samples
        //
        CAT_FASTQ (
            ch_fastq.multiple
        )
        .reads
        .mix(ch_fastq.single)
        .set { ch_cat_fastq }

        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)


        //
        // MODULE: Run FastQC
        //
        ch_fastqc_multiqc = Channel.empty()
        if ( ! params.skip_qc ) {
            FASTQC (
                ch_cat_fastq
            )
            ch_versions         = ch_versions.mix(FASTQC.out.versions)
            ch_fastqc_multiqc   = FASTQC.out.zip
        }

        //
        // MODULE: Run trimgalore
        //
        if ( ! params.skip_trimming ) {
            TRIMGALORE (
                ch_cat_fastq
            )
            ch_trimmed_reads            = TRIMGALORE.out.reads
            ch_trimgalore_multiqc       = TRIMGALORE.out.zip
            ch_trimgalore_multiqc_log   = TRIMGALORE.out.log
            ch_versions                 = ch_versions.mix(TRIMGALORE.out.versions)
        } else {
            ch_trimmed_reads            = INPUT_CHECK.out.reads
            ch_trimgalore_multiqc       = Channel.empty()
            ch_trimgalore_multiqc_log   = Channel.empty()
        }

        //
        // MODULE: Run bwa index
        //
        if (!bwa_index_exists & (run_ampliconarchitect | run_circexplorer2 |
                                run_circle_finder | run_circle_map_realign |
                                run_circle_map_repeats)) {
            BWA_INDEX (
                ch_fasta_meta
            )
            ch_bwa_index = BWA_INDEX.out.index.map{ meta, index -> ["bwa_index", index] }.collect()
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
        }


        //
        // MODULE: BWA MEM ALIGNMENT
        //
        if (run_ampliconarchitect | run_circexplorer2 | run_circle_finder |
            run_circle_map_realign | run_circle_map_repeats) {
            BWA_MEM (
                ch_trimmed_reads,
                ch_bwa_index,
                Channel.value(true)
            )
            ch_bam_sorted   = BWA_MEM.out.bam
            ch_full_bam_sorted   = BWA_MEM.out.bam
            ch_bwa_sorted   = BWA_MEM.out.bam
            ch_versions = ch_versions.mix(BWA_MEM.out.versions)

            // SAMTOOLS INDEX SORTED BAM
            SAMTOOLS_INDEX_BAM (
                ch_bam_sorted
            )
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_BAM.out.versions)
        }
    } else if (params.input_format == "BAM") {
    // Use BAM Files as input
        INPUT_CHECK (
            ch_input
        )
        if (!params.bam_sorted){
            SAMTOOLS_SORT_BAM (
                INPUT_CHECK.out.reads
            )
            ch_versions         = ch_versions.mix(SAMTOOLS_SORT_BAM.out.versions)
            ch_bam_sorted       = SAMTOOLS_SORT_BAM.out.bam
        } else {
            ch_bam_sorted       = INPUT_CHECK.out.reads
            ch_full_bam_sorted  = INPUT_CHECK.out.reads
            ch_bwa_sorted       = INPUT_CHECK.out.reads
        }
        // SAMTOOLS INDEX SORTED BAM
        SAMTOOLS_INDEX_BAM (
            ch_bam_sorted
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_BAM.out.versions)
        ch_fastqc_multiqc           = Channel.empty()
        ch_trimgalore_multiqc       = Channel.empty()
        ch_trimgalore_multiqc_log   = Channel.empty()
    }



    // Define Index channel and additional bam sorted channels for Circle_finder - not usable with duplicates removed
    ch_bam_sorted_bai       = SAMTOOLS_INDEX_BAM.out.bai
    ch_full_bam_sorted      = ch_bam_sorted
    ch_full_bam_sorted_bai  = SAMTOOLS_INDEX_BAM.out.bai

    if (run_ampliconarchitect | run_circexplorer2 | run_circle_finder |
        run_circle_map_realign | run_circle_map_repeats) {

        ch_fasta = ch_fasta_meta.map{ meta, index -> [index] }.collect()
        BAM_STATS_SAMTOOLS_RAW (
            ch_bam_sorted.join(ch_bam_sorted_bai).
                map { meta, bam, bai -> [meta, bam, bai] },
                ch_fasta
        )
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS_RAW.out.versions)
        ch_samtools_stats               = BAM_STATS_SAMTOOLS_RAW.out.stats
        ch_samtools_flagstat            = BAM_STATS_SAMTOOLS_RAW.out.flagstat
        ch_samtools_idxstats            = BAM_STATS_SAMTOOLS_RAW.out.idxstats

        // PICARD MARK_DUPLICATES
        if (!params.skip_markduplicates) {
            // Index Fasta File for Markduplicates
            SAMTOOLS_FAIDX ( ch_fasta_meta )
            ch_fai = SAMTOOLS_FAIDX.out.fai.map {meta, fai -> fai }.collect()

            // MARK DUPLICATES IN BAM FILE
            BAM_MARKDUPLICATES_PICARD (
                ch_bam_sorted, ch_fasta, ch_fai
            )

            // FILTER DUPLICATES IN BAM FILES USING SAMTOOLS VIEW
            if (!params.keep_duplicates) {
                SAMTOOLS_VIEW_FILTER (
                    ch_bam_sorted.join(ch_bam_sorted_bai),
                    ch_fasta,
                    []
                )
                ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FILTER.out.versions)

                // SORT FILTERED BAM FILE
                SAMTOOLS_SORT_FILTERED (
                    SAMTOOLS_VIEW_FILTER.out.bam
                )
                ch_versions = ch_versions.mix(SAMTOOLS_SORT_FILTERED.out.versions)

                // INDEX FILTERED BAM FILE
                SAMTOOLS_INDEX_FILTERED (
                    SAMTOOLS_SORT_FILTERED.out.bam
                )

                ch_bam_sorted = SAMTOOLS_SORT_FILTERED.out.bam
                ch_bam_sorted_bai = SAMTOOLS_INDEX_FILTERED.out.bai
                ch_versions = ch_versions.mix(SAMTOOLS_INDEX_FILTERED.out.versions)
            }
            else {
                ch_bam_sorted               = BAM_MARKDUPLICATES_PICARD.out.bam
                ch_bam_sorted_bai           = BAM_MARKDUPLICATES_PICARD.out.bai
                ch_markduplicates_stats     = BAM_MARKDUPLICATES_PICARD.out.stats
                ch_markduplicates_flagstat  = BAM_MARKDUPLICATES_PICARD.out.flagstat
                ch_markduplicates_idxstats  = BAM_MARKDUPLICATES_PICARD.out.idxstats
                ch_markduplicates_multiqc   = BAM_MARKDUPLICATES_PICARD.out.metrics
                ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
            }
        } else {
                ch_markduplicates_stats         = Channel.empty()
                ch_markduplicates_flagstat      = Channel.empty()
                ch_markduplicates_idxstats      = Channel.empty()
                ch_markduplicates_multiqc       = Channel.empty()
        }
    }

    if (run_ampliconarchitect) {
        CNVKIT_BATCH (
            ch_bam_sorted.join(ch_bam_sorted_bai),
            ch_fasta,
            ch_cnvkit_reference
        )
        ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)

        CNVKIT_SEGMENT (
            CNVKIT_BATCH.out.cnr
        )
        ch_versions = ch_versions.mix(CNVKIT_SEGMENT.out.versions)

        COLLECT_SEEDS (
            CNVKIT_SEGMENT.out.cns
        )
        ch_versions = ch_versions.mix(COLLECT_SEEDS.out.versions)

        ch_aa_seeds = COLLECT_SEEDS.out.bed
        AMPLIFIED_INTERVALS (
            ch_aa_seeds.join(ch_bam_sorted).join(ch_bam_sorted_bai)
        )
        ch_versions = ch_versions.mix(AMPLIFIED_INTERVALS.out.versions)

        AMPLICONARCHITECT_AMPLICONARCHITECT (
            ch_bam_sorted.join(ch_bam_sorted_bai).
                join(AMPLIFIED_INTERVALS.out.bed)
        )
        ch_versions = ch_versions.mix(AMPLICONARCHITECT_AMPLICONARCHITECT.out.versions)

        ch_aa_cycles = AMPLICONARCHITECT_AMPLICONARCHITECT.out.cycles
        ch_aa_graphs = AMPLICONARCHITECT_AMPLICONARCHITECT.out.graph

//        AMPLICONCLASSIFIER_MAKEINPUT (
//            ch_aa_cycles,
//            ch_aa_graphs
//        )

//        AMPLICONCLASSIFIER_AMPLICONCLASSIFIER (
//            ch_aa_cycles.join(ch_aa_graphs)
//        )
//        ch_versions = ch_versions.mix(AMPLICONCLASSIFIER_AMPLICONCLASSIFIER.out.versions)
//        AMPLICONCLASSIFIER_AMPLICONSIMILARITY (
//            ch_aa_cycles.join(ch_aa_graphs)
//        )
//        aa_summary_ch = AMPLICONARCHITECT_AMPLICONARCHITECT.out.summary
//        ch_versions = ch_versions.mix(AMPLICONCLASSIFIER_AMPLICONSIMILARITY.out.versions)

//        AMPLICONCLASSIFIER_MAKERESULTSTABLE (
//            aa_summary_ch.join(AMPLICONCLASSIFIER_AMPLICONCLASSIFIER.out.class_tsv)
//        )
//        ch_versions = ch_versions.mix(SUMMARISE_AA.out.versions)
    }


    //
    // SUBWORKFLOW - RUN CIRCLE_FINDER PIPELINE
    //
    if (run_circle_finder) {
        SAMTOOLS_SORT_QNAME_CF (
            ch_full_bam_sorted
        )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_QNAME_CF.out.versions)

        SAMBLASTER (
            SAMTOOLS_SORT_QNAME_CF.out.bam
        )
        ch_versions = ch_versions.mix(SAMBLASTER.out.versions)

        BEDTOOLS_SPLITBAM2BED (
            SAMBLASTER.out.split_bam
        )
        ch_versions = ch_versions.mix(BEDTOOLS_SPLITBAM2BED.out.versions)

        BEDTOOLS_SORTEDBAM2BED (
            ch_full_bam_sorted.join(ch_full_bam_sorted_bai)
        )
        ch_versions = ch_versions.mix(BEDTOOLS_SORTEDBAM2BED.out.versions)

        ch_b2b_sorted = BEDTOOLS_SORTEDBAM2BED.out.conc_txt
        ch_b2b_split = BEDTOOLS_SPLITBAM2BED.out.split_txt
        CIRCLEFINDER (
            ch_b2b_split.join(ch_b2b_sorted)
        )
    }

    //
    // SUBWORKFLOW: RUN CIRCLE-MAP REALIGN or REPEATS PIPELINE
    //

    if (run_circle_map_realign ||
            run_circle_map_repeats) {
        ch_bam_sorted_bai.view()
        ch_bam_sorted.view()
        SAMTOOLS_SORT_QNAME_CM (
            ch_bam_sorted
        )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_QNAME_CM.out.versions)

        CIRCLEMAP_READEXTRACTOR (
            SAMTOOLS_SORT_QNAME_CM.out.bam
        )
        ch_versions = ch_versions.mix(CIRCLEMAP_READEXTRACTOR.out.versions)

        SAMTOOLS_SORT_RE (
            CIRCLEMAP_READEXTRACTOR.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_RE.out.versions)

        SAMTOOLS_INDEX_RE (
            SAMTOOLS_SORT_RE.out.bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_RE.out.versions)

        // DEFINE CHANNELS FOR REALIGN AND REPEATS
        ch_qname_sorted_bam = SAMTOOLS_SORT_QNAME_CM.out.bam
        ch_re_sorted_bam = SAMTOOLS_SORT_RE.out.bam
        ch_re_sorted_bai = SAMTOOLS_INDEX_RE.out.bai

        //
        // MODULE: RUN CIRCLE_MAP REPEATS
        //
        if (run_circle_map_repeats) {
            CIRCLEMAP_REPEATS (
                ch_re_sorted_bam.join(ch_re_sorted_bai)
            )
            ch_versions = ch_versions.mix(CIRCLEMAP_REPEATS.out.versions)
        }

        //
        // MODULE: Run Circle-Map Realign
        //
        if (run_circle_map_realign) {

            CIRCLEMAP_REALIGN (
                ch_re_sorted_bam.join(ch_re_sorted_bai).
                    join(ch_qname_sorted_bam).
                    join(ch_bam_sorted).
                    join(ch_bam_sorted_bai),
                ch_fasta
            )
            ch_versions = ch_versions.mix(CIRCLEMAP_REALIGN.out.versions)
        }
    }


    if (run_circexplorer2) {
        CIRCEXPLORER2_PARSE (
            ch_bam_sorted.join(ch_bam_sorted_bai)
        )
        ch_versions = ch_versions.mix(CIRCEXPLORER2_PARSE.out.versions)
    }

    if (run_unicycler && params.input_format == "FASTQ") {

        UNICYCLER (
            ch_trimmed_reads
        )
        ch_versions = ch_versions.mix(UNICYCLER.out.versions)

        SEQTK_SEQ (
            UNICYCLER.out.scaffolds
        )
        ch_versions = ch_versions.mix(SEQTK_SEQ.out.versions)

        GETCIRCULARREADS (
            SEQTK_SEQ.out.fastq
        )

        GETCIRCULARREADS.out.fastq
            .map { meta, fastq -> [ meta + [single_end: true], fastq ] }
            .set { ch_circular_fastq }

        MINIMAP2_ALIGN (
            ch_circular_fastq,
            ch_fasta,
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    }
    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowCircdna.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_fastqc_multiqc.collect{it[1]}.ifEmpty([]),
            ch_trimgalore_multiqc.collect{it[1]}.ifEmpty([]),
            ch_trimgalore_multiqc_log.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_flagstat.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_stats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([]),
        )
        multiqc_report       = MULTIQC.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
