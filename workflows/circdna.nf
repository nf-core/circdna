/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCircdna.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Fasta reference genome not specified!' }

software = params.circle_identifier.split(",")
run_circexplorer2 = ("circexplorer2" in software)
run_circle_map_realign = ("circle_map_realign" in software)
run_circle_map_repeats = ("circle_map_repeats" in software)
run_circle_finder = ("circle_finder" in software)
run_ampliconarchitect = ("ampliconarchitect" in software)
run_unicycler = ("unicycler" in software)

if (!(run_unicycler | run_circle_map_realign | run_circle_map_repeats | run_circle_finder | run_ampliconarchitect)) {
    exit 1, 'circle_identifier param not valid. Please check!'
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
    if (!params.mosek_license_dir) { exit 1, "Mosek Missing" }
    if (!params.aa_data_repo) { exit 1, "AmpliconArchitect Data Repository Missing" }
    if (params.reference_build != "hg19" & params.reference_build != "GRCh38" & params.reference_build != "GRCh37"){
        exit 1, "Reference Build not given! Please specify --reference_build 'hg19', 'GRCh38', or 'GRCh37'."
    }
    ch_cnvkit_reference = Channel.fromPath(params.aa_data_repo + "/" + params.reference_build + "/" + params.reference_build + "_cnvkit_filtered_ref.cnn")
}


/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK           } from '../subworkflows/local/input_check'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS & LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'


// CONCATENATE FASTQ
include { CAT_FASTQ }     from '../modules/nf-core/modules/cat/fastq/main'

// QUALITY CONTROL
include { FASTQC }     from '../modules/nf-core/modules/fastqc/main'

// TRIMMING
include { TRIMGALORE }    from '../modules/nf-core/modules/trimgalore/main'

// Genome Preparation
include { BWA_INDEX }   from '../modules/nf-core/modules/bwa/index/main'

// Alignment
include { BWA_MEM   }   from '../modules/nf-core/modules/bwa/mem/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_BAM        }   from '../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BAM      }   from '../modules/nf-core/modules/samtools/index/main'

// PICARD
include { MARK_DUPLICATES_PICARD    } from '../subworkflows/nf-core/mark_duplicates_picard'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FILTER }   from '../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FILTERED   }   from '../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED }   from '../modules/nf-core/modules/samtools/index/main'

// SAMTOOLS SORT & INDEX
include { SAMTOOLS_FAIDX                            }   from '../modules/nf-core/modules/samtools/faidx/main'


// SAMTOOLS STATISTICS
include { SAMTOOLS_STATS                        }   from '../modules/nf-core/modules/samtools/stats/main'

// CIRCLE-MAP
include { CIRCLEMAP_READEXTRACTOR   }   from '../modules/local/circlemap/readextractor.nf'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_RE         }   from '../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RE       }   from '../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QNAME_CM   }   from '../modules/nf-core/modules/samtools/sort/main'
include { CIRCLEMAP_REALIGN         }   from '../modules/local/circlemap/realign.nf'
include { CIRCLEMAP_REPEATS         }   from '../modules/local/circlemap/repeats.nf'

// CIRCLE_FINDER
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QNAME_CF   }   from '../modules/nf-core/modules/samtools/sort/main'
include { SAMBLASTER                }     from '../modules/local/samblaster.nf'
include { BEDTOOLS_SORTEDBAM2BED    }     from '../modules/local/bedtools/sortedbam2bed.nf'
include { BEDTOOLS_SPLITBAM2BED     }     from '../modules/local/bedtools/splitbam2bed.nf'
include { CIRCLEFINDER              }     from '../modules/local/circlefinder.nf'

// CIRCexplorer2
include { CIRCEXPLORER2_PARSE       }     from '../modules/local/circexplorer2/parse.nf'

// AmpliconArchitect
include { CNVKIT_BATCH                              }     from '../modules/nf-core/modules/cnvkit/batch/main.nf'
include { CNVKIT_SEGMENT                            }     from '../modules/local/cnvkit/segment.nf'
include { AMPLICONARCHITECT_PREPAREAA               }     from '../modules/local/ampliconarchitect/prepareaa.nf'
include { AMPLICONARCHITECT_AMPLICONARCHITECT       }     from '../modules/local/ampliconarchitect/ampliconarchitect.nf'
include { AMPLICONARCHITECT_AMPLICONCLASSIFIER      }     from '../modules/local/ampliconarchitect/ampliconclassifier.nf'

// Unicycler
include { UNICYCLER           }     from '../modules/nf-core/modules/unicycler/main.nf'
include { SEQTK_SEQ           }     from '../modules/local/seqtk/seq.nf'
include { GETCIRCULARREADS    }     from '../modules/local/getcircularreads.nf'
include { MINIMAP2_ALIGN      }     from '../modules/nf-core/modules/minimap2/align/main.nf'


// MULTIQC
include { MULTIQC }     from '../modules/nf-core/modules/multiqc/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CIRCDNA {

    ch_software_versions = Channel.empty()

    // Define channels of fastqc, trimgalore, bwa stats for multiqc
    ch_fastqc_report     = Channel.empty()
    ch_trimgalore_report = Channel.empty()

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
        // ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

        CAT_FASTQ (
            ch_fastq.multiple
        )
        .reads
        .mix(ch_fastq.single)
        .set { ch_cat_fastq }


        //
        // MODULE: Run FastQC
        //
        if ( ! params.skip_qc ) {
            FASTQC (
                ch_cat_fastq
            )
            ch_software_versions = ch_software_versions.mix(FASTQC.out.versions.first().ifEmpty(null))
            ch_fastqc_report = FASTQC.out.zip
        }

        //
        // MODULE: Run trimgalore
        //
        if ( ! params.skip_trimming ) {
            TRIMGALORE (
                ch_cat_fastq
            )
            ch_trimmed_reads        = TRIMGALORE.out.reads
            ch_trimgalore_report    = TRIMGALORE.out.zip
        } else {
            ch_trimmed_reads = INPUT_CHECK.out.reads
        }

        //
        // MODULE: Run bwa index
        //
        if (!bwa_index_exists) {
            BWA_INDEX (
                ch_fasta
            )
            ch_bwa_index = BWA_INDEX.out.index
        }


        //
        // MODULE: BWA MEM ALIGNMENT
        //
        BWA_MEM (
            ch_trimmed_reads,
            ch_bwa_index,
            true
        )
        ch_bwa_sorted_bam = BWA_MEM.out.bam
        ch_bam_sorted     = BWA_MEM.out.bam

        // SAMTOOLS INDEX SORTED BAM
        SAMTOOLS_INDEX_BAM(
            ch_bwa_sorted_bam
        )

    } else if (params.input_format == "BAM") {
    // Use BAM Files as input
        INPUT_CHECK (
            ch_input
        )
        if (!params.bam_sorted){
            SAMTOOLS_SORT_BAM (
                INPUT_CHECK.out.reads
            )
            ch_bam_sorted       = SAMTOOLS_SORT_BAM.out.bam
            ch_bwa_sorted_bam   = SAMTOOLS_SORT_BAM.out.bam
        } else {
            ch_bwa_sorted_bam   = INPUT_CHECK.out.reads
            ch_bam_sorted       = INPUT_CHECK.out.reads
            // SAMTOOLS INDEX SORTED BAM
        }
        SAMTOOLS_INDEX_BAM (
            ch_bwa_sorted_bam
        )
    }
    ch_bwa_sorted_bai       = SAMTOOLS_INDEX_BAM.out.bai

    // PICARD MARK_DUPLICATES
    if (!params.skip_markduplicates) {
        MARK_DUPLICATES_PICARD (
            ch_bwa_sorted_bam
        )
        ch_bwa_sorted_bam         = MARK_DUPLICATES_PICARD.out.bam
        ch_bwa_sorted_bai         = MARK_DUPLICATES_PICARD.out.bai
        ch_bam_stats              = MARK_DUPLICATES_PICARD.out.stats
        ch_bam_flagstat           = MARK_DUPLICATES_PICARD.out.flagstat
        ch_bam_idxstats           = MARK_DUPLICATES_PICARD.out.idxstats
        ch_markduplicates_multiqc = MARK_DUPLICATES_PICARD.out.metrics
        // FILTER BAM FILES USING SAMTOOLS VIEW
        SAMTOOLS_VIEW_FILTER (
            ch_bwa_sorted_bam, ch_fasta
        )
        SAMTOOLS_SORT_FILTERED (
            SAMTOOLS_VIEW_FILTER.out.bam
        )
        ch_bwa_sorted_bam = SAMTOOLS_SORT_FILTERED.out.bam
        SAMTOOLS_INDEX_FILTERED (
            ch_bwa_sorted_bam
        )
        ch_bwa_sorted_bai = SAMTOOLS_INDEX_FILTERED.out.bai
    } else {
        ch_bam_stats              = Channel.empty()
        ch_bam_flagstat           = Channel.empty()
        ch_bam_idxstats           = Channel.empty()
        ch_markduplicates_multiqc = Channel.empty()
    }

    if (run_ampliconarchitect) {
        CNVKIT_BATCH (
            ch_bwa_sorted_bam.join(ch_bwa_sorted_bai),
            ch_fasta,
            ch_cnvkit_reference
        )

        CNVKIT_SEGMENT (
            CNVKIT_BATCH.out.cnr
        )

        AMPLICONARCHITECT_PREPAREAA (
            ch_bwa_sorted_bam.join(ch_bwa_sorted_bai).
            join(CNVKIT_SEGMENT.out.cns)
        )
        ch_prepareaa_bed = AMPLICONARCHITECT_PREPAREAA.out.bed
        AMPLICONARCHITECT_AMPLICONARCHITECT (
            ch_bwa_sorted_bam.join(ch_bwa_sorted_bai).
                join(ch_prepareaa_bed)
        )
        ch_aa_cycles = AMPLICONARCHITECT_AMPLICONARCHITECT.out.cycles
        ch_aa_graphs = AMPLICONARCHITECT_AMPLICONARCHITECT.out.graph
        AMPLICONARCHITECT_AMPLICONCLASSIFIER (
            ch_aa_cycles.join(ch_aa_graphs)
        )
    }


    //
    // SUBWORKFLOW - RUN CIRCLE_FINDER PIPELINE
    //
    if (run_circle_finder) {
        SAMTOOLS_SORT_QNAME_CF (
            ch_bam_sorted
        )
        SAMBLASTER (
            SAMTOOLS_SORT_QNAME_CF.out.bam
        )

        BEDTOOLS_SPLITBAM2BED (
            SAMBLASTER.out.split_bam
        )
        BEDTOOLS_SORTEDBAM2BED (
            ch_bwa_sorted_bam.join(ch_bwa_sorted_bai)
        )

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
        //
        // MODULE: Run samtools faidx
        //
        SAMTOOLS_FAIDX (
            ch_fasta
        )
        SAMTOOLS_SORT_QNAME_CM (
            ch_bwa_sorted_bam
        )
        CIRCLEMAP_READEXTRACTOR (
            SAMTOOLS_SORT_QNAME_CM.out.bam
        )

        SAMTOOLS_SORT_RE (
            CIRCLEMAP_READEXTRACTOR.out.circ_bam
        )

        SAMTOOLS_INDEX_RE (
            SAMTOOLS_SORT_RE.out.bam
        )


        // DEFINE CHANNELS FOR REALIGN AND REPEATS
        ch_qname_sorted_bam = SAMTOOLS_SORT_QNAME_CM.out.bam
        ch_re_sorted_bam = SAMTOOLS_SORT_RE.out.bam
        ch_re_sorted_bai = SAMTOOLS_INDEX_RE.out.bai
        ch_fasta_index = SAMTOOLS_FAIDX.out.fai

        //
        // MODULE: RUN CIRCLE_MAP REPEATS
        //
        if (run_circle_map_repeats) {
            CIRCLEMAP_REPEATS (
                ch_re_sorted_bam.join(ch_re_sorted_bai)
            )
        }

        //
        // MODULE: Run Circle-Map Realign
        //
        if (run_circle_map_realign) {
            CIRCLEMAP_REALIGN (
                ch_re_sorted_bam.join(ch_re_sorted_bai).
                    join(ch_qname_sorted_bam).
                    join(ch_bwa_sorted_bam).
                    join(ch_bwa_sorted_bai),
                ch_fasta,
                ch_fasta_index
            )
        }
    }


    if (run_circexplorer2) {
        CIRCEXPLORER2_PARSE (
            ch_bwa_sorted_bam.join(ch_bwa_sorted_bai)
        )
    }

    if (run_unicycler && params.input_format == "FASTQ") {
        UNICYCLER (
            ch_trimmed_reads
        )
        SEQTK_SEQ (
            UNICYCLER.out.scaffolds
        )
        GETCIRCULARREADS (
            SEQTK_SEQ.out.fastq
        )
        MINIMAP2_ALIGN (
            GETCIRCULARREADS.out.fastq,
            ch_fasta
        )
    } else if (run_unicycler && !params.input_format == "FASTQ") {
        exit 1, 'Unicycler needs FastQ input. Please specify input_format == "FASTQ", if possible, or don`t run unicycler.'
    }

    //
    // MODULE: Pipeline reporting
    //
//    CUSTOM_DUMPSOFTWAREVERSIONS (
//        ch_versions.unique().collectFile(name: 'collated_versions.yml')
//    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowCircdna.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

//        MULTIQC (
//            ch_multiqc_files.collect(),
//            ch_multiqc_custom_config.collect().ifEmpty([]),
////  CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
//            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
//
//        )
        //multiqc_report       = MULTIQC.out.report.toList()
    }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
