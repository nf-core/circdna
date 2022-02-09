/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCircleseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Fasta reference genome not specified!' }

if (params.circle_identifier != "circle_map_realign" &
    params.circle_identifier != "circle_map_repeats" &
    params.circle_identifier != "circle_finder" &
    params.circle_identifier != "circexplorer2" &
    params.circle_identifier != "ampliconarchitect" &
    params.circle_identifier != "all" ) {exit 1, 'Circle Identifier Software/Algorithm not specified!' }

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
if (params.circle_identifier == "ampliconarchitect") {
    if (!params.mosek_license_dir) { exit 1, "Mosek Missing" }
    if (!params.aa_data_repo) { exit 1, "AmpliconArchitect Data Repository Missing" }
    if (params.reference_build != "hg19" & params.reference_build != "GRCh38" & params.reference_build != "GRCh37"){
        exit 1, "Reference Build not given! Please specify --reference_build 'hg19', 'GRCh38', or 'GRCh37'."
    }
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

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK           } from '../subworkflows/local/input_check'                  addParams( options: [:] )
include { BAM_STATS_SAMTOOLS    } from '../subworkflows/nf-core/bam_stats_samtools/main'    addParams( options: modules['samtools_stats_options'] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS & LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

// CONCATENATE FASTQ
include { CAT_FASTQ }     from '../modules/nf-core/modules/cat/fastq/main'       addParams( options: modules['cat_fastq']   )

// QUALITY CONTROL
include { FASTQC as FASTQC_RAW  }     from '../modules/nf-core/modules/fastqc/main'       addParams( options: modules['fastqc']   )
include { FASTQC as FASTQC_TRIM }     from '../modules/nf-core/modules/fastqc/main'       addParams( options: modules['fastqc']   )

// TRIMMING
include { TRIMGALORE }    from '../modules/nf-core/modules/trimgalore/main'     addParams( options: modules['trimgalore'] )

// ALIGNMENT
include { BWA_INDEX }   from '../modules/nf-core/modules/bwa/index/main'    addParams( options: modules['bwa_index'])
include { BWA_MEM   }   from '../modules/nf-core/modules/bwa/mem/main'      addParams( options: modules['bwa_mem']  )

// PICARD
include { PICARD_MARKDUPLICATES     }   from '../modules/nf-core/modules/picard/markduplicates/main'      addParams( options: modules['picard_markduplicates']  )
include { MARK_DUPLICATES_PICARD    } from '../subworkflows/nf-core/mark_duplicates_picard'     addParams( markduplicates_options: modules['picard_markduplicates'], samtools_index_options: modules['picard_markduplicates_samtools'], samtools_stats_options:  modules['picard_markduplicates_samtools'] )


// SAMTOOLS SORT & INDEX
include { SAMTOOLS_FAIDX                            }   from '../modules/nf-core/modules/samtools/faidx/main'   addParams( options: modules['samtools_faidx']       )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BAM      }   from '../modules/nf-core/modules/samtools/index/main'   addParams( options: modules['samtools_index_bwa']   )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RE       }   from '../modules/nf-core/modules/samtools/index/main'   addParams( options: modules['samtools_index_re']    )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED }   from '../modules/nf-core/modules/samtools/index/main'   addParams( options: modules['samtools_index_re']    )
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QNAME_CM   }   from '../modules/nf-core/modules/samtools/sort/main'    addParams( options: modules['samtools_sort_qname']  )
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QNAME_CF   }   from '../modules/nf-core/modules/samtools/sort/main'    addParams( options: modules['samtools_sort_qname']  )
include { SAMTOOLS_SORT as SAMTOOLS_SORT_RE         }   from '../modules/nf-core/modules/samtools/sort/main'    addParams( options: modules['samtools_sort_re']     )
include { SAMTOOLS_SORT                             }   from '../modules/nf-core/modules/samtools/sort/main'    addParams( options: modules['samtools_sort']        )
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FILTERED   }   from '../modules/nf-core/modules/samtools/sort/main'    addParams( options: modules['samtools_sort_filtered']        )


// FILTER BAM FILE USING SAMTOOLS VIEW
def samtools_view_filter_options     = modules['samtools_view_filter']
samtools_view_filter_options.args   += params.keep_duplicates ? '' : Utils.joinModuleArgs(['-F 0x0400'])
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FILTER }   from '../modules/nf-core/modules/samtools/view/main'    addParams( options: samtools_view_filter_options )

// SAMTOOLS STATISTICS
include { SAMTOOLS_STATS                        }   from '../modules/nf-core/modules/samtools/stats/main'    addParams( options: modules['samtools_stats']   )

// CIRCLE-MAP
include { CIRCLEMAP_READEXTRACTOR   }   from '../modules/local/circlemap/readextractor.nf'  addParams( options: modules['circlemap_readextractor']  )
include { CIRCLEMAP_REALIGN         }   from '../modules/local/circlemap/realign.nf'        addParams( options: modules['circlemap_realign']        )
include { CIRCLEMAP_REPEATS         }   from '../modules/local/circlemap/repeats.nf'        addParams( options: modules['circlemap_repeats']        )

// CIRCLE_FINDER
include { SAMBLASTER                }     from '../modules/local/samblaster.nf'                 addParams( options: modules['samblaster']               )
include { BEDTOOLS_SORTEDBAM2BED    }     from '../modules/local/bedtools/sortedbam2bed.nf'     addParams( options: modules['bedtools_sortedbam2bed']   )
include { BEDTOOLS_SPLITBAM2BED     }     from '../modules/local/bedtools/splitbam2bed.nf'      addParams( options: modules['bedtools_splitbam2bed']    )
include { CIRCLEFINDER              }     from '../modules/local/circlefinder.nf'               addParams( options: modules['circlefinder']             )

// CIRCexplorer2
include { CIRCEXPLORER2_PARSE       }     from '../modules/local/circexplorer2/parse.nf'               addParams( options: modules['circexplorer2_parse']             )

// AmpliconArchitect
include { CNVKIT_BATCH                          }     from '../modules/nf-core/modules/cnvkit/batch/main.nf'            addParams( options: modules['cnvkit_batch']          )
include { CNVKIT_SEGMENT                        }     from '../modules/local/cnvkit/segment.nf'            addParams( options: modules['cnvkit_batch']          )
include { AMPLICONARCHITECT_PREPAREAA           }     from '../modules/local/ampliconarchitect/prepareaa.nf'            addParams( options: modules['ampliconarchitect_prepareaa']          )
include { AMPLICONARCHITECT_AMPLICONARCHITECT   }     from '../modules/local/ampliconarchitect/ampliconarchitect.nf'    addParams( options: modules['ampliconarchitect_ampliconarchitect']  )

// Unicycler
include { UNICYCLER           }     from '../modules/nf-core/modules/unicycler/main.nf'            addParams( options: modules['unicycler']          )
include { SEQTK_SEQ           }     from '../modules/local/seqtk/seq.nf'                           addParams( options: modules['seqtk_seq']          )
include { MINIMAP2_ALIGN      }     from '../modules/nf-core/modules/minimap2/align/main.nf'                           addParams( options: modules['minimap2_align']          )


// MULTIQC
include { MULTIQC }     from '../modules/nf-core/modules/multiqc/main'      addParams( options: multiqc_options     )


// TRIMGALORE OPTIONS
def trimgalore_options    = modules['trimgalore']
trimgalore_options.args  += params.trim_nextseq > 0 ? Utils.joinModuleArgs(["--nextseq ${params.trim_nextseq}"]) : ''
if (params.save_trimmed)  { trimgalore_options.publish_files.put('fq.gz','') }

// CONCATENATE FASTQ OPTIONS
def cat_fastq_options          = modules['cat_fastq']
if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

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
            FASTQC_RAW (
                ch_cat_fastq
            )
            ch_software_versions = ch_software_versions.mix(FASTQC_RAW.out.version.first().ifEmpty(null))
            ch_fastqc_report = FASTQC_RAW.out.zip
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
        // MODULE: Run samtools faidx
        //
        SAMTOOLS_FAIDX (
            ch_fasta
        )

        //
        // MODULE: BWA MEM ALIGNMENT
        //
        BWA_MEM (
            ch_trimmed_reads,
            ch_bwa_index
        )
        ch_bwa_sorted_bam = BWA_MEM.out.sorted_bam

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
            SAMTOOLS_SORT (
                INPUT_CHECK.out.reads
            )
            ch_bwa_sorted_bam = SAMTOOLS_SORT.out.bam
        } else {
            ch_bwa_sorted_bam = INPUT_CHECK.out.reads
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


    if (params.circle_identifier == "ampliconarchitect" || params.circle_identifier == "all") {
        ch_cnvkit_targets = params.cnvkit_targets ? Channel.fromPath(params.cnvkit_targets) : Channel.empty()
        ch_cnvkit_reference = params.cnvkit_reference ? Channel.fromPath(params.cnvkit_reference) : Channel.empty()

        CNVKIT_BATCH (
            ch_bwa_sorted_bam.join(ch_bwa_sorted_bai),
            ch_fasta
        )

        CNVKIT_SEGMENT (
            CNVKIT_BATCH.out.cnr
        )

        AMPLICONARCHITECT_PREPAREAA (
            ch_bwa_sorted_bam.join(ch_bwa_sorted_bai),
            CNVKIT_SEGMENT.out.cns
        )
        ch_prepareaa_bed = AMPLICONARCHITECT_PREPAREAA.out.bed
        AMPLICONARCHITECT_AMPLICONARCHITECT (
            ch_bwa_sorted_bam.join(ch_bwa_sorted_bai).
                join(ch_prepareaa_bed)
        )
    }

    //
    // SUBWORKFLOW - RUN CIRCLE_FINDER PIPELINE
    //
    if (params.circle_identifier == "circle_finder" || params.circle_identifier == "all") {
        SAMTOOLS_SORT_QNAME_CF (
            ch_bwa_sorted_bam
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

    if (params.circle_identifier == "circle_map_realign" ||
            params.circle_identifier == "circle_map_repeats" || params.circle_identifier == "all") {
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
        if (params.circle_identifier == "circle_map_repeats" || params.circle_identifier == "all") {
            CIRCLEMAP_REPEATS (
                ch_re_sorted_bam.join(ch_re_sorted_bai)
            )
        }

        //
        // MODULE: Run Circle-Map Realign
        //
        if (params.circle_identifier == "circle_map_realign" || params.circle_identifier == "all") {
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


    if (params.circle_identifier == "circexplorer2" || params.circle_identifier == "all") {
        CIRCEXPLORER2_PARSE (
            ch_bwa_sorted_bam.join(ch_bwa_sorted_bai)
        )
    }

    if (!params.skip_denovo_assembly && params.input_format == "FASTQ") {
        ch_trimmed_reads.map{ meta, file -> [meta, file, []] }.set{ch_unicycler_input}
        UNICYCLER (
            ch_unicycler_input
        )
        SEQTK_SEQ (
            UNICYCLER.out.scaffolds
        )
        SEQTK_SEQ.out.fastq_circular.view()
        MINIMAP2_ALIGN (
            ch_fasta,
            SEQTK_SEQ.out.fastq_circular
        )
    }

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCircleseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_trimgalore_report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_bam_flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_bam_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_bam_idxstats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
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
