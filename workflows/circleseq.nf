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
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS & LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

// QUALITY CONTROL
include { FASTQC as FASTQC_RAW  }     from '../modules/nf-core/modules/fastqc/main'       addParams( options: modules['fastqc']   )
include { FASTQC as FASTQC_TRIM }     from '../modules/nf-core/modules/fastqc/main'       addParams( options: modules['fastqc']   )

// TRIMMING
include { TRIMGALORE }    from '../modules/nf-core/modules/trimgalore/main'     addParams( options: modules['trimgalore'] )

// ALIGNMENT
include { BWA_INDEX }   from '../modules/nf-core/modules/bwa/index/main'    addParams( options: modules['bwa_index'])
include { BWA_MEM   }   from '../modules/nf-core/modules/bwa/mem/main'      addParams( options: modules['bwa_mem']  )

// SAMTOOLS SORT & INDEX
include { SAMTOOLS_FAIDX                        }   from '../modules/nf-core/modules/samtools/faidx/main'   addParams( options: modules['samtools_faidx']   )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BWA  }   from '../modules/nf-core/modules/samtools/index/main'   addParams( options: modules['samtools_index_bwa']   )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RE   }   from '../modules/nf-core/modules/samtools/index/main'   addParams( options: modules['samtools_index_re']   )
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QNAME  }   from '../modules/nf-core/modules/samtools/sort/main'    addParams( options: modules['samtools_sort_qname']   )
include { SAMTOOLS_SORT as SAMTOOLS_SORT_RE     }   from '../modules/nf-core/modules/samtools/sort/main'    addParams( options: modules['samtools_sort_re']   )


// CIRCLE-MAP
include { CIRCLEMAP_READEXTRACTOR   }   from '../modules/local/circlemap/readextractor.nf'  addParams( options: modules['circlemap_readextractor']  )
include { CIRCLEMAP_REALIGN         }   from '../modules/local/circlemap/realign.nf'        addParams( options: modules['circlemap_realign']        )
include { CIRCLEMAP_REPEATS         }   from '../modules/local/circlemap/repeats.nf'        addParams( options: modules['circlemap_repeats']        )

// CIRCLE_FINDER
include { SAMBLASTER                }     from '../modules/local/samblaster.nf'                 addParams( options: modules['samblaster']               )
include { BEDTOOLS_SORTEDBAM2BED    }     from '../modules/local/bedtools/sortedbam2bed.nf'     addParams( options: modules['bedtools_sortedbam2bed']   )
include { BEDTOOLS_SPLITBAM2BED     }     from '../modules/local/bedtools/splitbam2bed.nf'      addParams( options: modules['bedtools_splitbam2bed']    )
include { CIRCLEFINDER              }     from '../modules/local/circlefinder.nf'               addParams( options: modules['circlefinder']             )

// MULTIQC
include { MULTIQC }     from '../modules/nf-core/modules/multiqc/main'      addParams( options: multiqc_options     )

// TRIMGALORE OPTIONS
def trimgalore_options    = modules['trimgalore']
trimgalore_options.args  += params.trim_nextseq > 0 ? Utils.joinModuleArgs(["--nextseq ${params.trim_nextseq}"]) : ''
if (params.save_trimmed)  { trimgalore_options.publish_files.put('fq.gz','') }

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CIRCLESEQ {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    //
    // MODULE: Run FastQC
    //

    if ( ! params.skip_qc ) {
        FASTQC_RAW (
            INPUT_CHECK.out.reads
        )
        ch_software_versions = ch_software_versions.mix(FASTQC_RAW.out.version.first().ifEmpty(null))
    }

    //
    // MODULE: Run trimgalore
    //
    if ( ! params.skip_trimming ) {
        TRIMGALORE (
            INPUT_CHECK.out.reads
        ) 
        ch_trimmed_reads = TRIMGALORE.out.reads
    } else {
        ch_trimmed_reads = INPUT_CHECK.out.reads
    }

    //
    // MODULE: Run bwa index
    //
    BWA_INDEX (
        ch_fasta
    )

    //
    // MODULE: Run samtools faidx
    //
    SAMTOOLS_FAIDX (
        ch_fasta
    )

    //
    // MODULE: BWA MEM ALIGNMENT
    BWA_MEM (
        ch_trimmed_reads,
        BWA_INDEX.out.index
    )
    
    SAMTOOLS_SORT_QNAME (
        BWA_MEM.out.sorted_bam
    )


    SAMTOOLS_INDEX_BWA (
        BWA_MEM.out.sorted_bam
    )

    //
    // SUBWORKFLOW - RUN CIRCLE_FINDER PIPELINE
    //
    if (params.circle_identifier == "circle_finder") {
        SAMBLASTER (
            SAMTOOLS_SORT_QNAME.out.bam
        )
        ch_bwa_sorted_bam = BWA_MEM.out.sorted_bam
        ch_bwa_sorted_bai = SAMTOOLS_INDEX_BWA.out.bai

        BEDTOOLS_SPLITBAM2BED (
            SAMBLASTER.out.split_bam
        )
        BEDTOOLS_SORTEDBAM2BED (
            ch_bwa_sorted_bam.join(ch_bwa_sorted_bai)
        )

        ch_b2b_sorted = BEDTOOLS_SORTEDBAM2BED.out.conc_txt
        ch_b2b_split = BEDTOOLS_SPLITBAM2BED.out.split_txt
        CIRCLEFINDER (
            BEDTOOLS_BAM2BED.out.split_conc_txt
        )
    }

    //
    // SUBWORKFLOW: RUN CIRCLE-MAP REALIGN or REPEATS PIPELINE
    //

    if (params.circle_identifier == "circle_map_realign" ||
            params.circle_identifreier == "circle_map_repeats") {
        CIRCLEMAP_READEXTRACTOR (
            SAMTOOLS_SORT_QNAME.out.bam
        )

        SAMTOOLS_SORT_RE (
            CIRCLEMAP_READEXTRACTOR.out.circ_bam
        )

        SAMTOOLS_INDEX_RE (
            SAMTOOLS_SORT_RE.out.bam
        )


        // DEFINE CHANNELS FOR REALIGN AND REPEATS
        ch_bam_sorted = BWA_MEM.out.sorted_bam
        ch_bam_sorted_bai = SAMTOOLS_INDEX_BWA.out.bai
        ch_qname_sorted_bam = SAMTOOLS_SORT_QNAME.out.bam
        ch_re_sorted_bam = SAMTOOLS_SORT_RE.out.bam
        ch_re_sorted_bai = SAMTOOLS_INDEX_RE.out.bai
        ch_fasta_index = SAMTOOLS_FAIDX.out.fai

        // 
        // MODULE: RUN CIRCLE_MAP REPEATS
        //

        if (params.circle_identifier == "circle_map_repeats") {
            ch_sorted_bam = SAMTOOLS_SORT_RE.out.bam
            ch_sorted_bai = SAMTOOLS_INDEX_RE.out.bai
            CIRCLEMAP_REPEATS (
                ch_re_sorted_bam.join(ch_re_sorted_bai)
            )
        }

        //
        // MODULE: Run Circle-Map Realign
        //
        if (params.circle_identifier == "circle_map_realign") {
            CIRCLEMAP_REALIGN (
                ch_re_sorted_bam.join(ch_re_sorted_bai).
                    join(ch_qname_sorted_bam).
                    join(ch_bam_sorted).
                    join(ch_bam_sorted_bai),
                ch_fasta,
                ch_fasta_index
            )
        }
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
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.zip.collect{it[1]}.ifEmpty([]))

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
