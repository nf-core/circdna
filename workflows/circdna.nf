/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS & LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_circdna_pipeline'


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
include { PICARD_ADDORREPLACEREADGROUPS             }   from '../modules/nf-core/picard/addorreplacereadgroups/main'

// PICARD
include { SAMTOOLS_FAIDX                            }   from '../modules/nf-core/samtools/faidx/main'
include { BAM_MARKDUPLICATES_PICARD                 }   from '../subworkflows/nf-core/bam_markduplicates_picard/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FILTER     }   from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FILTERED   }   from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED }   from '../modules/nf-core/samtools/index/main'

// BAM STATS
include { BAM_STATS_SAMTOOLS                        }   from '../subworkflows/nf-core/bam_stats_samtools/main'

// CIRCLE-MAP
include { CIRCLEMAP_READEXTRACTOR                   }   from '../modules/local/circlemap/readextractor/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_RE         }   from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RE       }   from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QNAME_CM   }   from '../modules/nf-core/samtools/sort/main'
include { CIRCLEMAP_REALIGN                         }   from '../modules/local/circlemap/realign/main'
include { CIRCLEMAP_REPEATS                         }   from '../modules/local/circlemap/repeats/main'

// CIRCLE_FINDER
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QNAME_CF   }   from '../modules/nf-core/samtools/sort/main'
include { SAMBLASTER                                }     from '../modules/local/samblaster/main'
include { BEDTOOLS_SORTEDBAM2BED                    }     from '../modules/local/bedtools/sortedbam2bed/main'
include { BEDTOOLS_SPLITBAM2BED                     }     from '../modules/local/bedtools/splitbam2bed/main'
include { CIRCLEFINDER                              }     from '../modules/local/circlefinder/main'

// CIRCexplorer2
include { CIRCEXPLORER2_PARSE       }     from '../modules/nf-core/circexplorer2/parse/main.nf'

// AmpliconArchitect
include { AMPLICONSUITE                                 }     from '../modules/local/ampliconsuite/main'

// Unicycler
include { UNICYCLER           }     from '../modules/nf-core/unicycler/main'
include { SEQTK_SEQ           }     from '../modules/local/seqtk/seq/main'
include { GETCIRCULARREADS    }     from '../modules/local/getcircularreads/main'
include { MINIMAP2_ALIGN      }     from '../modules/nf-core/minimap2/align/main.nf'


// MULTIQC
include { MULTIQC }     from '../modules/nf-core/multiqc/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CIRCDNA {
    take:
    samplesheet
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:
    ch_versions = Channel.empty()
    ch_versions_topic = Channel.topic('versions')
        .map { process, tool, version ->
            "${process}:\n    ${tool}: ${version.toString().trim()}"
        }
    multiqc_report = Channel.empty()

    // FASTA reference channel
    if (params.fasta) {
        ch_fasta = Channel.fromPath(params.fasta)
    } else {
        error 'Fasta reference genome not specified!'
    }

    // Modify fasta channel to include meta data
    ch_fasta_meta = ch_fasta.map { it -> [[id: it[0].baseName], it] }.collect()

    // Circle identifier branches
    def branch_list = params.circle_identifier.split(",")
    run_circexplorer2 = ("circexplorer2" in branch_list)
    run_circle_map_realign = ("circle_map_realign" in branch_list)
    run_circle_map_repeats = ("circle_map_repeats" in branch_list)
    run_circle_finder = ("circle_finder" in branch_list)
    run_ampliconarchitect = ("ampliconarchitect" in branch_list)
    run_unicycler = ("unicycler" in branch_list)

    if (!(run_unicycler | run_circle_map_realign | run_circle_map_repeats | run_circle_finder | run_ampliconarchitect | run_circexplorer2)) {
        error 'circle_identifier param not valid. Please check!'
    }

    // Check if BWA Index is given
    if (params.bwa_index) {
        ch_bwa_index = Channel.fromPath(params.bwa_index, type: 'dir').collect()
        ch_bwa_index = ch_bwa_index.map { index -> ["bwa_index", index] }.collect()
        bwa_index_exists = true
    } else {
        ch_bwa_index = Channel.empty()
        bwa_index_exists = false
    }

    // AmpliconArchitect input validation
    if (run_ampliconarchitect) {
        if (!params.mosek_license_dir) {
            error "Mosek License Directory is missing! Please specify directory containing mosek license using --mosek_license_dir and rename license to 'mosek.lic'."
        }
        mosek_license_dir = file(params.mosek_license_dir)
        if (!params.aa_data_repo) {
            error "AmpliconArchitect Data Repository Missing! Please see https://github.com/jluebeck/AmpliconArchitect for more information and specify its absolute path using --aa_data_repo."
        }
        if (!(params.reference_build in ["hg19", "GRCh38", "GRCh37", "mm10"])) {
            error "Reference Build not given! Please specify --reference_build 'mm10', 'hg19', 'GRCh38', or 'GRCh37'."
        }
        if (!params.cnvkit_cnn) {
            ch_cnvkit_reference = file(params.aa_data_repo + "/" + params.reference_build + "/" + params.reference_build + "_cnvkit_filtered_ref.cnn", checkIfExists: true)
        } else {
            ch_cnvkit_reference = file(params.cnvkit_cnn)
        }
    }

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
        samplesheet
        .map { meta, fastq ->
            meta.id = meta.id.replaceFirst(/_T\\d+$/, '')
            def files = fastq instanceof List ? fastq : [ fastq ]
            def expected = meta.single_end ? 1 : 2
            if (files.size() < expected) {
                error("Unexpected number of FASTQ files for sample ${meta.id}: ${files.size()}")
            }
            [ meta, files.flatten() ]
        }
        .branch {
            meta, files ->
                def expected = meta.single_end ? 1 : 2
                single   : files.size() == expected
                multiple : files.size() > expected
        }
        .set { ch_fastq }

        //
        // MODULE: Concatenate FASTQs from the same samples
        //
        CAT_FASTQ (
            ch_fastq.multiple
        )
        .reads
        .mix(ch_fastq.single)
        .set { ch_cat_fastq }


        //
        // MODULE: Run FastQC
        //
        ch_fastqc_multiqc = Channel.empty()
        if ( ! params.skip_qc ) {
            FASTQC (
                ch_cat_fastq
            )
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
        } else {
            ch_trimmed_reads            = ch_cat_fastq
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
        }


        //
        // MODULE: BWA MEM ALIGNMENT
        //
        if (run_ampliconarchitect | run_circexplorer2 | run_circle_finder |
            run_circle_map_realign | run_circle_map_repeats) {
            BWA_MEM (
                ch_trimmed_reads,
                ch_bwa_index,
                ch_fasta_meta,
                Channel.value(true)
            )
            ch_bam_sorted = BWA_MEM.out.bam
            ch_full_bam_sorted   = BWA_MEM.out.bam
            ch_bwa_sorted   = BWA_MEM.out.bam

            // SAMTOOLS INDEX SORTED BAM
            SAMTOOLS_INDEX_BAM (
                ch_bam_sorted
            )
        }
    } else if (params.input_format == "BAM") {
        // Use BAM Files as input
        samplesheet
        .map { meta, bams ->
            def bam_list = bams instanceof List ? bams : [ bams ]
            if (bam_list.size() != 1) {
                error("Multiple BAMs per sample are not supported: ${meta.id}")
            }
            [ meta, bam_list[0] ]
        }
        .set { ch_bam_input }
        if (!params.bam_sorted){
            SAMTOOLS_SORT_BAM (
                ch_bam_input,
                [[id: null], [], []],
                []
            )
            ch_bam_sorted       = SAMTOOLS_SORT_BAM.out.bam
        } else {
            ch_bam_sorted       = ch_bam_input
            ch_full_bam_sorted  = ch_bam_input
            ch_bwa_sorted       = ch_bam_input
        }
        // SAMTOOLS INDEX SORTED BAM
        SAMTOOLS_INDEX_BAM (
            ch_bam_sorted
        )
        ch_fastqc_multiqc           = Channel.empty()
        ch_trimgalore_multiqc       = Channel.empty()
        ch_trimgalore_multiqc_log   = Channel.empty()
    }




    if (run_ampliconarchitect | run_circexplorer2 | run_circle_finder |
        run_circle_map_realign | run_circle_map_repeats) {

        // Define Index channel and additional bam sorted channels for Circle_finder - not usable with duplicates removed
        ch_bam_sorted_bai       = SAMTOOLS_INDEX_BAM.out.index
        ch_full_bam_sorted      = ch_bam_sorted
        ch_full_bam_sorted_bai  = SAMTOOLS_INDEX_BAM.out.index

        ch_fasta = ch_fasta_meta.map{ meta, index -> [index] }.collect()

        // Stub run is not yet implemented into BAM_STATS_SAMTOOLS subworkflow -> Will be skipped when stub is active
        ch_bam_with_bai = ch_bam_sorted.join(ch_bam_sorted_bai).map { meta, bam, bai -> [meta, bam, bai] }

        BAM_STATS_SAMTOOLS(
            ch_bam_with_bai,
            ch_fasta_meta.map { meta, fasta -> [meta, fasta, []] }
        )
        ch_samtools_stats    = BAM_STATS_SAMTOOLS.out.stats
        ch_samtools_flagstat = BAM_STATS_SAMTOOLS.out.flagstat
        ch_samtools_idxstats = BAM_STATS_SAMTOOLS.out.idxstats

        // PICARD MARK_DUPLICATES
        if (!params.skip_markduplicates) {
            // Index Fasta File for Markduplicates
            SAMTOOLS_FAIDX (
                ch_fasta_meta.map { meta, fasta -> [meta, fasta, []] },
                false
            )

            // Combine fasta and fai into [meta, fasta, fai] tuple required by picard modules
            ch_fasta_fai = SAMTOOLS_FAIDX.out.fa.join(SAMTOOLS_FAIDX.out.fai)

            PICARD_ADDORREPLACEREADGROUPS (
                ch_bam_sorted,
                ch_fasta_fai
            )

            ch_bam_md_input = PICARD_ADDORREPLACEREADGROUPS.out.bam

            // MARK DUPLICATES IN BAM FILE
            BAM_MARKDUPLICATES_PICARD (
                ch_bam_md_input,
                ch_fasta_fai
            )

            // FILTER DUPLICATES IN BAM FILES USING SAMTOOLS VIEW
            if (!params.keep_duplicates) {
                SAMTOOLS_VIEW_FILTER (
                    BAM_MARKDUPLICATES_PICARD.out.bam.join(BAM_MARKDUPLICATES_PICARD.out.index),
                    ch_fasta_fai,
                    [[], []],
                    [[], []],
                    []
                )

                // SORT FILTERED BAM FILE
                SAMTOOLS_SORT_FILTERED (
                    SAMTOOLS_VIEW_FILTER.out.bam,
                    [[id: null], [], []],
                    []
                )

                // INDEX FILTERED BAM FILE
                SAMTOOLS_INDEX_FILTERED (
                    SAMTOOLS_SORT_FILTERED.out.bam
                )

                ch_bam_sorted = SAMTOOLS_SORT_FILTERED.out.bam
                ch_bam_sorted_bai = SAMTOOLS_INDEX_FILTERED.out.index
            }
            else {
                ch_bam_sorted               = BAM_MARKDUPLICATES_PICARD.out.bam
                ch_bam_sorted_bai           = BAM_MARKDUPLICATES_PICARD.out.index
                ch_markduplicates_stats     = BAM_MARKDUPLICATES_PICARD.out.stats
                ch_markduplicates_flagstat  = BAM_MARKDUPLICATES_PICARD.out.flagstat
                ch_markduplicates_idxstats  = BAM_MARKDUPLICATES_PICARD.out.idxstats
                ch_markduplicates_multiqc   = BAM_MARKDUPLICATES_PICARD.out.metrics
            }
        } else {
                ch_markduplicates_stats         = Channel.empty()
                ch_markduplicates_flagstat      = Channel.empty()
                ch_markduplicates_idxstats      = Channel.empty()
                ch_markduplicates_multiqc       = Channel.empty()
        }
    }

    if (run_ampliconarchitect) {
        AMPLICONSUITE (
            ch_bam_sorted,
            file(params.mosek_license_dir),
            file(params.aa_data_repo)
        )
    }

    //
    // SUBWORKFLOW - RUN CIRCLE_FINDER PIPELINE
    //
    if (run_circle_finder) {
        SAMTOOLS_SORT_QNAME_CF (
            ch_full_bam_sorted,
            [[id: null], [], []],
            []
        )

        SAMBLASTER (
            SAMTOOLS_SORT_QNAME_CF.out.bam
        )

        BEDTOOLS_SPLITBAM2BED (
            SAMBLASTER.out.split_bam
        )

        BEDTOOLS_SORTEDBAM2BED (
            ch_full_bam_sorted.join(ch_full_bam_sorted_bai)
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
        SAMTOOLS_SORT_QNAME_CM (
            ch_bam_sorted,
            [[id: null], [], []],
            []
        )

        CIRCLEMAP_READEXTRACTOR (
            SAMTOOLS_SORT_QNAME_CM.out.bam
        )

        SAMTOOLS_SORT_RE (
            CIRCLEMAP_READEXTRACTOR.out.bam,
            [[id: null], [], []],
            []
        )

        SAMTOOLS_INDEX_RE (
            SAMTOOLS_SORT_RE.out.bam
        )

        // DEFINE CHANNELS FOR REALIGN AND REPEATS
        ch_qname_sorted_bam = SAMTOOLS_SORT_QNAME_CM.out.bam
        ch_re_sorted_bam = SAMTOOLS_SORT_RE.out.bam
        ch_re_sorted_bai = SAMTOOLS_INDEX_RE.out.index

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
                ch_re_sorted_bam
                    .join(ch_re_sorted_bai)
                    .join(ch_qname_sorted_bam)
                    .join(ch_bam_sorted)
                    .join(ch_bam_sorted_bai),
                ch_fasta
            )
        }
    }


    if (run_circexplorer2) {
        ch_bam_sorted
            .join(ch_bam_sorted_bai)
            .map { meta, bam, bai -> [meta, bam] }
            .set { ch_circexplorer2_input }

        CIRCEXPLORER2_PARSE (
            ch_circexplorer2_input
        )
    }

    if (run_unicycler && params.input_format == "FASTQ") {

        UNICYCLER (
            ch_trimmed_reads.map { meta, reads -> [meta, reads, []] }
        )

        SEQTK_SEQ (
            UNICYCLER.out.scaffolds
        )

        GETCIRCULARREADS (
            SEQTK_SEQ.out.fastq
        )

        GETCIRCULARREADS.out.fastq
            .map { meta, fastq -> [ meta + [single_end: true], fastq ] }
            .set { ch_circular_fastq }

        MINIMAP2_ALIGN (
            ch_circular_fastq,
            ch_fasta_meta,
            false,
            false,
            false,
            false
        )
    }

    //
    // MODULE: Pipeline reporting
    //
//    CUSTOM_DUMPSOFTWAREVERSIONS (
//        ch_versions.unique().collectFile(name: 'collated_versions.yml')
//    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'circdna_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
        def ch_multiqc_custom_methods_description = multiqc_methods_description
            ? file(multiqc_methods_description, checkIfExists: true)
            : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

        def ch_multiqc_files = ch_fastqc_multiqc.collect{it[1]}.ifEmpty([])
            .mix(ch_trimgalore_multiqc.collect{it[1]}.ifEmpty([]))
            .mix(ch_trimgalore_multiqc_log.collect{it[1]}.ifEmpty([]))
            .mix(ch_samtools_stats.collect{it[1]}.ifEmpty([]))
            .mix(ch_samtools_flagstat.collect{it[1]}.ifEmpty([]))
            .mix(ch_samtools_idxstats.collect{it[1]}.ifEmpty([]))
            .mix(ch_markduplicates_stats.collect{it[1]}.ifEmpty([]))
            .mix(ch_markduplicates_flagstat.collect{it[1]}.ifEmpty([]))
            .mix(ch_markduplicates_idxstats.collect{it[1]}.ifEmpty([]))
            .mix(ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([]))
            .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
            .mix(ch_collated_versions)
            .mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

        MULTIQC(
            ch_multiqc_files.flatten().collect().map { files ->
                [
                    [id: 'circdna'],
                    files,
                    multiqc_config
                        ? file(multiqc_config, checkIfExists: true)
                        : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                    multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
        )
        multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
    }

    emit:
    multiqc_report = multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
