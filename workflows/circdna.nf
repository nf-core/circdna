/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-schema'
include { INPUT_CHECK           } from '../subworkflows/local/input_check/main'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_circdna_pipeline/main'
include { CAT_FASTQ     }     from '../modules/nf-core/cat/fastq/main'
include { FASTQC        }     from '../modules/nf-core/fastqc/main'
include { TRIMGALORE    }    from '../modules/nf-core/trimgalore/main'
include { BWA_INDEX     }   from '../modules/nf-core/bwa/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_BAM        }   from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BAM      }   from '../modules/nf-core/samtools/index/main'
include { CIRCEXPLORER2_PARSE       }     from '../modules/nf-core/circexplorer2/parse/main'
include { MULTIQC }     from '../modules/nf-core/multiqc/main'

include { BAM_PREPROCESSING              } from '../subworkflows/local/bam_preprocessing/main'
include { CIRCLE_MAP_PIPELINE            } from '../subworkflows/local/circle_map_pipeline/main'
include { CIRCLE_FINDER_PIPELINE         } from '../subworkflows/local/circle_finder_pipeline/main'
include { AMPLICONARCHITECT_PIPELINE     } from '../subworkflows/local/ampliconarchitect_pipeline/main'
include { UNICYCLER_PIPELINE             } from '../subworkflows/local/unicycler_pipeline/main'
include { LONG_READ_PREPROCESSING        } from '../subworkflows/local/long_read_preprocessing/main'
include { LONG_READ_REFERENCE            } from '../subworkflows/local/long_read_reference/main'
include { CRESIL_PIPELINE               } from '../subworkflows/local/cresil_pipeline/main'
include { FLED_PIPELINE                 } from '../subworkflows/local/fled_pipeline/main'
include { FLYE_PIPELINE                 } from '../subworkflows/local/flye_pipeline/main'
include { REMAP_ASSEMBLED_CIRCLES       } from '../subworkflows/local/remap_assembled_circles/main'
include { LONG_READ_FILTERING as LONG_READ_FILTERING_CRESIL } from '../subworkflows/local/long_read_filtering/main'
include { LONG_READ_FILTERING as LONG_READ_FILTERING_FLED   } from '../subworkflows/local/long_read_filtering/main'
include { LONG_READ_FILTERING as LONG_READ_FILTERING_FLYE   } from '../subworkflows/local/long_read_filtering/main'
include { LONG_READ_FILTERING as LONG_READ_FILTERING_CIRCLESEEKER } from '../subworkflows/local/long_read_filtering/main'
include { LONG_READ_FILTERING as LONG_READ_FILTERING_ECCFINDER } from '../subworkflows/local/long_read_filtering/main'
include { LONG_READ_FILTERING as LONG_READ_FILTERING_ECCFINDER_ASM } from '../subworkflows/local/long_read_filtering/main'
include { CIRCLESEEKER_PIPELINE         } from '../subworkflows/local/circleseeker_pipeline/main'
include { CRESIL_TO_BED                   } from '../modules/local/cresil_to_bed/main'
include { COLLECT_LONG_READ_BED          } from '../modules/local/collect_long_read_bed/main'
include { CIDERSEQ_PIPELINE             } from '../subworkflows/local/ciderseq_pipeline/main'
include { ORGANELLE_TAG                 } from '../subworkflows/local/organelle_tag/main'
include { REFERENCE_MODE                } from '../subworkflows/local/reference_mode/main'
include { ECCDNA_MODE                   } from '../subworkflows/local/eccdna_mode/main'
include { ECCSPLORER_PIPELINE            } from '../subworkflows/local/eccsplorer_pipeline/main'
include { ECCSPLORER_SLIM_PIPELINE       } from '../subworkflows/local/eccsplorer_slim_pipeline/main'
include { ECCSPLORER_CLU_SLIM            } from '../subworkflows/local/eccsplorer_clu_slim/main'
include { ECCSPLORER_ALL_SLIM            } from '../subworkflows/local/eccsplorer_all_slim/main'
include { ECC_FINDER_SLIM_PIPELINE       } from '../subworkflows/local/ecc_finder_slim_pipeline/main'
include { ECC_FINDER_ONT_SLIM            } from '../subworkflows/local/ecc_finder_ont_slim/main'

workflow CIRCDNA {
    // Mode validation (§1.1: mode is NOT a master switch anymore)
    def valid_modes = ['reference', 'eccdna']
    if (!valid_modes.contains(params.mode)) {
        exit 1, "Invalid mode: ${params.mode}. Valid modes: ${valid_modes.join(', ')}"
    }

    if (params.containsKey('fasta') && params.fasta) {
        ch_fasta = channel.fromPath(params.fasta).collect()
    } else {
        def genome_fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
        if (genome_fasta) {
            params.fasta = genome_fasta
            ch_fasta = channel.fromPath(genome_fasta).collect()
        } else {
            exit 1, 'Fasta reference genome not specified!'
        }
    }
    if (!(params.input_format == "FASTQ" | params.input_format == "BAM")) {
    exit 1, 'Please specifiy --input_format "FASTQ" or "BAM" in capital letters, depending on the input file format.'
    }
    // ch_fasta 经 .collect() 已是 value channel，无需 .first()（否则触发
    // "first is useless on value channel" 警告）
    ch_fasta_meta = ch_fasta.map{ fasta -> [[id: fasta.baseName], fasta] }

    // ---------------- short-read identifier parsing (legacy / slim / blackbox) ----------------
    def use_legacy_mode = false
    def run_circexplorer2 = false
    def run_circle_map_realign = false
    def run_circle_map_repeats = false
    def run_circle_finder = false
    def run_ampliconarchitect = false
    def run_unicycler = false
    def run_eccsplorer_blackbox = false
    def run_ecc_finder_map_sr_blackbox = false
    def run_ecc_finder_asm_sr_blackbox = false
    def run_ecc_finder_map_ont_blackbox = false
    def run_ecc_finder_asm_ont_blackbox = false

    if (params.circle_identifier) {
        def branch = params.circle_identifier.split(",")
        run_circexplorer2 = ("circexplorer2" in branch)
        run_circle_map_realign = ("circle_map_realign" in branch)
        run_circle_map_repeats = ("circle_map_repeats" in branch)
        run_circle_finder = ("circle_finder" in branch)
        run_ampliconarchitect = ("ampliconarchitect" in branch)
        run_unicycler = ("unicycler" in branch)
        run_eccsplorer_blackbox = ("eccsplorer" in branch)
        run_ecc_finder_map_sr_blackbox = ("ecc_finder_map_sr" in branch)
        run_ecc_finder_asm_sr_blackbox = ("ecc_finder_asm_sr" in branch)
        run_ecc_finder_map_ont_blackbox = ("ecc_finder_map_ont" in branch)
        run_ecc_finder_asm_ont_blackbox = ("ecc_finder_asm_ont" in branch)
        def use_blackbox_any = (run_eccsplorer_blackbox | run_ecc_finder_map_sr_blackbox | run_ecc_finder_asm_sr_blackbox | run_ecc_finder_map_ont_blackbox | run_ecc_finder_asm_ont_blackbox)
        def use_legacy_any = (run_unicycler | run_circle_map_realign | run_circle_map_repeats | run_circle_finder | run_ampliconarchitect | run_circexplorer2)
        // slim 标识符（原子化链，不走 legacy 路径）
        def use_slim_any = branch.any { ident -> ident in ['eccsplorer_map_slim','eccsplorer_clu_slim','eccsplorer_all_slim','ecc_finder_map_sr_slim','ecc_finder_asm_sr_slim','ecc_finder_map_ont_slim','ecc_finder_asm_ont_slim'] }
        use_legacy_mode = use_legacy_any
        if (!use_legacy_any && !use_slim_any && !use_blackbox_any) {
            exit 1, 'circle_identifier param not valid. Please check!'
        }
        if (run_unicycler && !params.input_format == "FASTQ") {
            exit 1, 'Unicycler needs FastQ input. Please specify input_format == "FASTQ", if possible, or don`t run unicycler.'
        }
    }

    if (!params.input) { exit 1, 'Input samplesheet not specified!' }
    def bwa_index_exists = false
    def ch_bwa_index = channel.empty()
    if (params.bwa_index) {
    ch_bwa_index = channel.fromPath(params.bwa_index, type: 'dir').map{ index -> [[id: 'bwa_index'], index] }
    bwa_index_exists = true
    } else {
    ch_bwa_index = channel.empty()
    bwa_index_exists = false
    }
    def mosek_license_dir = null
    def ch_cnvkit_reference = channel.empty()
    if (run_ampliconarchitect) {
    mosek_license_dir = params.mosek_license_dir
    if (!params.mosek_license_dir) {
        exit 1, "Mosek License Directory is missing! Please specifiy directory containing mosek license using --mosek_license_dir and rename license to 'mosek.lic'."
    } else {
        mosek_license_dir = file(params.mosek_license_dir)
    }
    if (!params.aa_data_repo) { exit 1, "AmpliconArchitect Data Repository Missing! Please see https://github.com/jluebeck/AmpliconArchitect for more information and specify its absolute path using --aa_data_repo." }
    if (params.reference_build != "hg19" & params.reference_build != "GRCh38" & params.reference_build != "GRCh37" & params.reference_build != "mm10"){
        exit 1, "Reference Build not given! Please specify --reference_build 'mm10', 'hg19', 'GRCh38', or 'GRCh37'."
    }
    if (!params.cnvkit_cnn) {
        ch_cnvkit_reference = file(params.aa_data_repo + "/" + params.reference_build + "/" + params.reference_build + "_cnvkit_filtered_ref.cnn", checkIfExists: true)
    } else {
        ch_cnvkit_reference = file(params.cnvkit_cnn)
    }
    }

    // ---------------- long-read identifier parsing (engine switches, §3) ----------------
    def lr_branch = (params.long_read_identifier ?: '').split(',').collect { ident -> ident.trim() }.findAll { ident -> ident }
    def run_cresil       = ("cresil" in lr_branch)
    def run_fled         = ("fled" in lr_branch)
    def run_flye         = ("flye" in lr_branch)
    def run_eccfinder    = ("eccfinder" in lr_branch)
    def run_circleseeker = ("circleseeker" in lr_branch)
    def run_ciderseq     = ("ciderseq" in lr_branch)

    if (run_ciderseq) {
        if (!params.ciderseq_config) {
            exit 1, 'ciderseq identifier requires --ciderseq_config (CIDER-Seq2 config JSON)'
        }
        if (!params.ciderseq_blastdb) {
            exit 1, 'ciderseq identifier requires --ciderseq_blastdb (directory of the blastn database)'
        }
        if (!params.ciderseq_align_targets) {
            exit 1, 'ciderseq identifier requires --ciderseq_align_targets (directory of align target fastas)'
        }
        if (!params.ciderseq_protein_db) {
            exit 1, 'ciderseq identifier requires --ciderseq_protein_db (directory of the tblastn protein database)'
        }
    }

    // CONFIG FILES
    def ch_multiqc_config          = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_custom_config   = params.multiqc_config ? channel.fromPath( params.multiqc_config, checkIfExists: true ) : channel.empty()
    def _ch_multiqc_logo            = params.multiqc_logo   ? channel.fromPath( params.multiqc_logo, checkIfExists: true ) : channel.empty()
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    // RUN MAIN WORKFLOW
    def multiqc_report = []

    ch_versions = channel.empty()

    // Define Empty Channels for MultiQC
    ch_samtools_stats           = channel.empty()
    ch_samtools_flagstat        = channel.empty()
    ch_samtools_idxstats        = channel.empty()
    ch_markduplicates_stats     = channel.empty()
    ch_markduplicates_flagstat  = channel.empty()
    ch_markduplicates_idxstats  = channel.empty()
    ch_markduplicates_multiqc   = channel.empty()
    ch_fastqc_multiqc           = channel.empty()
    ch_trimgalore_multiqc       = channel.empty()
    ch_trimgalore_multiqc_log   = channel.empty()
    ch_mosdepth_multiqc         = channel.empty()

    // Unified long-read candidate BED (eccdna.smk 消费)
    ch_long_read_bed = channel.empty()

    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_short_reads = INPUT_CHECK.out.reads_short   // [ meta, [reads] ] or [ meta, bam ]
    ch_long_reads  = INPUT_CHECK.out.reads_long    // [ meta, fastq, input_bam, entrypoint ]

    // =================================================================================
    // SHORT-READ PATH（illumina 行；按 meta.datatype 分流 gdna / eccdna，§2）
    // =================================================================================
    if ( params.input_format == "FASTQ" ) {
        ch_short_reads
        .map {
            meta, fastq ->
                if (!meta.containsKey('lane')) {
                    def parts = meta.id.split('_')
                    if (parts.size() > 1 && parts[-1] ==~ /T\d+/) {
                        meta.id = parts[0..-2].join('_')
                    }
                }
                [ meta.id, meta, fastq ] }
        .groupTuple(by: [0])
        .map { id, meta_list, fastq_list ->
            def merged_meta = meta_list[0].clone()
            merged_meta.id = id
            [ merged_meta, fastq_list ] }
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1 && meta.single_end
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1 || !meta.single_end
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }

        CAT_FASTQ (
            ch_fastq.multiple
        )
        .reads
        .mix(ch_fastq.single)
        .set { ch_cat_fastq }

        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions_cat)

        ch_fastqc_multiqc = channel.empty()
        if ( ! params.skip_qc ) {
            FASTQC (
                ch_cat_fastq
            )
            ch_versions         = ch_versions.mix(FASTQC.out.versions_fastqc)
            ch_fastqc_multiqc   = FASTQC.out.zip
        }

        if ( ! params.skip_trimming ) {
            TRIMGALORE (
                ch_cat_fastq
            )
            ch_trimmed_reads            = TRIMGALORE.out.reads
            ch_trimgalore_multiqc       = TRIMGALORE.out.zip
            ch_trimgalore_multiqc_log   = TRIMGALORE.out.log
            ch_versions                 = ch_versions.mix(TRIMGALORE.out.versions_trimgalore)
        } else {
            ch_trimmed_reads            = ch_cat_fastq
            ch_trimgalore_multiqc       = channel.empty()
            ch_trimgalore_multiqc_log   = channel.empty()
        }

        // BWA index：由短读行触发（Nextflow 26 中 if(ch) 恒 true，不能用空 channel 门控；
        // 纯长读 run 时 trigger 为空 → BWA_INDEX 无任务）
        if (!bwa_index_exists) {
            def ch_bwa_trigger = ch_short_reads.map { meta, _r -> meta }.first()
            BWA_INDEX (
                ch_bwa_trigger
                    .combine(ch_fasta_meta)
                    .map { _t, m_fa, fa -> [ m_fa, fa ] }
            )
            // .first() 必须：BWA_INDEX.out.index 是队列通道（1 元素），BWA_MEM 多队列输入
            // 会 one-to-one 配对只匹配 1 个样本；转 value 后广播给所有 trimmed_reads
            ch_bwa_index = BWA_INDEX.out.index.map{ _meta, index -> [[id: 'bwa_index'], index] }.first()
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions_bwa)
            bwa_index_exists = true
        }

        def ch_sr_gdna = ch_trimmed_reads.filter { meta, _reads -> meta.datatype == 'gdna' }
        def ch_sr_ecc  = ch_trimmed_reads.filter { meta, _reads -> meta.datatype == 'eccdna' }

        if (use_legacy_mode) {
            BAM_PREPROCESSING (
                ch_trimmed_reads,
                ch_bwa_index,
                ch_fasta_meta,
                true
            )
            ch_bam_sorted               = BAM_PREPROCESSING.out.bam_sorted
            ch_bam_sorted_bai           = BAM_PREPROCESSING.out.bam_sorted_bai
            ch_full_bam_sorted          = BAM_PREPROCESSING.out.full_bam_sorted
            ch_full_bam_sorted_bai      = BAM_PREPROCESSING.out.full_bam_sorted_bai
            ch_samtools_stats           = BAM_PREPROCESSING.out.samtools_stats
            ch_samtools_flagstat        = BAM_PREPROCESSING.out.samtools_flagstat
            ch_samtools_idxstats        = BAM_PREPROCESSING.out.samtools_idxstats
            ch_markduplicates_stats     = BAM_PREPROCESSING.out.markduplicates_stats
            ch_markduplicates_flagstat  = BAM_PREPROCESSING.out.markduplicates_flagstat
            ch_markduplicates_idxstats  = BAM_PREPROCESSING.out.markduplicates_idxstats
            ch_markduplicates_multiqc   = BAM_PREPROCESSING.out.markduplicates_multiqc
            ch_versions = ch_versions.mix(BAM_PREPROCESSING.out.versions)
        }

        // 短读 gdna → REFERENCE_MODE（BWA + mosdepth）
        if ( ch_sr_gdna && !use_legacy_mode ) {
            REFERENCE_MODE (
                ch_sr_gdna,
                ch_bwa_index,
                ch_fasta_meta,
                params.repeat_gff ? file(params.repeat_gff) : channel.empty()
            )
            ch_versions = ch_versions.mix(REFERENCE_MODE.out.versions)
            ch_mosdepth_multiqc = ch_mosdepth_multiqc.mix(REFERENCE_MODE.out.mosdepth_summary.map { _meta, summary -> summary })
        }

        // 短读 eccdna → ECCDNA_MODE（默认 Circle-Map）+ slim / blackbox / legacy identifier
        if ( ch_sr_ecc && !use_legacy_mode ) {
            ECCDNA_MODE (
                ch_sr_ecc,
                ch_bwa_index,
                ch_fasta_meta,
                true
            )
            ch_versions = ch_versions.mix(ECCDNA_MODE.out.versions)
            ch_mosdepth_multiqc = ch_mosdepth_multiqc.mix(ECCDNA_MODE.out.mosdepth_summary.map { _meta, summary -> summary })

            // ================================================================
            // SLIM 原子化链（仅作用于 sr_ecc，§3）
            // ================================================================
            if (params.circle_identifier) {
                def identifiers = params.circle_identifier.split(",")
                def use_eccsplorer_slim = identifiers.contains('eccsplorer_map_slim')
                def use_ecc_finder_map_sr_slim = identifiers.contains('ecc_finder_map_sr_slim')
                def use_ecc_finder_asm_sr_slim = identifiers.contains('ecc_finder_asm_sr_slim')
                def use_eccsplorer_clu_slim = identifiers.contains('eccsplorer_clu_slim')
                def use_ecc_finder_map_ont_slim = identifiers.contains('ecc_finder_map_ont_slim')
                def use_ecc_finder_asm_ont_slim = identifiers.contains('ecc_finder_asm_ont_slim')
                def use_eccsplorer_all_slim = identifiers.contains('eccsplorer_all_slim')
                def use_slim_any = (use_eccsplorer_slim | use_ecc_finder_map_sr_slim | use_ecc_finder_asm_sr_slim |
                    use_eccsplorer_clu_slim | use_ecc_finder_map_ont_slim | use_ecc_finder_asm_ont_slim | use_eccsplorer_all_slim)

                if (use_slim_any) {
                    // datatype 分流（§1.1：按行 datatype）
                    def ch_eccdna_reads = ch_sr_ecc
                    def ch_gdna_reads = ch_sr_gdna

                    if (use_eccsplorer_all_slim) {
                        def taxon_val = params.eccsplorer_taxon ?: 'vir'
                        ECCSPLORER_ALL_SLIM ( ch_eccdna_reads, ch_fasta_meta, ch_gdna_reads, taxon_val )
                        ch_versions = ch_versions.mix(ECCSPLORER_ALL_SLIM.out.versions)
                    }
                    if (use_eccsplorer_slim) {
                        ECCSPLORER_SLIM_PIPELINE ( ch_eccdna_reads, ch_fasta_meta, ch_gdna_reads )
                        ch_versions = ch_versions.mix(ECCSPLORER_SLIM_PIPELINE.out.versions)
                    }
                    if (use_eccsplorer_clu_slim) {
                        // clu join 仅 illumina×circleseq×eccdna ⋈ illumina×wgs×gdna 且 pair 相同（§2）
                        def ch_pair_treatment = ch_eccdna_reads.map { meta, reads -> [ meta.pair, meta, reads ] }
                        def ch_pair_control  = ch_gdna_reads.map { meta, reads -> [ meta.pair, meta, reads ] }
                        def ch_clu_paired = ch_pair_treatment.join(ch_pair_control, remainder: false)
                            .map { _pair, e_meta, e_reads, _c_meta, _c_reads -> [ e_meta, e_reads[0], e_reads[1] ] }
                        def ch_clu_control = ch_pair_treatment.join(ch_pair_control, remainder: false)
                            .map { _pair, e_meta, _e_reads, _c_meta, c_reads -> [ e_meta, c_reads[0], c_reads[1] ] }
                        def taxon_val = params.eccsplorer_taxon ?: 'vir'
                        ECCSPLORER_CLU_SLIM ( ch_clu_paired, ch_clu_control, taxon_val )
                        ch_versions = ch_versions.mix(ECCSPLORER_CLU_SLIM.out.versions)
                    }
                    if (use_ecc_finder_map_sr_slim || use_ecc_finder_asm_sr_slim) {
                        def use_legacy_any = (run_unicycler | run_circle_map_realign | run_circle_map_repeats | run_circle_finder | run_ampliconarchitect | run_circexplorer2)
                        if (!use_legacy_any) {
                            // slim-only：生成 sorted BAM（复用 BAM_PREPROCESSING）
                            BAM_PREPROCESSING ( ch_eccdna_reads, ch_bwa_index, ch_fasta_meta, true )
                            ch_versions = ch_versions.mix(BAM_PREPROCESSING.out.versions)
                        }
                        def ch_slim_bam = use_legacy_any ? ch_bam_sorted : BAM_PREPROCESSING.out.bam_sorted
                        def ch_slim_bai = use_legacy_any ? ch_bam_sorted_bai : BAM_PREPROCESSING.out.bam_sorted_bai
                        ECC_FINDER_SLIM_PIPELINE (
                            ch_slim_bam, ch_slim_bai, ch_fasta_meta, ch_eccdna_reads,
                            use_ecc_finder_map_sr_slim, use_ecc_finder_asm_sr_slim
                        )
                        ch_versions = ch_versions.mix(ECC_FINDER_SLIM_PIPELINE.out.versions)
                    }
                    // 短读 ONT slim identifier 无意义（sr 仅 illumina），此处保留调用但走长读链由 reads_long 处理
                    if (use_ecc_finder_map_ont_slim || use_ecc_finder_asm_ont_slim) {
                        exit 1, "ecc_finder_map_ont_slim / ecc_finder_asm_ont_slim 仅用于长读（platform=pacbio|ont）行；短读表请改用 map_sr/asm_sr slim。"
                    }
                }

                // ================================================================
                // 黑盒链（circle_identifier 含 eccsplorer / ecc_finder_map_sr / ecc_finder_asm_sr）
                // ================================================================
                def use_blackbox_sr_any = (run_eccsplorer_blackbox | run_ecc_finder_map_sr_blackbox | run_ecc_finder_asm_sr_blackbox)
                if (use_blackbox_sr_any) {
                    ECCSPLORER_PIPELINE (
                        ch_sr_ecc,
                        ch_bwa_index,
                        ch_fasta_meta,
                        run_eccsplorer_blackbox,
                        false,
                        run_ecc_finder_map_sr_blackbox,
                        run_ecc_finder_asm_sr_blackbox,
                        false,
                        params.input_format,
                        ch_sr_gdna
                    )
                    ch_versions = ch_versions.mix(ECCSPLORER_PIPELINE.out.versions)
                }
            }
        }
    } else if ( params.input_format == "BAM" ) {
        if (ch_short_reads){
            SAMTOOLS_SORT_BAM (
                ch_short_reads,
                channel.value([]),
                'bai'
            )
            ch_versions         = ch_versions.mix(SAMTOOLS_SORT_BAM.out.versions_samtools)
            ch_bam_sorted       = SAMTOOLS_SORT_BAM.out.bam
        } else {
            ch_bam_sorted       = channel.empty()
        }
        ch_fastqc_multiqc           = channel.empty()
        ch_trimgalore_multiqc       = channel.empty()
        ch_trimgalore_multiqc_log   = channel.empty()

        // BAM 输入 + legacy identifier：BAM_PREPROCESSING（仅索引，不重比对）
        if (use_legacy_mode && (run_ampliconarchitect | run_circexplorer2 | run_circle_finder |
            run_circle_map_realign | run_circle_map_repeats)) {
            BAM_PREPROCESSING (
                ch_bam_sorted,
                ch_bwa_index,
                ch_fasta_meta,
                false
            )
            ch_bam_sorted               = BAM_PREPROCESSING.out.bam_sorted
            ch_bam_sorted_bai           = BAM_PREPROCESSING.out.bam_sorted_bai
            ch_full_bam_sorted          = BAM_PREPROCESSING.out.full_bam_sorted
            ch_full_bam_sorted_bai      = BAM_PREPROCESSING.out.full_bam_sorted_bai
            ch_samtools_stats           = BAM_PREPROCESSING.out.samtools_stats
            ch_samtools_flagstat        = BAM_PREPROCESSING.out.samtools_flagstat
            ch_samtools_idxstats        = BAM_PREPROCESSING.out.samtools_idxstats
            ch_markduplicates_stats     = BAM_PREPROCESSING.out.markduplicates_stats
            ch_markduplicates_flagstat  = BAM_PREPROCESSING.out.markduplicates_flagstat
            ch_markduplicates_idxstats  = BAM_PREPROCESSING.out.markduplicates_idxstats
            ch_markduplicates_multiqc   = BAM_PREPROCESSING.out.markduplicates_multiqc
            ch_versions = ch_versions.mix(BAM_PREPROCESSING.out.versions)
        }
    }

    // 短读 legacy identifier（circle_finder / circle_map / circexplorer2 / unicycler / ampliconarchitect）
    if (use_legacy_mode && (run_ampliconarchitect | run_circexplorer2 | run_circle_finder |
        run_circle_map_realign | run_circle_map_repeats)) {
        if (run_ampliconarchitect) {
            AMPLICONARCHITECT_PIPELINE (
                ch_bam_sorted,
                ch_bam_sorted_bai,
                BAM_PREPROCESSING.out.fasta_fai,
                ch_cnvkit_reference,
                file(params.mosek_license_dir),
                file(params.aa_data_repo)
            )
            ch_versions = ch_versions.mix(AMPLICONARCHITECT_PIPELINE.out.versions)
        }

        if (run_circle_finder) {
            CIRCLE_FINDER_PIPELINE (
                ch_bam_sorted,
                ch_bam_sorted_bai,
                ch_full_bam_sorted,
                ch_full_bam_sorted_bai,
                BAM_PREPROCESSING.out.fasta_fai
            )
            ch_versions = ch_versions.mix(CIRCLE_FINDER_PIPELINE.out)
        }

        if (run_circle_map_realign || run_circle_map_repeats) {
            CIRCLE_MAP_PIPELINE (
                ch_bam_sorted,
                ch_bam_sorted_bai,
                BAM_PREPROCESSING.out.fasta_fai,
                run_circle_map_realign,
                run_circle_map_repeats
            )
            ch_versions = ch_versions.mix(CIRCLE_MAP_PIPELINE.out)
        }

        if (run_circexplorer2) {
            CIRCEXPLORER2_PARSE (
                ch_bam_sorted
            )
            ch_versions = ch_versions.mix(CIRCEXPLORER2_PARSE.out.versions_circexplorer2)
        }

        if (run_unicycler && params.input_format == "FASTQ") {
            UNICYCLER_PIPELINE (
                ch_trimmed_reads,
                ch_fasta_meta
            )
            ch_versions = ch_versions.mix(UNICYCLER_PIPELINE.out)
        }
    }

    // =================================================================================
    // LONG-READ PATH（pacbio / ont 行；按 assay × datatype × concatemer × read_type 分流，§2/§3）
    // =================================================================================
    if ( ch_long_reads ) {
        def ch_lr_gdna  = ch_long_reads.filter { meta, _f, _b, _ep -> meta.datatype == 'gdna' }
        def ch_lr_det   = ch_long_reads.filter { meta, _f, _b, _ep -> meta.datatype == 'eccdna' }

        LONG_READ_PREPROCESSING ( ch_long_reads )
            .preprocessed_fastq
            .set { ch_preprocessed }
        ch_versions = ch_versions.mix(LONG_READ_PREPROCESSING.out.versions)

        // 长读 WGS 背景（gdna）→ LONG_READ_REFERENCE（minimap2 + mosdepth，无检测引擎，§4.4）
        if ( ch_lr_gdna ) {
            LONG_READ_REFERENCE (
                ch_preprocessed.filter { meta, _f -> meta.datatype == 'gdna' },
                ch_fasta_meta
            )
            ch_versions = ch_versions.mix(LONG_READ_REFERENCE.out.versions)
            ch_mosdepth_multiqc = ch_mosdepth_multiqc.mix(LONG_READ_REFERENCE.out.mosdepth_summary.map { _meta, summary -> summary })
        }

        if ( ch_lr_det ) {
            def ch_lr_rca_enriched = ch_preprocessed.filter { meta, _f -> meta.assay in ['rca', 'enriched'] }
            def ch_lr_cider        = ch_preprocessed.filter { meta, _f -> meta.assay == 'ciderseq' }

            // ---- CRESIL（rca / enriched；§3 允许） ----
            if ( run_cresil && ch_lr_rca_enriched ) {
                CRESIL_PIPELINE (
                    ch_lr_rca_enriched,
                    ch_fasta
                )
                .eccdna_candidates
                .map { meta, file -> [ meta, file ] }
                .set { ch_cresil_raw }

                // CReSIL 原始表（id / merge_region / numreads ...）必须先用
                // convert_cresil_to_bed.py 转成统一 BED 契约再进 LONG_READ_FILTERING。
                // 旧实现直接把 10 列原始表喂进去：filter_by_read_support.py 既识别不到
                // 支持度列（回退到第 5 列 ctc，int("False") 抛错后整行保留）、也无法产出
                // BED，导致长读终集几乎为空。
                CRESIL_TO_BED ( ch_cresil_raw )
                ch_versions = ch_versions.mix(CRESIL_TO_BED.out.versions)

                LONG_READ_FILTERING_CRESIL ( CRESIL_TO_BED.out.bed )
                ch_long_read_bed = ch_long_read_bed.mix(LONG_READ_FILTERING_CRESIL.out.filtered_candidates)
                ch_versions = ch_versions.mix(CRESIL_PIPELINE.out.versions)
                ch_versions = ch_versions.mix(LONG_READ_FILTERING_CRESIL.out.versions)
            }

            // ---- FLED（rca / enriched；junctions → 统一 BED，§4.6） ----
            if ( run_fled && ch_lr_rca_enriched ) {
                FLED_PIPELINE (
                    ch_lr_rca_enriched,
                    ch_fasta
                )
                .bed
                .set { ch_fled_candidates }

                LONG_READ_FILTERING_FLED ( ch_fled_candidates )
                ch_long_read_bed = ch_long_read_bed.mix(LONG_READ_FILTERING_FLED.out.filtered_candidates)
                ch_versions = ch_versions.mix(FLED_PIPELINE.out.versions)
                ch_versions = ch_versions.mix(LONG_READ_FILTERING_FLED.out.versions)
            }

            // ---- Flye（rca / enriched；组装 FASTA → remap → 统一 BED，§4.5/§8） ----
            if ( run_flye && ch_lr_rca_enriched ) {
                FLYE_PIPELINE ( ch_lr_rca_enriched )
                ch_versions = ch_versions.mix(FLYE_PIPELINE.out.versions)

                REMAP_ASSEMBLED_CIRCLES ( FLYE_PIPELINE.out.assembly, ch_fasta_meta )
                ch_versions = ch_versions.mix(REMAP_ASSEMBLED_CIRCLES.out.versions)

                LONG_READ_FILTERING_FLYE ( REMAP_ASSEMBLED_CIRCLES.out.bed )
                ch_long_read_bed = ch_long_read_bed.mix(LONG_READ_FILTERING_FLYE.out.filtered_candidates)
                ch_versions = ch_versions.mix(LONG_READ_FILTERING_FLYE.out.versions)
            }

            // ---- ecc_finder（RCA 默认路径：map → 统一 BED；asm → remap → 统一 BED，§4.1/§4.5） ----
            if ( run_eccfinder && ch_lr_rca_enriched ) {
                def eccfinder_mode = params.eccfinder_mode
                def run_map = eccfinder_mode in ['map', 'both']
                def run_asm = eccfinder_mode in ['asm', 'both']

                ECC_FINDER_ONT_SLIM (
                    ch_lr_rca_enriched,
                    ch_fasta_meta,
                    run_map,
                    run_asm
                )
                ch_versions = ch_versions.mix(ECC_FINDER_ONT_SLIM.out.versions)

                if (run_map) {
                    LONG_READ_FILTERING_ECCFINDER ( ECC_FINDER_ONT_SLIM.out.map_bed )
                    ch_long_read_bed = ch_long_read_bed.mix(LONG_READ_FILTERING_ECCFINDER.out.filtered_candidates)
                    ch_versions = ch_versions.mix(LONG_READ_FILTERING_ECCFINDER.out.versions)
                }
                if (run_asm) {
                    REMAP_ASSEMBLED_CIRCLES ( ECC_FINDER_ONT_SLIM.out.asm_fasta, ch_fasta_meta )
                    ch_versions = ch_versions.mix(REMAP_ASSEMBLED_CIRCLES.out.versions)

                    LONG_READ_FILTERING_ECCFINDER_ASM ( REMAP_ASSEMBLED_CIRCLES.out.bed )
                    ch_long_read_bed = ch_long_read_bed.mix(LONG_READ_FILTERING_ECCFINDER_ASM.out.filtered_candidates)
                    ch_versions = ch_versions.mix(LONG_READ_FILTERING_ECCFINDER_ASM.out.versions)
                }
            }

            // ---- CircleSeeker（仅 hifi + rca + concatemer=true，§3） ----
            if ( run_circleseeker ) {
                def ch_cs_input = ch_lr_rca_enriched.filter { meta, _f -> meta.assay == 'rca' && meta.concatemer && meta.read_type == 'hifi' }
                if ( ch_cs_input ) {
                    CIRCLESEEKER_PIPELINE ( ch_cs_input, ch_fasta )
                        .bed
                        .set { ch_circleseeker_candidates }

                    LONG_READ_FILTERING_CIRCLESEEKER ( ch_circleseeker_candidates )
                    ch_long_read_bed = ch_long_read_bed.mix(LONG_READ_FILTERING_CIRCLESEEKER.out.filtered_candidates)
                    ch_versions = ch_versions.mix(CIRCLESEEKER_PIPELINE.out.versions)
                    ch_versions = ch_versions.mix(LONG_READ_FILTERING_CIRCLESEEKER.out.versions)
                }
            }

            // ---- CIDER-Seq2（仅 assay=ciderseq，§4.2；可选宿主锚定） ----
            if ( run_ciderseq && ch_lr_cider ) {
                def ciderseq_cfg = new groovy.json.JsonSlurper().parseText(file(params.ciderseq_config).text)
                def ciderseq_align_genomes = ciderseq_cfg.align.targets.keySet() as List
                def ciderseq_phase_genomes = ciderseq_cfg.phase.phasegenomes.keySet() as List
                def ch_host_genome = params.ciderseq_host_genome
                    ? channel.fromPath(params.ciderseq_host_genome).map { fa -> [[ id: 'host' ], fa] }
                    : channel.empty()

                CIDERSEQ_PIPELINE (
                    ch_lr_cider,
                    file(params.ciderseq_config, checkIfExists: true),
                    file(params.ciderseq_blastdb, checkIfExists: true),
                    file(params.ciderseq_align_targets, checkIfExists: true),
                    file(params.ciderseq_protein_db, checkIfExists: true),
                    ciderseq_align_genomes,
                    ciderseq_phase_genomes,
                    ch_host_genome
                )
                // 有宿主坐标才进统一 BED（§4.2）
                ch_long_read_bed = ch_long_read_bed.mix(CIDERSEQ_PIPELINE.out.host_bed)
                ch_versions = ch_versions.mix(CIDERSEQ_PIPELINE.out.versions)
            }
        }
    }

    // 可选 organelle_tag：打 origin 列，默认不丢（§4.7）
    if ( params.filter_organelle && params.organelle_genome && ch_long_read_bed ) {
        ORGANELLE_TAG ( ch_long_read_bed, ch_fasta_meta )
        ch_long_read_bed = ORGANELLE_TAG.out.tagged
        ch_versions = ch_versions.mix(ORGANELLE_TAG.out.versions)
    }

    // 统一长读候选 BED：把该样本所有引擎的 filtered BED 合并成单一清单。
    // 引擎专属文件仍各自发布，这里额外产出一份跨引擎汇总（§4.6 BED 契约）。
    if ( ch_long_read_bed ) {
        COLLECT_LONG_READ_BED (
            ch_long_read_bed.groupTuple(by: [0])
        )
        ch_versions = ch_versions.mix(COLLECT_LONG_READ_BED.out.versions)
    }

    //
    // MODULE: Pipeline reporting
    //
    ch_versions_for_multiqc = softwareVersionsToYAML(ch_versions)
        .collectFile(name: 'software_versions_mqc.yml')

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        ch_multiqc_config                     = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config              = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
        ch_multiqc_logo                       = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()
        summary_params                        = paramsSummaryMap(workflow)
        ch_workflow_summary                   = channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = channel.empty()
            .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
            .mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
            .mix(ch_versions_for_multiqc)
            .mix(ch_fastqc_multiqc.map { _meta, zip -> zip })
            .mix(ch_trimgalore_multiqc.map { _meta, zip -> zip })
            .mix(ch_trimgalore_multiqc_log.map { _meta, log -> log })
            .mix(ch_samtools_stats.map { _meta, stats -> stats })
            .mix(ch_samtools_flagstat.map { _meta, flagstat -> flagstat })
            .mix(ch_samtools_idxstats.map { _meta, idxstats -> idxstats })
            .mix(ch_markduplicates_stats.map { _meta, stats -> stats })
            .mix(ch_markduplicates_flagstat.map { _meta, flagstat -> flagstat })
            .mix(ch_markduplicates_idxstats.map { _meta, idxstats -> idxstats })
            .mix(ch_markduplicates_multiqc.map { _meta, metrics -> metrics })
            .mix(ch_mosdepth_multiqc)

        MULTIQC (
            ch_multiqc_files.collect().ifEmpty([]),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// No top-level statements allowed here in strict DSL2.
// The onComplete and onError event handlers must be moved to main.nf or be placed outside module scope if not causing errors.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
