
workflow CIRCDNA_ASSEMBLY {
    take:
        reads
        fasta

    main:
    ch_trimmed_reads.map{ meta, file -> [meta, file, []] }.set{ch_unicycler_input}
}
