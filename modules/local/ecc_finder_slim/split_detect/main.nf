process ECC_FINDER_SPLIT_DETECT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ecc_finder_slim:1.0.0--0' :
        'quay.io/bioinfortools/ecc_finder_slim:1.0.0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.split.bed"), emit: split_bed
    tuple val(meta), path("${prefix}.disc.bed"),  emit: disc_bed
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # BAM → BED (6 cols: chrom, start, end, name, score, strand)
    bedtools bamtobed -i ${bam} > ${prefix}.bam.bed

    # Split-read + discordant-read detection (replicates ecc_finder map-sr.py run_split/run_disc)
    python <<PYEOF
    import pandas as pd
    from pybedtools import BedTool

    df = pd.read_csv('${prefix}.bam.bed', sep='\t', header=None,
                    names=['chrom','start','end','name','score','strand'])
    # name is 'read/1' or 'read/2' → split into read + pair number
    tmp = df['name'].astype(str).str.rsplit('/', n=1, expand=True)
    df['read'] = tmp[0]
    df['pair'] = tmp[1]
    df = df.drop(columns=['name','score'])
    # key by chrom__read
    df['chrom'] = df['chrom'] + '__' + df['read']
    df = df[['chrom','start','end','pair','strand']]

    # Merge overlapping alignments per read (collapse pair + strand)
    merged = BedTool.merge(
        BedTool.sort(BedTool.from_dataframe(df)),
        c='4,5', o='collapse,collapse'
    ).groupby(g=[1], c=[2,3,4,5], o=['min','max','collapse','collapse'])
    dM = BedTool.to_dataframe(merged, names=['chrom','start','end','pair','strand'])

    def write_bed(sel, outfile):
        if sel.empty:
            pd.DataFrame(columns=['chrom','start','end','read','length']).to_csv(
                outfile, header=False, index=False, sep='\t')
            return
        sel[['chrom','read']] = sel['chrom'].str.split('__', n=1, expand=True)
        sel['length'] = sel['end'] - sel['start']
        sel = sel.sort_values(by=['chrom','start'])
        sel[['chrom','start','end','read','length']].to_csv(
            outfile, header=False, index=False, sep='\t')

    # Split reads: 3 unique hits, same chromosome, specific orientation pattern
    splitM = dM[
        ((dM.pair == '1,2,1') & (dM.strand == '+,-,+')) |
        ((dM.pair == '2,1,2') & (dM.strand == '-,+,-')) |
        ((dM.pair == '1,2,1') & (dM.strand == '-,+,-')) |
        ((dM.pair == '2,1,2') & (dM.strand == '+,-,+'))
    ].copy()
    write_bed(splitM, '${prefix}.split.bed')

    # Discordant reads: outward facing read pair
    discM = dM[
        ((dM.pair == '1,2') & (dM.strand == '-,+')) |
        ((dM.pair == '2,1') & (dM.strand == '+,-'))
    ].copy()
    write_bed(discM, '${prefix}.disc.bed')
    PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ecc_finder_slim: 1.0.0
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.split.bed
    touch ${prefix}.disc.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ecc_finder_slim: 1.0.0
    END_VERSIONS
    """
}
