#!/usr/bin/env python
"""Patch CRESIL identify_wgls runtime bugs in the copied module.

1. CRESIL trim writes strand as numeric -1/1 while identify_wgls and the
    final graph logic expect '+'/'-'; normalize both numeric and '+'/'-' input.
2. Mosdepth region files are read without a header, so numeric contig names
    are ints while the summary keys are strings; cast contigs before joining.
3. If no 10bp windows pass the depth threshold, CRESIL still tries to open a
    missing bedgraph and crashes with FileNotFoundError. Emit the same ABORT
    message the module already treats as an empty result.
4. If no depth regions are within 1bp of a breakpoint window, pybedtools
    passes an empty BED file to `bedtools groupby`, which then raises
    BEDToolsError. Treat that as a valid no-detection abort too.
5. identify_wgls runs minimap2 against the full genome and uses BAI for the
    resulting BAM; use --split-prefix and CSI so large references retain @SQ
    records and can be indexed.

The module file in site-packages is read-only inside the container, so we
copy it next to the module and shadow it via PYTHONPATH. This script
performs the in-place replacement on the copy.
"""
import pathlib
import sys

REPLACEMENTS = [
    (
        "    readTrim = readTrim.loc[:,ord_header]\n",
        "    readTrim = readTrim.loc[:,ord_header]\n"
        "    readTrim['strand'] = readTrim['strand'].map(lambda x: '+' if str(x) == '1' else '-' if str(x) == '-1' else str(x))\n",
    ),
    (
        "            df[4] = df[0].apply(lambda x: dict_chrom_avg_depth.get(x, 0.0))\n",
        "            df[0] = df[0].astype(str)\n"
        "            df[4] = df[0].apply(lambda x: dict_chrom_avg_depth.get(x, 0.0))\n",
    ),
    (
        "        filtered_region_depth_bed = bt.BedTool(filtered_region_depth_path).sort()\n",
        "        if not os.path.exists(filtered_region_depth_path) or os.path.getsize(filtered_region_depth_path) == 0:\n"
        "            sys.exit(\"[ABORT] no potential merge region was detected\\n\")\n"
        "\n"
        "        filtered_region_depth_bed = bt.BedTool(filtered_region_depth_path).sort()\n",
    ),
    (
        "        merged_window = sub_breaks2x_bed.window(break_min2x_merge_bed, w=1)\n"
        "\n"
        "        df_count_merged_window = merged_window.groupby(g=[1, 2, 3], c=4, o=['count']).to_dataframe()\n",
        "        merged_window = sub_breaks2x_bed.window(break_min2x_merge_bed, w=1)\n"
        "\n"
        "        if not os.path.exists(merged_window.fn) or os.path.getsize(merged_window.fn) == 0:\n"
        "            sys.exit(\"[ABORT] no potential merge region was detected\\n\")\n"
        "\n"
        "        df_count_merged_window = merged_window.groupby(g=[1, 2, 3], c=4, o=['count']).to_dataframe()\n",
    ),
    (
        "        cmd = \"minimap2 -t {} --no-long-join -a {} {} | samtools sort -o {}/ref_aln.bam; samtools index {}/ref_aln.bam\".format(threads, fref, fastaName, tmpDir, tmpDir)\n",
        "        cmd = \"minimap2 -t {} --no-long-join --split-prefix {}/ref_aln -a {} {} | samtools sort -o {}/ref_aln.bam; samtools index -c {}/ref_aln.bam\".format(threads, tmpDir, fref, fastaName, tmpDir, tmpDir)\n",
    ),
]


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit("usage: patch_wgls.py <cresil/cli/identify_wgls.py>")
    target = pathlib.Path(sys.argv[1])
    src = target.read_text()
    for old, new in REPLACEMENTS:
        if old not in src:
            sys.exit(f"error: pattern not found in {target}: {old!r}")
        src = src.replace(old, new)
    target.write_text(src)
    print(f"patched {target}")


if __name__ == "__main__":
    main()
