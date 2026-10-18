#!/usr/bin/env python
"""Patch FLED 1.7.0 for long-read eccDNA runs.

FLED can write the same contig name twice into temp.PseudoReference.fa when
multiple reads support the same multi-segment junction. minimap2 then emits a
duplicate SAM header and samtools sort aborts before the FLED process finishes.

Upstream problems patched here:

1. multisegs_cluster skips the first group while clustering, so identical
    segment combinations can create duplicate groups.
2. MS_PseudoReference builds one FASTA record per Segreads key but the record
    id omits strand, so strand variants of the same junction collapse to the
    same FASTA name.
3. When a non-empty multiseg dict merges to zero usable candidates, FLED still
    realigns an empty pseudo-reference and pysam rejects the headerless BAM.
    Open alignment files with check_sq=False so this is handled as an empty
    multi-segment result.

4. FLED's internal minimap2 command can build a multi-part index on very
    large references. Pass --split-prefix so the SAM retains @SQ records.
5. ecDNAdetection.py uses BAI for the full-genome BAM. Use CSI when
    reference contigs exceed the ~512 Mb BAI limit.

The module copies the pinned FLED package into a writable directory, applies
this patch there, and shadows it via PYTHONPATH.
"""
import pathlib
import sys


REPLACEMENTS = [
    (
        "        for group in groups[1:] :\n",
        "        for group in groups :\n",
    ),
    (
        "            if inGroup == True :\n"
        "                newGroup = False\n"
        "                group.append(x)\n",
        "            if inGroup == True :\n"
        "                newGroup = False\n"
        "                group.append(x)\n"
        "                break\n",
    ),
    (
        "    ecclist = []\n"
        "    for record in Segreads:\n",
        "    ecclist = []\n"
        "    seen_pseudo = set()\n"
        "    for record in Segreads:\n",
    ),
    (
        "        ecclist.append(eccref)\n",
        "        if eccid in seen_pseudo:\n"
        "            print(datetime.datetime.now().strftime(\"\\n%Y-%m-%d %H:%M:%S:\"), "
        "\"Skipped duplicate PseudoReference junction\", eccid)\n"
        "        else:\n"
        "            seen_pseudo.add(eccid)\n"
        "            ecclist.append(eccref)\n",
    ),
    (
        "    bamFile = ps.AlignmentFile(\"%s\" % ont_bam, \"rb\")\n",
        "    bamFile = ps.AlignmentFile(\"%s\" % ont_bam, \"rb\", check_sq=False)\n",
    ),
]


ECDNA_REPLACEMENTS = [
    (
        "        samtoolsIndex = subprocess.call([\"samtools\", \"index\",bamfile, baifile], shell=False)\n",
        "        samtoolsIndex = subprocess.call([\"samtools\", \"index\", \"-c\", bamfile], shell=False)\n",
    ),
    (
        "        alignment = subprocess.call([\"minimap2\", \"-t\", str(self.threads), \"-ax\", \"map-ont\", self.reffa, self.input_fq, \"-o\", samfile], shell=False)\n",
        "        split_prefix = self.out_dir + '/MappingResult/' + self.label + \".mmi\"\n"
        "        alignment = subprocess.call([\"minimap2\", \"-t\", str(self.threads), \"-ax\", \"map-ont\", \"--split-prefix\", split_prefix, self.reffa, self.input_fq, \"-o\", samfile], shell=False)\n",
    ),
]

def apply_replacements(target: pathlib.Path, replacements: list) -> None:
    src = target.read_text()
    for old, new in replacements:
        if old not in src:
            sys.exit(f"error: pattern not found in {target}: {old!r}")
        src = src.replace(old, new)
    target.write_text(src)
    print(f"patched {target}")


def main() -> None:
    if len(sys.argv) < 2:
        sys.exit("usage: patch_fled.py <FLED/_utils.py> [FLED/ecDNAdetection.py]")
    for target_arg in sys.argv[1:]:
        target = pathlib.Path(target_arg)
        if target.name == "_utils.py":
            apply_replacements(target, REPLACEMENTS)
        elif target.name == "ecDNAdetection.py":
            apply_replacements(target, ECDNA_REPLACEMENTS)
        else:
            sys.exit(f"error: unsupported patch target: {target}")


if __name__ == "__main__":
    main()
