#!/usr/bin/env python
"""
collapse_circle_alignments.py

Collapse per-query alignments (from `bedtools bamtobed`, 6-column default:
chr start end query_name mapq strand) into circular-DNA candidate loci and
emit the unified plant eccDNA BED contract:

    chr  start  end  name  score  strand  read_count  [origin]  [pair]  [engine]

Rules:
    * drop alignments with MAPQ below --min-mapq (default 20);
    * per query: cluster blocks on the same chromosome within --max-gap into a
    single locus (handles chimeric alignments / SA splits that appear as
    multiple BAM records under the same query name);
    * merge identical loci across queries (overlapping / within --max-gap) and
    sum read_count;
    * score = mean MAPQ of supporting alignments (rounded).

Usage:
    collapse_circle_alignments.py --bed in.bed --out out.bed [--min-mapq 20] [--max-gap 100]
"""
import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bed", required=True, help="Input 6-column BED from bedtools bamtobed")
    p.add_argument("--out", required=True, help="Output BED6+read_count")
    p.add_argument("--min-mapq", type=int, default=20, help="Drop alignments below this MAPQ (default: 20)")
    p.add_argument("--max-gap", type=int, default=100, help="Merge blocks/loci within this gap (default: 100)")
    return p.parse_args()


def read_bed(path, min_mapq):
    """Yield (chr, start, end, qname, mapq, strand) tuples with low-MAPQ dropped."""
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, start, end, qname, score, strand = parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4], parts[5]
            try:
                mapq = int(score)
            except ValueError:
                mapq = 0
            if mapq < min_mapq:
                continue
            yield chrom, start, end, qname, mapq, strand


def cluster_blocks(blocks, max_gap):
    """Cluster (start, end) intervals on one chromosome with max_gap; return merged intervals."""
    if not blocks:
        return []
    blocks = sorted(blocks)
    merged = [[blocks[0][0], blocks[0][1]]]
    for s, e in blocks[1:]:
        if s <= merged[-1][1] + max_gap:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def main():
    args = parse_args()
    # per-query blocks grouped by (qname, chr)
    queries = {}
    for chrom, start, end, qname, mapq, strand in read_bed(args.bed, args.min_mapq):
        queries.setdefault(qname, {}).setdefault(chrom, []).append((start, end, mapq, strand))

    # per-query merged loci
    loci = []  # (chr, start, end, qname, mean_mapq, strand)
    for qname, by_chr in queries.items():
        for chrom, blocks in by_chr.items():
            for s, e in cluster_blocks([(b[0], b[1]) for b in blocks], args.max_gap):
                maps = [b[2] for b in blocks if b[0] >= s and b[1] <= e]
                strands = [b[3] for b in blocks]
                strand = "+" if strands.count("+") >= strands.count("-") else "-"
                mean_mapq = int(round(sum(maps) / len(maps))) if maps else 0
                loci.append((chrom, s, e, qname, mean_mapq, strand))

    # merge identical loci across queries -> read_count
    loci_by_pos = {}
    for chrom, s, e, qname, mapq, strand in loci:
        key = (chrom, s, e)
        if key not in loci_by_pos:
            loci_by_pos[key] = {"read_count": 0, "mapqs": [], "qnames": [], "strand": strand}
        loci_by_pos[key]["read_count"] += 1
        loci_by_pos[key]["mapqs"].append(mapq)
        loci_by_pos[key]["qnames"].append(qname)

    # sort by chr then start, and collapse within max_gap summing read_count
    ordered = sorted(loci_by_pos.items(), key=lambda kv: (kv[0][0], kv[0][1]))
    out_rows = []
    cur = None
    for (chrom, s, e), info in ordered:
        if cur is None or chrom != cur["chrom"] or s > cur["end"] + args.max_gap:
            if cur is not None:
                out_rows.append(cur)
            cur = {"chrom": chrom, "start": s, "end": e,
                    "read_count": info["read_count"], "mapqs": list(info["mapqs"]),
                    "strand": info["strand"]}
        else:
            cur["end"] = max(cur["end"], e)
            cur["read_count"] += info["read_count"]
            cur["mapqs"].extend(info["mapqs"])
    if cur is not None:
        out_rows.append(cur)

    with open(args.out, "w") as fh:
        for i, r in enumerate(out_rows, start=1):
            mean_mapq = int(round(sum(r["mapqs"]) / len(r["mapqs"]))) if r["mapqs"] else 0
            name = "circle_{}".format(i)
            fh.write("\t".join([
                r["chrom"], str(r["start"]), str(r["end"]), name,
                str(mean_mapq), r["strand"], str(r["read_count"]),
            ]) + "\n")


if __name__ == "__main__":
    sys.exit(main())
