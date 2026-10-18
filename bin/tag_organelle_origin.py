#!/usr/bin/env python
"""
tag_organelle_origin.py

Annotate eccDNA candidate BED entries with an origin column:
    origin = nuclear | pt | mt | ambiguous

Inputs:
    --candidates  unified eccDNA BED (chr start end name score strand [read_count ...])
    --organelle   BED of organelle-origin regions in the SAME coordinate system as the
                candidates (e.g. produced by mapping the organelle genome to the
                nuclear reference). Optional column 4 may label the type: mt / pt /
                ambiguous (defaults to 'ambiguous').
    --out         output BED = candidate rows + '\t' + origin
    --min-overlap fraction of the candidate that must overlap organelle regions to be
                tagged (default 0.5)

Defaults to NOT dropping organelle candidates; the caller decides whether to
filter them out afterwards (--drop handled by the workflow param).
"""
import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidates", required=True)
    p.add_argument("--organelle", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--min-overlap", type=float, default=0.5)
    return p.parse_args()


def load_bed(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            rows.append(parts)
    return rows


def main():
    args = parse_args()
    candidates = load_bed(args.candidates)
    organelle = load_bed(args.organelle)

    # index organelle regions by chromosome (sorted for sweep)
    org_by_chr = {}
    for parts in organelle:
        chrom = parts[0]
        try:
            s, e = int(parts[1]), int(parts[2])
        except ValueError:
            continue
        label = "ambiguous"
        if len(parts) >= 4 and parts[3]:
            lab = parts[3].strip().lower()
            if lab in ("mt", "pt"):
                label = lab
        org_by_chr.setdefault(chrom, []).append((s, e, label))
    for chrom in org_by_chr:
        org_by_chr[chrom].sort()

    def organelle_overlap(chrom, s, e):
        """Return (overlap_bp, label) of best-overlapping organelle region."""
        best = (0, "ambiguous")
        for os_, oe, label in org_by_chr.get(chrom, []):
            if oe <= s:
                continue
            if os_ >= e:
                break
            ov = min(e, oe) - max(s, os_)
            if ov > best[0]:
                best = (ov, label)
        return best

    with open(args.out, "w") as fh:
        for parts in candidates:
            chrom, s, e = parts[0], int(parts[1]), int(parts[2])
            ov, label = organelle_overlap(chrom, s, e)
            if ov >= (e - s) * args.min_overlap:
                origin = label
            else:
                origin = "nuclear"
            fh.write("\t".join(parts + [origin]) + "\n")


if __name__ == "__main__":
    sys.exit(main())
