#!/usr/bin/env python
"""Merge per-engine long-read candidate BEDs into one unified BED.

The long-read branch produces one filtered BED per detection engine
(CReSIL / FLED / Flye / ecc_finder / CircleSeeker). They all follow the
shared contract

    chr  start  end  name  score  strand  read_count [te_overlap] [origin]

but they are separate files. Downstream consumers (scoring, cross-species
comparison) want a single per-sample candidate list, so concatenate them,
keep exactly one header, and sort by position.

Rows are copied verbatim: the optional trailing columns are engine-specific
and are preserved rather than normalised.
"""
import argparse
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, help='Unified BED to write')
    parser.add_argument('beds', nargs='+', help='Per-engine filtered BED files')
    return parser.parse_args()


def is_header(line):
    """Header is either a '#' comment or the column-name row."""
    if line.startswith('#'):
        return True
    first = line.split('\t', 1)[0].strip().lower()
    return first in ('chr', 'chrom', 'chromosome', '#chr')


def natural_key(text):
    """Sort 1, 2, 10 before 1A/2B-style names and keep 'Un' last."""
    parts = re.split(r'(\d+)', text)
    key = []
    for part in parts:
        if part.isdigit():
            key.append((0, int(part), ''))
        elif part:
            key.append((1, 0, part))
    # Contigs without a numeric component (e.g. 'Un') sort after numbered ones.
    return key


def row_key(fields):
    try:
        start = int(fields[1])
    except (IndexError, ValueError):
        start = 0
    try:
        end = int(fields[2])
    except (IndexError, ValueError):
        end = 0
    chrom = fields[0] if fields else ''
    return (natural_key(chrom), start, end)


def main():
    args = parse_args()
    header = None
    rows = []

    for path in args.beds:
        with open(path) as handle:
            for line in handle:
                if not line.strip():
                    continue
                if is_header(line):
                    if header is None:
                        # Keep the canonical column names, not '#' comments.
                        if not line.startswith('#'):
                            header = line.rstrip('\n')
                    continue
                rows.append(line.rstrip('\n'))

    if header is None:
        header = '\t'.join(['chr', 'start', 'end', 'name', 'score',
                            'strand', 'read_count'])

    rows.sort(key=lambda line: row_key(line.split('\t')))

    with open(args.out, 'w') as out:
        out.write(header + '\n')
        for row in rows:
            out.write(row + '\n')

    print('unified BED: {} rows from {} file(s) -> {}'.format(
        len(rows), len(args.beds), args.out), file=sys.stderr)


if __name__ == '__main__':
    main()
