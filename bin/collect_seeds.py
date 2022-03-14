#!/usr/bin/env python

# Code adapted from PrepareAA (https://github.com/jluebeck/PrepareAA)
# commit/92b81ba55356af85958985b8f80308c8f88921ac


import argparse
from datetime import datetime
from subprocess import call

# Read the CNVkit .cns files
def collect_seeds(sample, cns):
    with open(cns) as infile, open(sample + "_CNV_GAIN.bed", 'w') as outfile:
        head = next(infile).rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            s, e = int(fields[1]), int(fields[2])
            cn_r = float(fields[4])
            cn = 2 ** (cn_r + 1)
            if cn >= args.cngain:  # do not filter on size since amplified_intervals.py will merge small ones.
                outline = "\t".join(fields[0:3] + ["CNVkit", str(cn)]) + "\n"
                outfile.write(outline)
    return sample + "_CNV_GAIN.bed"

# MAIN #
if __name__ == '__main__':
    # Parses the command line arguments
    parser = argparse.ArgumentParser(
        description="Collect AmpliconArchitect Copy Number Seeds")
    parser.add_argument("-s", "--sample", help="sample name", required=True)
    parser.add_argument("--cns", help="CNVKit .cns file of CNV changes.", default="")
    parser.add_argument("--cngain", type=float, help="CN gain threshold to consider for AA seeding", default=4.5)
    args = parser.parse_args()
    collect_seeds(args.sample, args.cns)

