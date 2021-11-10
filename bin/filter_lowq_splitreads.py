#!/usr/bin/env python
import pysam as ps
import sys
import argparse
import os
 
parser = argparse.ArgumentParser()
 
parser.add_argument('-i', '--input', metavar='',
                                  help="Input bam file")
parser.add_argument('-o', '--output', metavar='',
                                  help="Output bam file")
parser.add_argument('-q', '--quality', type=int, metavar='',
                                  help="bwa-mem mapping quality cutoff. Default value 20",
                                  default=20)
parser.parse_args()
args = parser.parse_args(sys.argv[1:])
bam = ps.AlignmentFile("%s" % args.input,"rb")
filtered = ps.AlignmentFile("%s"% args.output, "wb", template=bam)

quality_score = args.quality
 
for read in bam:
    if read.mapq < quality_score:
        continue
    else:
        if read.has_tag('SA'):
            if int(read.get_tag('SA').split(',')[4]) < quality_score:
                continue
            else:
                #write to disk
                filtered.write(read)
        else:
            #write to disk
            filtered.write(read)