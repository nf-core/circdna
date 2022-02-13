#!/usr/bin/env python
# requires circle map installed
from circlemap.Coverage import coverage
import pysam as ps
import os
import pybedtools as bt
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', metavar='',
                    help="Input bam file")
# parser.add_argument('-o', '--output', metavar='',
#                                   help="Output bam file")
parser.add_argument('-b', '--bed', metavar='',
                    help="CIRCexplorer parse output in bed format")
parser.add_argument('-q', '--mapq', metavar='',
                    help="Mapping quality threshold for split reads",
                    default = 20)
parser.add_argument('--bases', metavar='',
                    help = "Number of bases for computing the coverage at junction",
                    default = 200)
parser.add_argument('-e', '--extension', metavar='',
                    help="Number of bases inside the circDNA breakpoint to comput the coverage ration",
                    default = 100)


parser.parse_args()
args = parser.parse_args(sys.argv[1:])
sorted_bam = args.input
directory = os.path.dirname(args.input)
bam_file = os.path.basename(args.input)

input_file, input_file_extension = os.path.splitext(os.path.basename(args.input))
output = directory + "/" + input_file + "_coverage.bed"

coverage_object = coverage(bam_file, bt.BedTool(args.bed), args.bases, args.mapq,
                        args.extension, directory)
bed_coverage = coverage_object.compute_coverage(coverage_object.get_wg_coverage())
bt.BedTool(bed_coverage).saveas( output )