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
parser.add_argument('-d', '--directory', metavar='',
                                  help="Directory in which output is saved", 
                                  default = "./")
parser.add_argument('-o', '--output', metavar='',
                                  help="Output bed file with coverage information")

parser.parse_args()
args = parser.parse_args(sys.argv[1:])
sorted_bam = args.input
directory = os.path.dirname(args.input)
bam_file = os.path.basename(args.input)

coverage_object = coverage(bam_file, bt.BedTool(args.bed),300,20,150, directory)
bed_coverage = coverage_object.compute_coverage(coverage_object.get_wg_coverage())
bt.BedTool(bed_coverage).saveas( args.output )