# Simple script to convert bam file into fragments file format
# TODO: additional options, e.g.  mapQ threshold, chromosome filtering

import sys
import pysam

infile_name = sys.argv[1]
outfile_name = sys.argv[2]

bamfile = pysam.AlignmentFile(infile_name, "rb")
outfile = open(outfile_name, 'w')

for read in bamfile.fetch():
	# Filter by mapQ
	# This should be a parameter; for now we use 0 by default
  # Note: this is different from sinto's default (30) since we are much more likely to get multi-mappers with single-end reads

	if (read.mapq == 0):
		continue

	# Here we assume that the number of duplicate reads per fragment is always 1 since we filter out duplicate UMIs
	outfile.write(f"{read.reference_name}\t{read.reference_start}\t{read.reference_end+1}\t{read.get_tag('CB')}\t1\n")

outfile.close()