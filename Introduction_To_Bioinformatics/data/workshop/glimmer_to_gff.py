# glimmer_to_gff.py
#
# Ad hoc script that converts the default Glimmer output to a GFF format file
#
# L.Pritchard 2011

import sys
import os

# We take the filename as the only argument, and generate a gff file with name based
# on the input filename

filename = sys.argv[-1]

data = [l.strip().split() for l in open(filename, 'rU').readlines()]

# The first line is a header, from which we obtain the seqid, by removing the right
# angled bracket (>)

seqid = data[0][0][1:]

# The remaining lines are in columns:
# ID, start, end, strand, score
# We need to reorder start/end in numerical order, and convert strand to +/- only

outfilename = os.path.splitext(filename)[0] + '.gff'
outfh = open(outfilename, 'w')

for line in data[1:]:
    start, end = int(line[1]), int(line[2])
    if start > end:
        start, end = end, start
    strand = line[3][0]
    print >> outfh, '\t'.join([seqid, 'Glimmer', 'CDS', str(start), str(end),
                              line[4], strand, '0', 'ID=%s;Name=%s' % (line[0], line[0])])

outfh.close()
