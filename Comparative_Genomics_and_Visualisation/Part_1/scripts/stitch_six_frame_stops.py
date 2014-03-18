#!/usr/bin/env python
#
# stitch_six_frame_stops.py
#
# Takes an input (multiple) FASTA sequence file, and replaces all runs of
# N with the sequence NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN, which
# contains start and stop codons in all frames.  All the sequences in the
# input file are then stitched together with the same sequence.
#
# Overall, the effect is to replace all regions of base uncertainty with
# the insert sequence, forcing stops and starts for gene-calling, to avoid
# chimeras or frame-shift errors due to placement of Ns.
#
# This script is intended for use in assembly pipelines, where contigs are
# provided in the correct (or, at least, an acceptable) order.
#
# Updated 21/7/11 to produce a GFF file describing contig locations on the
# stitched assembly
#
# (C) Leighton Pritchard 2011

###
# IMPORTS

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from optparse import OptionParser

import matplotlib

import logging
import logging.handlers
import os
import re
import sys
import time

separator = 'NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN'

###
# FUNCTIONS

# Parse command-line
def parse_cmdline(args):
    """ Parse command-line arguments
    """
    usage = "usage: %prog [options] <organism_type> <infile>"
    parser = OptionParser(usage)
    parser.add_option("-o", "--outfile", dest="outfilename",
                      action="store", default=None,
                      help="Output filename")
    parser.add_option("-i", "--infile", dest="infilename",
                      action="store", default=None,
                      help="Input filename")
    parser.add_option("--id", dest="seqid",
                      action="store", default="stitched",
                      help="ID/Accession for the output stitched sequence")
    parser.add_option("-v", "--verbose", dest="verbose",
                      action="store_true", default=False,
                      help="Give verbose output")
    return parser.parse_args()

# Replace runs of N in each sequence with the separator
def stitch_ns(sequences):
    """ Loop over each input sequence in sequences, and replace any
        occurrences of N with the separator sequence.
    """
    new_sequences = []
    for s in sequences:
        seqdata = str(s.seq)
        ncount = seqdata.count('n') + seqdata.count('N')
        logger.info("%d Ns in %s" % (ncount, s.id))
        if 0 == ncount:
            logger.info("Keeping unmodified %s" % s.id)
            new_sequences.append(s)
            continue
        new_seq, repcount = re.subn('[nN]{1,}', separator, seqdata)
        logger.info("Replaced %d N runs in %s" % (repcount, s.id))
        new_seqrecord = SeqRecord(Seq(new_seq.upper()), id=s.id + "_N-replaced",
                                       name=s.name + "_N-replaced",
                                       description=s.description)
        logger.info("New SeqRecord created:\n%s" % new_seqrecord)
        new_sequences.append(new_seqrecord)
    return new_sequences

# Stitch passed sequences together with the separator
def stitch_seqs(sequences, seqid):
    """ Stitch the passed sequences together, giving the result the passed ID,
        but a description that derives from each of the passed sequences
    """
    new_seq = separator.join([str(s.seq) for s in sequences])
    new_id = seqid
    new_name = seqid+"_stitched"
    new_desc = '+'.join([s.id for s in sequences])
    stitched_seq = SeqRecord(Seq(new_seq), id=new_id, name=new_name,
                             description=new_desc)
    logger.info("Created stitched sequence (len:%d):\n%s" % \
                (len(stitched_seq), stitched_seq))
    return stitched_seq

# Generate GFF file corresponding to the stitched sequence
def build_gff(sequences, seqid):
    """ Loop over the passed sequences, and generate a GFF file
        describing the locations of the contigs on the stitched
        chromosome.  We do this last, to account for the possible
        introduction of internal separator sequence where (PE)
        assemblies introduce Ns.
    """
    gff_out = ['##gff-version 3']
    start = 1
    for s in sequences:
        if start == 1:
            end = len(s)
        else:
            end = start + len(s)
        gff_out.append('\t'.join([seqid, 'stitching', 'contig', str(start),
                                  str(end), '.', '.', '.',
                                  "ID=%s;Name=%s" % (s.id, s.id)]))
        start = end + len(separator)
    return '\n'.join(gff_out) + '\n'

###
# SCRIPT

if __name__ == '__main__':

    # Parse command-line
    # options are options, arguments are the .sff files
    options, args = parse_cmdline(sys.argv)
    
    # We set up logging, and modify loglevel according to whether we need
    # verbosity or not
    logger = logging.getLogger('stitch_six_frame_stops.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = \
                  logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    if options.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info(options)
    logger.info(args)

    # Load data from file (STDIN if no filename provided)
    if options.infilename is None:
        inhandle = sys.stdin
        logger.info("Using STDIN for input")
    else:
        inhandle = open(options.infilename, 'rU')
    data = list(SeqIO.parse(inhandle, 'fasta'))
    inhandle.close()

    # Stitch individual sequences
    data = stitch_ns(data)

    # Stitch all sequences together
    stitchedseq = stitch_seqs(data, options.seqid)

    # Generate GFF output
    gff_data = build_gff(data, options.seqid)

    # Write the stitched sequence to file (or STDOUT if no filename provided)
    if options.outfilename is None:
        outhandle = sys.stdout
    else:
        outhandle = open(options.outfilename, 'w')
    SeqIO.write([stitchedseq], outhandle, 'fasta')
    outhandle.close()

    # Write GFF output
    if options.outfilename is None:
        gff_filename = stitchedseq.id+'_contigs.gff'
    else:
        gff_filename = os.path.splitext(options.outfilename)[0]+'_contigs.gff'
    logger.info("Writing GFF to %s" % gff_filename)
    gffhandle = open(gff_filename, 'w')
    gffhandle.write(gff_data)
    gffhandle.close()
    

    

