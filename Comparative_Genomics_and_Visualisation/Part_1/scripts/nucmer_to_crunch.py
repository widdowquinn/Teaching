#!/usr/bin/env python
#
# nucmer_to_crunch.py
#
# USAGE: Usage: nucmer_to_crunch.py [-h] [-o OUTFILENAME] [-i INFILENAME] [-v]
#
# optional arguments:
#   -h, --help            show this help message and exit
#   -o OUTFILENAME, --outfile OUTFILENAME
#                         Output .crunch file
#   -i INFILENAME, --infile INFILENAME
#                         Input .coords file
#   -v, --verbose         Give verbose output
#
# A short script that converts the output of MUMmer's show-coords package
# to a .crunch file that can be used in Sanger's ACT comparative genomics
# visualisation tool.
#
# Copyright (C) 2014 The James Hutton Institute
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@hutton.ac.uk
#
# Leighton Pritchard,
# Information and Computing Sciences,
# James Hutton Institute,
# Errol Road,
# Invergowrie,
# Dundee,
# DD6 9LH,
# Scotland,
# UK
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

###
# IMPORTS

from argparse import ArgumentParser

import logging
import logging.handlers
import os
import re
import sys


###
# FUNCTIONS

# Parse cmd-line
def parse_cmdline(args):
    """ Parse command-line arguments. Note that the input and output
        directories are positional arguments
    """
    parser = ArgumentParser(prog="nucmer_to_crunch.py")
    parser.add_argument("-o", "--outfile", dest="outfilename",
                        action="store", default=None,
                        help="Output .crunch file")
    parser.add_argument("-i", "--infile", dest="infilename",
                        action="store", default=None,
                        help="Input .coords file")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true",
                        help="Give verbose output")
    return parser.parse_args()


# Report last exception as string
def last_exception():
    """ Returns last exception as a string
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# Process the input stream
def process_stream(infh, outfh):
    """ Processes the input stream, assuming show-coords output, with
        five header lines, and whitespace separation.

        show-coords output has the following columns (post-processing)

        1: reference sequence start
        2: reference sequence end
        3: subject sequence start
        4: subject sequence end
        5: reference alignment length
        6: subject alignment length
        7: alignment percentage identity
        8: reference sequence ID
        9: subject sequence ID

        This is converted to .crunch (MSPcrunch) format output with the
        following columns, separated by whitespace.

        1: score (reference alignment length)
        2: alignment percentage identity
        3: reference sequence start
        4: reference sequence end
        5: reference sequence ID
        6: subject sequence start
        7: subject sequence end
        8: subject sequence ID
    """
    # Read in the input stream into a list of lines
    try:
        tbldata = list(infh.readlines())
    except:
        logger.error("Could not process input (exiting)")
        logger.error(last_exception())
        sys.exit(1)
    logger.info("Read %d lines from input" % len(tbldata))

    tbldata = tbldata[5:]
    logger.info("Skipped header lines.")

    for line in [l.strip().split() for l in tbldata if
                 len(l.strip())]:
        # Due to the column marker symbols, there are offsets for the
        # columns, relative to those in the text above.
        outfh.write(' '.join([line[6], line[9], line[0], line[1],
                              line[11], line[3], line[4], line[12]]) + '\n')


###
# SCRIPT

if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline(sys.argv)

    # We set up logging, and modify loglevel according to whether we need
    # verbosity or not
    logger = logging.getLogger('nucmer_to_crunch.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info(args)

    # Do we have an input file?  No? Then use stdin
    if args.infilename is None:
        infhandle = sys.stdin
        logger.info("Using stdin for input")
    else:
        logger.info("Using %s for input" % args.infilename)
        try:
            infhandle = open(args.infilename, 'rU')
        except:
            logger.error("Could not open input file: %s (exiting)" %
                         args.infilename)
            logger.error(''.join(
                traceback.format_exception(sys.last_type,
                                           sys.last_value,
                                           sys.last_traceback)))
            sys.exit(1)

    # Do we have an output file?  No? Then use stdout
    if args.outfilename is None:
        outfhandle = sys.stdout
        logger.info("Using stdout for output")
    else:
        logger.info("Using %s for output" % args.outfilename)
        try:
            outfhandle = open(args.outfilename, 'w')
        except:
            logger.error("Could not open output file: %s (exiting)" %
                         args.outfilename)
            logger.error(''.join(
                traceback.format_exception(sys.last_type,
                                           sys.last_value,
                                           sys.last_traceback)))
            sys.exit(1)

    # Process input stream
    process_stream(infhandle, outfhandle)
