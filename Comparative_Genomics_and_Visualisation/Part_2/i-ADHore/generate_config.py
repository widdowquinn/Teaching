#!/usr/bin/env python
#
# generate_config.py
#
# This script takes GenBank files representing bacterial chromosomes, and
# writes a list of genes for each with strand indicated, suitable for input
# to i-ADHoRe.
#
# This script also writes the BLAST table for i-ADHoRe, and we use RBBH
# results, rather than one-way BLAST hits, taking data from an earlier
# exercise
 
###
# IMPORTS
from Bio import SeqIO
import os
 
###
# GLOBALS
#
# These are mostly parameters that you can change, for example when undertaking
# the course exercise.

# parameters for the i-ADHoRe config file
parameters = ['output_path= i-ADHoRe_activity',
              'prob_cutoff= 0.001',
              'write_stats= true',
              'number_of_threads= 7',
              'multiple_hypothesis_correction= bonferroni',
              'gap_size= 5',
              'cluster_gap= 5',
              'anchor_points= 5',
              'alignment_method= gg2',
              'visualizeGHM= true',
              'visualizeAlignment= true',
              'q_value= 0.9',
              'max_gaps_in_alignment= 5',
              'tandem_gap= 5']

# directory containing GenBank files for analysis
gbkdir = 'data'

# accession numbers and file locations
genbank_files = {'ECA': os.path.join(gbkdir, 'NC_004547.gbk'),
                 'Pecwa': os.path.join(gbkdir, 'NC_013421.gbk'),
                 'ETA': os.path.join(gbkdir, 'NC_010694.gbk')
                 }

# file containing all RBBH data, e.g. .crunch output suitable for ACT
rbbh_data = 'data/rbbh_data.tab'

# file containing the RBBH data to be input to i-ADHoRe
blast_table = 'rbbh_input.tab'

# name of the i-ADHoRe config file to be produced
ini_filename = 'i-ADHoRe_config.ini'
 
# Open ini file for writing'
inifh = open(ini_filename, 'w')
 
# Loop over each GBK file, and write out a genome list for each to 
# the genome_lists directory (which is created if it doesn't exist), keyed by 
# the genbank_files key for that file
if not os.path.isdir('genome_lists'):
    os.mkdir('genome_lists')
 
###
# GENERATE FILES

# Parse the GenBank files and write the corresponding .lst file for 
# each genome
ftdict = {}
for k, v in genbank_files.items():
    data = SeqIO.read(v, 'genbank')
    outfilename = os.path.join('genome_lists', k + '.lst')
    print "Writing genome list for %s to %s" % (k, outfilename)
    print >> inifh, "genome= %s" % k
    outfh = open(outfilename, 'w')
    print >> inifh, "%s %s\n" % (k, outfilename)
    for feature in [f for f in data.features if f.type == 'CDS']:
        strand = '-' if feature.strand == -1 else '+'
        try:
            print >> outfh, "%s%s" % (feature.qualifiers['locus_tag'][0], strand)
            ftdict[feature.qualifiers['protein_id'][0]] = \
                   feature.qualifiers['locus_tag'][0]
        except KeyError:
            print >> outfh, "%s%s" % (feature.qualifiers['protein_id'][0], strand)            
    outfh.close()
  
def trim_ft(ft):
    if not '|' in ft:
        return ft
    if ft.endswith('|'):
        return ftdict[ft.split('|')[-2]]
    else:
        return ftdict[ft.split('|')[-1]]

# Process the RBBH data to make identifiers match the GenBank file
with open(rbbh_data, 'rU') as infh:
    with open(blast_table, 'w') as outfh:
        for ft1, ft2 in [tuple(l.strip().split()) for l in infh]:
            outfh.write('\t'.join([trim_ft(ft1), trim_ft(ft2)]) + '\n')

# Write parameters to ini file
print >> inifh, "\nblast_table= %s" % blast_table
print >> inifh, '\n' + '\n'.join(parameters)
 
# Close ini file
inifh.close()
