#!/usr/bin/env bash
#
# run_BLAST.sh
#
# Short schell script to run all-vs-all BLAST job for example MCL
# clustering to find putative orthologues in three bacteria. This
# is for an activity in a UoD training workshop/course on 
# comparative genomics.

# Grab the protein sequences from the find_rbbh activity
cat ../find_rbbh/data/*.faa > data/proteins.faa

# Build the BLAST protein database
makeblastdb -in data/proteins.faa -dbtype prot

# Run the BLASTP search (default parameters)
blastp -query data/proteins.faa -db data/proteins.faa -out data/all-vs-all.tab -outfmt 6 -evalue 1e-30
