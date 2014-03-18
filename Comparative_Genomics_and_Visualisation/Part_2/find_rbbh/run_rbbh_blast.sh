#!/usr/bin/env bash
#
# run_rbbh_blast.sh
#
# Short script to run example RBBH BLAST comparisons for UoD teaching course
# materials. Reciprocal BLASTP searches are run on three input protein 
# sequence sets.
#
# (c) The James Hutton Institute 2014
# Author: Leighton Pritchard

# Define input files and input/output directories
filenames="NC_004547.faa NC_013421.faa NC_010694.faa"
indir="data"
outdir="rbbh_output"

# Make the output directory, if needed
mkdir -p ${outdir}

# Make BLAST databases
echo "Making databases"
for f in ${filenames}
do
    cmd="time makeblastdb -in ${indir}/${f} -dbtype prot -out ${outdir}/${f}"
    echo ${cmd}
    eval ${cmd}
done

# Loop over files and run BLASTP on all pairs of files in all combinations
for f in ${filenames}
do
    for g in ${filenames}
    do
	if [ "${f}" != "${g}" ]; then
	    echo "Running BLASTP..."
	    cmd="time blastp -query ${indir}/${f} -db ${outdir}/${g} -outfmt \"6 qseqid sseqid qlen slen length nident pident evalue bitscore\" -out ${outdir}/${f%.faa}_vs_${g%.faa}.tab -max_target_seqs 1"
	    echo ${cmd}
	    eval ${cmd}
	fi
    done
done