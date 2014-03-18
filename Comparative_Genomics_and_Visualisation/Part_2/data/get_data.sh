#!/usr/bin/env bash
#
# Run this script to download data files for the 
# Part 2 acivities.

REMDIR="ftp://ftp.ncbi.nih.gov/genomes/Bacteria"

# GenBank files are downloaded to the current directory, and 
# symbolic links placed in the mcl_orthologues/data 
# i-ADHoRe/data and find_rbbh/data directories
gbk="Pectobacterium_atrosepticum_SCRI1043_uid57957/NC_004547.gbk
     Erwinia_tasmaniensis_Et1_99_uid59029/NC_010694.gbk
     Pectobacterium_wasabiae_WPP163_uid41297/NC_013421.gbk"

LINKDIRS="../mcl_orthologues/data/
          ../find_rbbh/data/
          ../i-ADHore/data/"
for l in ${LINKDIRS}
do
    mkdir -p ${l}
done

for f in ${gbk}
do
    wget -nc ${REMDIR}/${f}
    for l in ${LINKDIRS}
    do
	ln -s ${f##.*/} ${l}
    done
done

# Protein FASTA files are downloaded to the current directory, and 
# symbolic links placed in the find_rbbh/data directory
faa="Pectobacterium_atrosepticum_SCRI1043_uid57957/NC_004547.faa
     Erwinia_tasmaniensis_Et1_99_uid59029/NC_010694.faa
     Pectobacterium_wasabiae_WPP163_uid41297/NC_013421.faa"

LINKDIR="../find_rbbh/data/"
mkdir -p ${LINKDIR}

for f in ${faa}
do
    wget -nc ${REMDIR}/${f}
    ln -s ${f##.*/} ${LINKDIR}
done
