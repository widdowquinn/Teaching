#!/usr/bin/env bash
#
# Run this script to download data files for the 
# Part 1 acivities.

REMDIR="ftp://ftp.ncbi.nih.gov/genomes/Bacteria"

# These bacteria are downloaded to the current directory only
bacteria="Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna
          Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.fna
          Escherichia_coli_O157_H7_Sakai_uid57781/NC_002695.fna
          Escherichia_coli_CFT073_uid57915/NC_004431.fna
          Escherichia_coli_ATCC_8739_uid58783/NC_010468.fna
          Nostoc_punctiforme_PCC_73102_uid57767/NC_010628.fna
          Dickeya_dadantii_3937_uid52537/NC_014500.fna
          Mycobacterium_tuberculosis_UT205_uid162183/NC_016934.fna
          Mycobacterium_tuberculosis_CCDC5079_uid161943/NC_017523.fna
          Mycobacterium_tuberculosis_Haarlem_uid54453/NC_022350.fna"

for f in ${bacteria}
do
    wget -nc ${REMDIR}/${f}
done

# The Mycoplasma are downloaded to the current directory, and 
# symbolic links also go into the average_nucleotide_identity/data
# directory
mycoplasma="Mycoplasma_pneumoniae_M129_uid57709/NC_000912.fna
            Mycoplasma_pneumoniae_309_uid85495/NC_016807.fna
            Mycoplasma_pneumoniae_FH_uid162027/NC_017504.fna
            Mycoplasma_genitalium_M2321_uid173373/NC_018495.fna
            Mycoplasma_genitalium_M6282_uid173371/NC_018496.fna
            Mycoplasma_genitalium_M6320_uid173370/NC_018497.fna
            Mycoplasma_genitalium_M2288_uid173372/NC_018498.fna
            Mycoplasma_pneumoniae_M129_B7_uid185759/NC_020076.fna"

LINKDIR="../average_nucleotide_identity/data/"
mkdir -p ${LINKDIR}

for f in ${mycoplasma}
do
    wget -nc ${REMDIR}/${f}
    ln -s ${f##.*/} ${LINKDIR}
done

# The Mycoplasma are downloaded to the current directory, and 
# symbolic links also go into the biopython_visualisation/data
# directory
pyrococcus="Pyrococcus_abyssi_GE5_uid62903/NC_000868.fna
            Pyrococcus_horikoshii_OT3_uid57753/NC_000961.fna
            Pyrococcus_furiosus_DSM_3638_uid57873/NC_003413.fna"

LINKDIR="../biopython_visualisation/data/"
mkdir -p ${LINKDIR}

for f in ${pyrococcus}
do
    wget -nc ${REMDIR}/${f}
    ln -s ${f##.*/} ${LINKDIR}
done
