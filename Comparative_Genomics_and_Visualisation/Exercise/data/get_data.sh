#!/usr/bin/env bash
#
# Run this script to download data files for the 
# Exercise

REMDIR="ftp://ftp.ncbi.nih.gov/genomes/Bacteria"

stems="Staphylococcus_epidermidis_ATCC_12228_uid57861/NC_004461
       Pseudomonas_syringae_tomato_DC3000_uid57967/NC_004578
       Pseudomonas_fluorescens_Pf0_1_uid57591/NC_007492
       Staphylococcus_aureus_NCTC_8325_uid57795/NC_007795
       Pseudomonas_putida_GB_1_uid58735/NC_010322
       Klebsiella_pneumoniae_342_uid59145/NC_011283
       Dickeya_dadantii_Ech703_uid59363/NC_012880
       Dickeya_zeae_Ech1591_uid59297/NC_012912
       Dickeya_dadantii_Ech586_uid42519/NC_013592
       Klebsiella_variicola_At_22_uid42113/NC_013850
       Staphylococcus_lugdunensis_HKU09_01_uid46233/NC_013893
       Dickeya_dadantii_3937_uid52537/NC_014500
       Klebsiella_oxytoca_KCTC_1686_uid83159/NC_016612
       Staphylococcus_aureus_ED133_uid159689/NC_017337
       Klebsiella_pneumoniae_KCTC_2242_uid162147/NC_017540
       Pseudomonas_aeruginosa_PA1_uid228931/NC_022808"

suffixes="fna faa gbk gff"

for x in ${suffixes}
do
    mkdir -p ../${x}
    for s in ${stems}
    do
	wget -nc ${REMDIR}/${s}.${x}
	ln -s ${s##.*/}.${x} ../${x}/
    done
done
