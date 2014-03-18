# README.md - average_nucleotide_identity

## Overview

This activity is intended to guide you through calculation of Average Nucleotide Identity measures, using **JSpecies** and a Python script.

The purpose of this activity is

1. to demonstrate and give practice in calculation of average nucleotide identity (ANI) measures for prokaryotic sequences.
2. to give practice in interpretation of ANI output.

When you have completed this activity you should be able 

1. to use JSpecies and/or a Python script to calculate ANI measures for prokaryotes.
2. to interpret ANI calculation output

## Prerequisites

This exercise was written and tested using the following software, but other versions may also be suitable

* **JSpecies** 1.2.1 <http://www.imedea.uib.es/jspecies/>
* **BLAST+** 2.2.29+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>
* legacy **BLAST** 2.2.26 <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/>
* **MUMMer** 3.0 <http://mummer.sourceforge.net/>
* **Python** 2.7 <http://www.python.org>
* **Biopython** 1.63 <http://www.biopython.org>
* **RPy2** 2.3.8 <http://rpy.sourceforge.net/rpy2.html>
* **R** 3.0.2 *shared libraries* <http://www.r-project.org/>

# 1. JSpecies

## Starting JSpecies

In this `average_nucleotide_identity` subdirectory, create a subdirectory called `jspecies_output`:

```
$ mkdir jspecies_output
```

to store activity files, then run **JSpecies** with the command

```
$ java -jar -Xms1024m -Xmx1024m ../scripts/jspecies1.2.1.jar
```

## Using JSpecies to calculate ANIb, ANIm and TETRA measures

### Setup

Once JSpecies has started, you will be asked to provide some file locations for **JSpecies** to run.

![JSpecies splash warning](images/jspecies1.png?raw=True =400x)

To follow this activity, you should make the "Workspace" be the `jspecies_output` subdirectory.

**NOTE:** **JSpecies** uses the legacy version of **BLAST**, not **BLAST+** so the `blastall`, `formatDB` and `fastacmd` options needs to be set to legacy **BLAST** executables. This should have been done for you. On my laptop, the default settings look like this:

![Completed JSpecies preferences screen](images/jspecies2.png?raw=True =300x)

![JSpecies ANIb preferences](images/jspecies3.png?raw=True =300x)

![JSpecies ANIm preferences](images/jspecies4.png?raw=True =300x)

### Loading sequences from file

Create a new Group, using the `File -> New` menu item, and give it a suitable name:

![JSpecies File dialogue](images/jspecies5.png?raw=True =200x)

Now you can use the `File -> Import FASTA(S) From File(s)` menu item to load sequences.

![JSpecies Import dialogue](images/jspecies6.png?raw=True =200x)

Select the following files (they are all *Mycoplasma* chromosomes):

* data/NC_000912.fna
* data/NC_016807.fna
* data/NC_017504.fna
* data/NC_018495.fna
* data/NC_018496.fna
* data/NC_018497.fna
* data/NC_018498.fna
* data/NC_020076.fna

You should see that the sequences load, and relevant information for each genome is displayed:

![JSpecies sequence data window](images/jspecies7.png?raw=True =400x)

### Running Analyses

To run ANIb, select the `Calculation -> Average Nucleotide Identity (BLAST)` menu item. A new window will appear:

![JSpecies Calculation menu option](images/jspecies8.png?raw=True =200x)

This will allow you to select which pairs of sequences will be compared. Select all the comparisons:

![JSpecies sequence selection window](images/jspecies9.png?raw=True =400x)

and click on the `Start` button. The calculations will take a couple of minutes to complete, but soon you will see a range of summary information regarding the comparisons made for each input sequence:

![JSpecies analysis result window](images/jspecies10.png?raw=True =400x)

The same procedure, but selecting `Calculation -> Tetra Nucleotide Distribution` and `Calculation -> Average Nucleotide Identity` will give you results for TETRA and ANIm analyses independently. Alternatively, you can select all three analyses in the check boxes at the top of the `Calculation Selection` window.

![JSpecies method selection options](images/jspecies11.png?raw=True =200x)

**ACTIVITY 1** Carry out TETRA and ANIm analyses on this dataset, and inspect the results.


### Viewing results

**JSpecies** will show you output results, if you use the `Results -> Show` menu item:

![JSpecies show results menu option](images/jspecies12.png?raw=True =200x)

![JSpecies analysis output](images/jspecies13.png?raw=True =400x)

These are in the form of a table, where green entries indicate that a comparison passes the nominal 95% species boundary threshold (i.e. the compared sequences are members of the same species), and red indicates that it does not.

# 2. calculate_ani.py script

The `scripts/calculate_ani.py` script was written as an alternative to the **JSpecies** software, for use in a bioinformatics pipeline. While GUI-based analyses can be intuitive for the user, they are not well-suited to automation, as you may have found when selecting sequences for analysis.

You can run this script to carry out the same analysis as above, with **JSpecies**, using the commands:

```
$ time ../scripts/calculate_ani.py -i data/ -g --format=png -m ANIm -o ANIm
$ time ../scripts/calculate_ani.py -i data/ -g --format=png -m ANIb -o ANIb
$ time ../scripts/calculate_ani.py -i data/ -g --format=png -m TETRA -o TETRA
```

This generates plain text tables and (when the `rpy2` library is available - version incompatibilities with CentOS 6 packages mean it is not for this course! Also running on the VM, only one processor is available) output images in each of the three directories `ANIb`, `ANIm` and `TETRA`, as seen below:

![ANIb script output](images/ANIb.png?raw=True =300x)

![ANIm script output](images/ANIm.png?raw=True =300x)

![TETRA script output](images/TETRA.png?raw=True =300x)

The ANIm tabular output in `perc_ids.tab`:

```
# calculate_ani.py Sat Mar  1 19:03:51 2014
# ANIm
        NC_000912       NC_016807       NC_017504       NC_018495       NC_018496       NC_018497       NC_018498       NC_020076
NC_000912       NA      0.9951  0.9949  0.9047  0.9043  0.9039  0.9045  0.9981
NC_016807       0.9951  NA      0.9976  0.9049  0.9045  0.9041  0.9047  0.9954
NC_017504       0.9949  0.9976  NA      0.9049  0.9044  0.9041  0.9046  0.9952
NC_018495       0.9047  0.9049  0.9049  NA      0.9890  0.9914  0.9898  0.9046
NC_018496       0.9043  0.9045  0.9044  0.9890  NA      0.9895  0.9896  0.9042
NC_018497       0.9039  0.9041  0.9041  0.9914  0.9895  NA      0.9893  0.9038
NC_018498       0.9045  0.9047  0.9046  0.9898  0.9896  0.9893  NA      0.9044
NC_020076       0.9981  0.9954  0.9952  0.9046  0.9042  0.9038  0.9044  NA
```

and TETRA output:

```
# calculate_ani.py Sat Mar  1 19:04:45 2014
# TETRA
	NC_000912	NC_016807	NC_017504	NC_018495	NC_018496	NC_018497	NC_018498	NC_020076
NC_000912	NA	0.9998	0.9999	0.7376	0.7358	0.7401	0.7363	1.0000
NC_016807	0.9998	NA	0.9999	0.7359	0.7342	0.7386	0.7347	0.9998
NC_017504	0.9999	0.9999	NA	0.7359	0.7342	0.7386	0.7346	0.9999
NC_018495	0.7376	0.7359	0.7359	NA	0.9994	0.9995	0.9992	0.7373
NC_018496	0.7358	0.7342	0.7342	0.9994	NA	0.9994	0.9994	0.7355
NC_018497	0.7401	0.7386	0.7386	0.9995	0.9994	NA	0.9994	0.7398
NC_018498	0.7363	0.7347	0.7346	0.9992	0.9994	0.9994	NA	0.7360
NC_020076	1.0000	0.9998	0.9999	0.7373	0.7355	0.7398	0.7360	NA
```