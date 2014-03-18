# README.md - predict_CDS

**NOTE:** Due to time constraints, this is an optional activity. The results of completing this activity should already be present, so that visualisation and analysis can proceed in later activities.

## Overview

This activity is intended to 

1. demonstrate *ab initio* gene-calling in prokaryotes
2. demonstrate visualisation of CDS calls using **Artemis**

When you have completed this activity, you should be able to

1. use the **Prodigal** genecaller to generate a set of genecalls on a prokaryotic genome
2. use **Artemis** to visualise and inspect the resulting genecalls

## Prerequisites

This activity was written and tested using the following software, though other versions may also be suitable:

**Prodigal** v2.50 <http://prodigal.ornl.gov/>

# 1. Predicting Prokaryotic CDS

**Prodigal** is a very fast prokaryotic gene-calling tool. It is used in a number of bacterial/archaeal annotation pipelines (e.g. **PROKKA** <http://www.vicbioinformatics.com/software.prokka.shtml>), and has a number of advantages, including that it can run in a single step (it doesn't *necessarily* need to be trained to perform well), and it is not strongly affected by the genome's GC content.

In our work, we found **Prodigal** to be the best of the bacterial gene-callers we tried, and it typically has near-100% 3' accuracy, and around 95% 5' accuracy on the genes it calls, using benchmark sets. On our data, we found that it had closer to 70-80% 5' accuracy, and 90-97% 3' accuracy. This was still better than the common alternatives we tried (including **RAST** <http://rast.nmpdr.org/>, **Glimmer** <http://ccb.jhu.edu/software/glimmer/index.shtml> and **GeneMark** <http://opal.biology.gatech.edu/>).

There is, however, [no free lunch](http://en.wikipedia.org/wiki/No_free_lunch_in_search_and_optimization), so you may wish to try and compare other genecallers.

In this activity, you will use **Prodigal** to annotate a complete bacterial genome sequence, and visualise the results.

**NOTE:** The **Prodigal** output is already present, to save some time in proceeding to later activities.

## Using Prodigal to Predict Prokaryotic CDS

Ensure that you are in the correct working directory:

```
$ pwd
[...]/Part_2/predict_CDS
```

The subdirectory `data` contains the file `genome.fna`, which describes the chromosomal genome sequence of a bacterium.

Use **Prodigal** to annotate this genome, by issuing:

```
$ prodigal -i data/genome.fna -a data/predicted_CDS.faa -d data/predicted_CDS.ffn -f gff -o data/predicted_CDS.gff
```

This will write the following files to the `data` directory:

* `predicted_CDS.ffn`: nucleotide sequences of predicted CDS in FASTA format
* `predicted_CDS.faa`: protein sequences of predicted CDS in FASTA format
* `predicted_CDS.gff`: the locations of the predicted CDS on the input genome, in [GFF](http://www.sequenceontology.org/gff3.shtml) format.

**ACTIVITY 1:** Inspect the **Prodigal** output. How many CDS are predicted? Does this seem reasonable?

## Visualising Predicted CDS

You will use **Artemis** (a companion program to **ACT**) to view **Prodigal**'s predictions on the bacterial genome.

Start **Artemis** from the command-line:

```
$ art &
```
You will see the landing window:

![Artemis landing window](images/art1.png?raw=True =300x)

Use the `File -> Open` menu option to get the file dialogue box, and use this to select the `genome.fna` file:

![Artemis menu option](images/art2.png?raw=True =200x)

![Artemis file selection dialogue](images/art3.png?raw=True =300x)

This will open the genome browsing window. This shows the genome sequence in a three reading frame view, six-frame conceptual translations, and stop codons (short vertical lines in each reading frame). Currently, no genome features are indicated, but potential ORFs are clearly visible as large gaps between successive stop codons in a single reading frame:

![Artemis browser view, no features present](images/art4.png?raw=True =400x)

Read in **Prodigal**'s CDS predictions, by using the `File -> Read An Entry...` menu option,  and choosing the `data/predicted_CDS.gff` file:

![Artemis entry menu option](images/art5.png?raw=True =200x)

![Artemis file dialogue](images/art6.png?raw=True =300x)

This will now show blue blocks in the genome browser window, and a table of feature locations in the lower window.

![Artemis browser view, with CDS features](images/art7.png?raw=True =400x)

You should be able to see that the predicted CDS have a good correspondence to the larger ORFs you should have been able to identify by eye.

**ACTIVITY 2:** Use **Artemis** to inspect the CDS predictions made by **Prodigal**.

*HINT:* Right-clicking on a feature will bring up a context menu that may be useful. You can use the editor (select a feature, and press `Ctrl-E`) to see and modify the feature annotation.