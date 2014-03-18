# README.md - mcl_orthologues

## Overview

This activity is intended to

1. demonstrate use of **MCL** (Markov Clustering) to cluster protein sequences into groups of putative orthologues
2. give an example of comparison of RBBH and MCL-based orthologue prediction, on the same dataset

When you have completed this activity you will

1. be able to use **MCL** to cluster protein sequences on the basis of **BLAST** comparisons
2. be able to make informed comparisons and choices between alternative orthologue prediction methods

Note that, although **OrthoMCL** <http://orthomcl.org/orthomcl/> is preferred over plain **MCL** for orthologue finding, the software installation required is prohibitively complex for a simple training example. 

## Prerequisites

These activities were written and tested using the following software versions, though others may also work:

* **MCL** 12-135 <http://micans.org/mcl/>
* **BLAST** 2.2.29+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>

# 1. Perform MCL Clustering
## All-vs-All BLASTP

The first step of **MCL** clustering analysis is to generate output describing the sequence similarities. This is done using **BLAST**.

As we want an all-against-all comparison, we concatenate all protein sequences from the `find_rbbh` activity into a single FASTA file, and build a database from that file. Then we use **BLASTP** to query those protein sequences against themselves.

In order to restrict the number of false positive associations between sequences, at the risk of failing to identify divergent pairs of sequences that are truly homologous, we're using an E-value cutoff of 1e-30.

This process is described in the script file `run_BLAST.sh`, in this directory. The search takes a while so, to save time, it has been run for you already, and the output should already be in the `data` directory.

## Create an `.abc` File

The `.abc` format contains three columns: `sequenceID1`, `sequenceID2`, `Evalue`. This could be obtained directly from the **BLASTP** search with an appropriate argument to `-outfmt`, but instead we keep the full **BLASTP** output for reference, and use the `cut` command to extract the appropriate columns:

```
$ cut -f 1,2,11 data/all-vs-all.tab > data/all-vs-all.abc
```

## Create Initial Network

Use the `mcxload` command to create the network file `seq.mci` and a corresponding dictionary file `seq.tab`, that will connect sequence identifiers with the values used internally by **MCL** in the network description.

```
$ mcxload -abc data/all-vs-all.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o data/seq.mci -write-tab data/seq.tab
```

This command is converting the **BLAST** output to a long list of which sequences make BLAST matches to which other sequences, and with what E-value, converting that E-value by making it -log10(E-value). Space is also being saved by representing every sequence as an integer, and only storing the names once, in `seq.tab`.

```
$ head data/seq.*
==> data/seq.mci <==
(mclheader
mcltype matrix
dimensions 12336x12336
)
(mclmatrix
begin
0     0:104 1:99.09691 2:64.09691 3:23 4:19.69897 5:19 6:14.52288
        7:14.52288 8:11.69897 9:11.52288 10:10.22185 11:9.045757 12:1.21467
        13:0.5086383 14:0.4685211 15:0.455932 16:0.455932 17:-0.04139268
        18:-0.2304489 19:-0.30103 20:-0.5185139 21:-0.6812412 22:-0.7242759

==> data/seq.tab <==
0	gi|50118966|ref|YP_048133.1|
1	gi|261823728|ref|YP_003261834.1|
2	gi|188532166|ref|YP_001905963.1|
3	gi|188534847|ref|YP_001908644.1|
4	gi|261822604|ref|YP_003260710.1|
5	gi|50119966|ref|YP_049133.1|
6	gi|188532464|ref|YP_001906261.1|
7	gi|50121000|ref|YP_050167.1|
8	gi|261822764|ref|YP_003260870.1|
9	gi|261820093|ref|YP_003258199.1|
```

Technically, our **BLAST** data is now a graph, where each the first column of each line in `seq.mci` represents a node in the graph (identified by an integer), and the second column describes the edges and their corresponding weights, in `node:weight` pairs.

Some information is being lost, here. A **BLAST** match between two sequences A and B may not produce the same E-value when run in each direction. The `--stream-mirror` argument is assigning the "best" E-value as the edge weight between each pair of sequences. 

## Cluster the Network

We will be using `mcl` to perform the clustering operation. This command takes a single argument to control the inflation value used in the algorithm. We will use one setting for example purposes though, in real use, you would want to ensure that clustering is robust for your chosen inflation value.

Here, we use an inflation value of 6, and generate the output file `out.seq.mci.I60`:

```
$ mcl data/seq.mci -I 6 -use-tab data/seq.tab -o data/out.seq.mci.I60
```

and you can inspect the output in a text editor or with 

```
$ head data/out.seq.mci.I60
```

In this output, each line represents a cluster of sequences, where the sequence IDs of members of the cluster are all given in tab-separated plain text format.

### Process the cluster data

Please now use the `mcl_orthologues.ipynb` iPython notebook to process the **MCL** cluster output into `.crunch` output, so you can use **ACT** to compare clustering and RBBH output.

Start the notebook by issuing:

```
$ ipython notebook --pylab inline
```

at the command-line, in the `mcl_orthologues` directory, and clicking on the notebook filename in your browser window.

# 2. Visualise MCL clustering output with ACT

After using the `mcl_orthologues.ipynb` notebook, you should have generated the cluster output `.crunch` files, and now the `data` subdirectory should look like this:

```
$ tree data
data
├── NC_004547.gbk
├── NC_004547_vs_NC_004547.crunch
├── NC_010694.gbk
├── NC_010694_vs_NC_004547.crunch
├── NC_010694_vs_NC_010694.crunch
├── NC_010694_vs_NC_013421.crunch
├── NC_013421.gbk
├── NC_013421_vs_NC_004547.crunch
├── NC_013421_vs_NC_013421.crunch
├── all-vs-all.abc
├── all-vs-all.tab
├── out.seq.mci.I60
├── proteins.faa
├── proteins.faa.phr
├── proteins.faa.pin
├── proteins.faa.psq
├── seq.mci
└── seq.tab
```

Now, using **ACT**, load the three GenBank files `NC_004547.gbk`, `NC_013421.gbk`, `NC_010694.gbk`, and the corresponding `.crunch` output:

**NOTE:** Be sure to use the files in the order shown: `NC_010694.gbk`, `NC_013421.gbk` and `NC_004547.gbk`, with the appropriate `.crunch` files.

![ACT file selection](images/act1.png?raw=True =300x)

This should give you output that resembles the image below. Note that by loading in annotated GenBank files, rather than plain sequences, we see green bars indicating genome features, such as genes:

![ACT comparison](images/act2.png?raw=True =400x)

On zooming out, you should be able to see that the two *Pectobacterium* genomes (top and middle, in the image below) are, on the whole, more similar to each other than to the *Erwinia* genome, which is consistent with their taxoomic classification. There is also a large inversion in the *Pectobacterium* comparison.

Both comparisons show some evidence of synteny, but the degree of rearrangement is much greater in the *Pectobacterium*:*Erwinia* comparison, than between the two *Pectobacterium* genomes.

![ACT zoomed out](images/act3.png?raw=True =400x)

# 3. Visual comparison of MCL and RBBH orthologue predictions

If you have completed the `find_rbbh` activity, you should by now have two alternative sets of orthologue predictions, which you can compare using **ACT**.

Select the appropriate files for comparison (as shown below), and visualise the putative RBBH predictions:

![ACT file selection](images/act4.png?raw=True =200x)

Once you have zoomed out a little, and aligned the genomes, you should see an image like this:

![ACT RBBH comparisons](images/act5.png?raw=True =400x)

The MCL predictions are at the top, and the RBBH predictions at the bottom. You should see that the overall agreement between the methods, especially regarding genome synteny and the rearrangment events, is good. However, the number of relationships identified by the RBBH method is fewer than that identified with MCL. On the whole, the MCL comparisons predictions seem to be "noisier", and contain more predictions where multiple sequences from the same genome contribute to a cluster.

**ACTIVITY 1:** Are multiple-member clusters consistent with the concept of orthology? How does this affect interpretation of the MCL output? What advantages and disadvantages does the MCL clustering approach have with respect to RBBH for identifying protein family members and evolutionary relationships, and how are they evident in these comparisons? What else could you do to "tidy up" or improve the MCL orthologue-finding approach used in this activity?

*HINT:* You may find it useful to zoom into the comparison image, and identify pairs of proteins that are in an MCL cluster, but do not make RBBH, and inspect the one-way BLAST output to make a judgement on whether they are likely to be related.