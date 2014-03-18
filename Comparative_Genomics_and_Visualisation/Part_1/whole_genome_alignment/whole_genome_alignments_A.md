# Whole Genome Alignments

## Overview

This activity is intended to guide you through simple examples of whole genome alignment using two alignment packages: **megaBLAST** and **MUMmer**. You should work with the alignment software at the command line, following the examples as indicated.

The purpose of this activity is

1. to demonstrate and give practice in use of **BLAST**, **megaBLAST**, and **MUMmer**/**NUCmer** at the command-line.
2. to demonstrate and give practice in **MUMMer**'s dotplot visual output.
3. to demonstrate and give practice in the use of Sanger's **ACT** package to visualise genome alignments.
4. to give practice in, and examples of, interpretation of whole-genome alignments

When you have completed this activity you should be able 

1. to use the appropriate choice of **BLAST**, **megaBLAST**, and **MUMmer**/**NUCmer** for a whole genome alignment.
2. to visualise and navigate around whole-genome alignments using **ACT**.
3. to carry out basic interpretation of whole-genome alignments.

# A. Alignment of two closely-related bacterial genomes

Firstly, we'll compare some closely-related complete bacterial genomes. These are small, so all alignments should be fast. Being closely-related, we should see little divergence between genomes. Also, being completely sequenced bacterial genomes, there should not be any awkward issues with fragmentation, rearrangements or scaffolding.

To begin with, we'll compare the *E. coli* O157:H7 Sakai chromosome (`NC_002695.fna`) with that of *E. coli* CFT 073 (`NC_004431.fna`), using **megaBLAST** and **MUMmer**. These two bacteria are different isolates of the same bacterial species. Although they have distinguishable phenotypes, it is reasonable to expect them to be closely related at the sequence level.

### Check we're in the correct directory

We should check that we're working in the current directory, using the `pwd` command. This should give a directory ending in `Teaching/2014-03-07_University_of_Dundee` (the precise location will depend on where you cloned the repository to):

```
$ pwd
[...]/Teaching/2014-03-07_University_of_Dundee/whole_genome_alignment
```

If this is not the case, change to that directory now, or adapt the commands below, accordingly. 

We'll be collecting program output in the `wga_output` subdirectory. If this does not already exist, you will need to make it using the command

```
$ mkdir wga_output
```

## 1. Run default **BLASTN**

As a starting point for comparison, run a default (legacy-style) **BLASTN** search on the two *E. coli* genomes:

```
$ time blastn -query ../data/NC_002695.fna -subject ../data/NC_004431.fna -outfmt 6 -out wga_output/E_coli_blastn.tab -task blastn
```

**NOTE:** The `time` command reports how long it takes to run the command-line instructions. It is useful for getting a rough measurement of the time taken to run your commands, so you can estimate efficiency. In this case, the reported `real` time - also known as *wall clock* time - indicates the time taken for all operations, including CPU waiting, and it does not accurately represent CPU usage time. If you need to profile algorithms directly for efficiency, it is usually better to use a dedicated profiling tool or module.

You can inspect the resulting output file with the `head` and `wc` commands (see below).

```
$ head wga_output/E_coli_blastn.tab 
$ wc wga_output/E_coli_blastn.tab 
```

## 2. Run default **megaBLAST**

At the command line, this requires only a small modification of the usual **BLASTN** syntax to run:

```
$ time blastn -query ../data/NC_002695.fna -subject ../data/NC_004431.fna -outfmt 6 -out wga_output/E_coli_megablast.tab -task megablast
```

Again, we can the `head` and `wc` commands can be used to get a quick overview of the data.

```
$ head wga_output/E_coli_megablast.tab 
$ wc wga_output/E_coli_megablast.tab 
```

## 3. Run discontinuous megaBLAST

This is a similarly small modification to the default **BLASTN** syntax as before:

```
$ time blastn -query ../data/NC_002695.fna -subject ../data/NC_004431.fna -outfmt 6 -out wga_output/E_coli_dc-megablast.tab -task dc-megablast
```

**ACTIVITY 1:** Were there differences in the time taken, or the final output alignments, between the three **BLAST** algorithms? On the basis of these results, which **BLAST** method would you be likely choose for closely-related bacterial genome alignments, and why?

*HINT*: You can use the `diff` command (see `man diff` for details) to list differences between output files

## 4. Run **MUMmer**

**MUMmer** works quite differently to any **BLAST** alignment, at the command-line, so you may find it useful to inspect the documentation ([link](http://mummer.sourceforge.net/manual/)) from time to time.

**NOTE:** **MUMmer** is both the name of an alignment program (`mummer`), as well as the suite of alignment programs.

We'll start with `mummer`, which is useful to get a global view of a whole genome alignment, but is not as well-suited as `nucmer` (see below) for more detailed analysis

Run the command:

```
$ time mummer -mum -b -c ../data/NC_002695.fna ../data/NC_004431.fna > wga_output/E_coli_mummer.mums
```

Here, the options mean:

* `-b`: compare both forward and reverse strand sequences
* `-c`: report reverse complement sequence match locations relative to the forward strand

**NOTE:** the *first* named sequence is the *reference* sequence. This is important to know, especially if you are making a choice of query and reference on the basis of sequence size, for speed (which is recommended under some circumstances - see **MUMmer** documentation for details).

A quick idea of the output format of **MUMmer**, and the number of alignments, can be gained with the `wc` and `head` commands, as before, but it should be clear that the **MUMmer** output format differs quite greatly from **BLAST** output, and would need to be handled differently in any script you write.

```
$ head wga_output/E_coli_mummer.mums 
$ wc wga_output/E_coli_mummer.mums 
```

The program `mummer` is useful for getting an overview of a global alignment, and the `mummerplot` command is provided so that you can visualise its output. `mummerplot` provides script output suitable for the plotting package [`GNUplot`](http://www.gnuplot.info/). `gnuplot` allows you to plot figures directly from a command-line, or a script, and also has an interactive environment for plotting. `mummerplot` will run `gnuplot` automatically to create a figure for you, so long as that package is installed.

To visualise the `mummer` output then, run:

```
$ mummerplot --png --prefix=wga_output/E_coli_mummer wga_output/E_coli_mummer.mums
```

This will render a PNG format figure of your alignment, like that displayed below.

![MUMmer pairwise graphical output](./images/E_coli_mummer.png?raw=True =300x)

In the figure above, you should see `mummer`'s output displayed with the reference genome along the X-axis, and the query genome on the Y-axis. You can view the figure you created by issuing `eog wga_output/E_coli_mummer.png` at the command-line.

Red points/lines indicate an alignment in forward orientation, and blue points/lines an alignment in reverse orientation (with respect to the reference sequence).

* The dominant feature of the alignment is the broken red line from bottom-left to top-right. This indicates a substantial full-length alignment between the two genomes.
* Breaks in the diagonal that look like 'jumps' in the X- or Y-axis represent insertions/deletions in either the query or reference sequence
* Diagonal runs of alignments that are off the main diagonal indicate repeated similar sequence, or sequence rearrangements.
* It is usual with `mummer` output to see many individual alignments as 'noise' off the main diagonal. 

It is possible to run `mummer` with more than one query file. For example, we could introduce a third *E. coli* chromosome into the alignment as follows:

```
$ time mummer -mum -b -c ../data/NC_002695.fna ../data/NC_004431.fna ../data/NC_010468.fna > wga_output/E_coli_mummer_threeway.mums
$ mummerplot --png --prefix=wga_output/E_coli_mummer_threeway wga_output/E_coli_mummer_threeway.mums
```

and the graphical output changes (see below):

![MUMmer three-way graphical output](./images/E_coli_mummer_threeway.png?raw=true =300x)

If you viewed the threeway alignment image above, you would have seen that more blue entries were introduced, notably two blue diagonal lines, corresponding to most of the reference genome sequence. This results from the overlay of two independent alignments, on the same plot. This is not easy to interpret, as there is no way to determine visually which points and lines derive from which alignment.

**ACTIVITY 2:** Using `mummer`, align the `../data/NC_010468.fna` sequence against the same reference as above, and generate a plot of the alignment using `mummerplot`. Interpret the image. Both sequences represent completely sequenced chromosomes of record, stored in the public repository at NCBI. What does this alignment suggest about the agreement between genome sequences of the same species in public repositories? 

## 5. Run NUCmer

`nucmer` is actually a script that implements a pipeline of analysis, using `mummer` as its aligner. This pipeline  allows multiple reference and multiple query sequences to be aligned in a many-to-many alignment. That makes `nucmer` particularly useful for scaffolding assembly contigs or other fragmented sequences together against a reference genome. It is also useful for the kinds of three (and more)-way alignments that lead to confusing figures if `mummer` is used on its own.

Applying `nucmer` to genome alignment is more involved than running a simple `mummer` alignment. There are several stages to the analysis, if we want to produce visualisations or readily-interpretable output:

* Use `nucmer` to identify matching regions, producing a `.delta` file.
* Use `delta-filter` to remove weaker matches from the `.delta` file output, and produce a `.filter` file.
* Use `show-coords` and/or `show-aligns` to generate human-readable output from the `.delta` and/or `.filter` files.
* Use `mummerplot` to visualise the filtered matches in the `.delta` or `.filter` files.

A typical analysis would involve running the following commands (which you should do at the command-line):

```
$ time nucmer --maxgap=500 --mincluster=100 --prefix=wga_output/E_coli_nucmer ../data/NC_002695.fna ../data/NC_004431.fna

$ show-coords -r wga_output/E_coli_nucmer.delta > wga_output/E_coli_nucmer.coords

$ show-aligns wga_output/E_coli_nucmer.delta "gi|15829254|ref|NC_002695.1|" "gi|26245917|ref|NC_004431.1|" > wga_output/E_coli_nucmer.aligns

$ delta-filter -q -r wga_output/E_coli_nucmer.delta > wga_output/E_coli_nucmer.filter

$ mummerplot --png wga_output/E_coli_nucmer.filter -R ../data/NC_002695.fna -Q ../data/NC_004431.fna --prefix=wga_output/E_coli_nucmer
```

Note that the `show-aligns` command requires you to know and use the FASTA sequence IDs of your query and reference sequences.

As before, a PNG format figure of your alignment will be rendered:

![NUCmer pairwise graphical output](./images/E_coli_nucmer.png?raw=true =300x)

You will notice several things about this image in comparison to the `mummer` alignment output:

* It is much less "noisy" - most of the off-diagonal alignments have been filtered out.
* The reference genome is named explicitly (it was referred to as "REF" in the `mummer` output).
* The graph axes limits correspond to the start and end points of the genomes, not the next round number.

`nucmer` output, once filtered using `delta-filter` is much more useful than `mummer`'s for whole genome alignment. It is obvious from the graph that there is good one-to-one alignment across nearly the complete length of both genomes. The vertical/horizontal arrangement of points around the 1.5Mbp mark is suggestive of repeat or other similar regions that may be worth further investigation. There do not appear to be large-scale genome rearrangements.

The graph is easy to interpret but however, a quick look at the output files with `head` or in an editor, will quickly show that the raw `nucmer` output is quite complex. The `show-coords` and `show-aligns` commands used above give more readily-interpretable human-readable output, in a familiar style.

```
# Raw NUCmer output
$ head output/E_coli_nucmer.delta

# Filtered NUCmer output (used for plotting)
$ head output/E_coli_nucmer.filter

# show-coords output
$ head output/E_coli_nucmer.coords

# show-aligns output
$ head -n 20 output/E_coli_nucmer.aligns
```

**ACTIVITY 3:** Generate a visualisation of the `E_coli_nucmer.delta` output file, using `mummerplot`. Compare this with the plots generated from the `E_coli_mummer.mums` and `E_coli_nucmer.filter` output, and comment on the differences.

## 6. Visualising BLAST comparisons using ACT

**ACT** is an extension of the **Artemis** genome browser ([homepage](http://www.sanger.ac.uk/resources/software/artemis/#download)) designed to visualise genome comparisons. The software is free, and was developed at The Wellcome Trust Sanger Institute, and is widely used for genome annotation.

* **Artemis** manual: ([link](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.pdf))
* **ACT** manual: ([link](ftp://ftp.sanger.ac.uk/pub/resources/software/act/act.pdf))

**ACT** and **Artemis** will take sequence data in `FASTA` format, and annotations in `GFF` format. They will also read combined sequence and annotation data in `GenBank` format. Genome comparison data can be read in two essentially equivalent formats: `.crunch` files, and **BLAST** tabular output.

You will use ACT to visualise and compare the outputs of **megaBLAST** and **MUMmer** comparisons, above.

To start **ACT** issue

```
$ act
```

at the command-line.

### 6.1 **BLAST** vs. **megaBLAST**

You should have found above that the outputs of `BLASTN` and `megaBLAST` produced different output. To become more familiar with **ACT**'s operation, you will go through the steps of loading in the two alignments.

![ACT splash screen](./images/act_fig1.png?raw=true =200x)

Use the `File -> Open` menu option to obtain the file selection dialogue box, then click on the `more files ...` button to obtain an option to enter a third sequence:

![ACT file selection dialogue](./images/act_fig2.png?raw=true =200x)

![ACT file selection dialogue (expanded)](./images/act_fig3.png?raw=true =200x)

Use the `Choose ...` buttons to select the `NC_004431.fna` FASTA file (in the `../data` subdirectory) for sequence files 1 and 3, and the `NC_002695.fna` FASTA file as sequence file 2:

![ACT file selection dialogue (sequences selected)](./images/act_fig4.png?raw=true =200x)

Then use the `Choose ...` buttons to select the `wga_output/E_coli_blastn.tab` and `wga_output/E_coli_megablast.tab` files as comparison files 1 and 2, respectively:

![ACT file selection dialogue (all files selected)](./images/act_fig5.png?raw=true =200x)

Then click on the `Apply` button. You will see a notification window to say that the **dc-megaBLAST** hits were flipped to match sequence orientation, and the main window showing a portion of the genome alignment:

![ACT initial alignment view](./images/act_fig6.png?raw=true =400x)

In this view, the top pair of grey bars, and the bottom pair of grey bars each indicate the forward and reverse strands of the `NC_004431` (*E. coli* CFT 073) sequence. The centre pair of grey bars indicate the forward and reverse strands of the `NC_002695` (*E. coli* O157:H7) sequence. If there were functional annotations loaded for these sequences, they would be placed on these bars. However, we are looking only at the nucleotide sequences.

The red bars connecting each of the genomes indicate the locations of the sequence alignments you calculated. The top set are those from the **BLAST** comparison, and the bottom set are those from the **megaBLAST** comparison. As you can see in the figure above, they are very similar but not identical.

Use the scrollbars on the right-hand side (next to the genome tracks) to zoom out, and show the complete alignment:

*Hint:* this can be quite fiddly on a small monitor - it may be helpful to maximise the **ACT** window to full-screen.

*Hint:* double-clicking on any of the red or blue alignment links will turn them yellow, and align the two connected genomes at that location, as you can see below. Selecting a link in this way also reports useful information about the alignment, on the left-hand side.

![ACT initial alignment, zoomed out](./images/act_fig7.png?raw=true =400x)

As you can see, both **BLAST** and **megaBLAST** report many matches and there is a confusing mass of lines criss-crossing between many regions of the genomes.

To simplify the view, low-scoring alignments can be filtered, using the scroll-bars to the right of the screen. Applying the maximum filtering level makes the genome similarities clearer:

![ACT alignment, filtered](./images/act_fig8.png?raw=true =400x)

Now we can see that there are many large regions of similarity, indicated by the red connections between genomes.

Red connections represent matches that run in the same direction on the two genomes being compared, and blue connections indicate alignments that run in opposite directions on each genome. (though note that **ACT** allows you to flip the orientation of any genome, interactively).

We can see from the 'wall' of red connections quite quickly that these two genomes are very similar across most of their lengths, and that the larger alignments are all in the same order.

However, there are many smaller alignments - mostly in blue - that appear to radiate from a single point on the `NC_004431` genome, to many points on the `NC_002695` genome, and *vice versa*. These are indicative of repeated regions of sequence similarity, which may suggest (for bacteria) phage integration, or other repetitive elements.

This kind of view also draws attention to the regions between alignments - where the two genomes differ. There are 'wedge-like' gaps - small in one genome, large in the other - that suggest an insertion or deletion in one or other genome. There are also gaps that are of approximately equal size in each genome, which may indicate sequence divergence at that location, or a common insertion site.

#### Genomic insertion

One way to identify a possible genomic insertion, such as a pathogenicity island, is to look at the nucleotide (GC) content in that region. It is often the case for bacteria that an inserted sequence may come from some other organism (*via* lateral/horizontal gene transfer: HGT/LGT) with a different balance of nucleotide usage. The inserted regions thus stand out when nucleotide use is plotted.

**ACT** allows us to show these graphs in conjunction with the genome comparison data.

Using **ACT**, zoom in to the region around 3,276,000bp in the `NC_004431` sequence, and 3,685,500bp in the `NC_002695` sequence

![ACT focused insertion](./images/act_fig9.png?raw=true =400x)

Use the `Graph -> NC_002695.fna -> GC Content (%)` option from the menu bar to render a %GC content graph for the central genome: `NC_002695` *E. coli* O157:H7 Sakai:

![ACT graph choices](./images/act_fig10.png?raw=true =200x)

The resulting graph shows us that, where there is a section of genome sequence present in `NC_002695`, but not `NC_004431` (3709877-3737122 in `NC_002695` co-ordinates), there is a corresponding variation in %GC content. This is potentially indicative of an insertion event that resulted in a genomic island.

In fact, this region is genomic island GI28, as described in [Roos & van Passel (2011)](http://www.biomedcentral.com/1471-2164/12/427).

**NOTE:** You can use the slider to the right of the graph to vary the window size over which the %GC content statistic is calculated.

![ACT showing graph](./images/act_fig11.png?raw=true =400x)

### 6.2 **megaBLAST** vs **MUMmer**

As noted above **MUMmer** (and `nucmer`) output can be complex, but it is possible to generate human-readable output from this complex output. However, the files produced by e.g. `show-coords` are not directly readable by **ACT**. This is a common issue with genome comparisons, as applications produce many different types of data, and there are few standards in the field.

The `.crunch` tabular file accepted by **ACT** is simple, and can be generated from the output of **MUMmer**'s `show-coords` package. The script `nucmer_to_crunch.py` in the `scripts` subdirectory can do this for you. To use it on the `nucmer` output you generated above, run the following at the command-line:

```
$ python ../scripts/nucmer_to_crunch.py -i wga_output/E_coli_nucmer.coords -o wga_output/E_coli_nucmer.crunch -v
```

and then proceed as above, to load in the same genome sequences as before, but replacing the **BLAST** alignment data with that from the `nucmer` alignment. Your file selection dialogue should resemble that below:

![ACT file dialogue](./images/act_fig12.png?raw=true =200x)

**ACTIVITY 4:** Click on `Apply` and compare the visualisations of the **megaBLAST** and `nucmer` whole genome comparisons for `NC_002695` and `NC_004431`, and comment on them. What are the differences between the two outputs? Do the differences appear to be significant? Is there any reason to choose one alignment tool over the other?
