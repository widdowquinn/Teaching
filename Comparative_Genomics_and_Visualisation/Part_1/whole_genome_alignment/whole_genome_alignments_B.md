# Whole Genome Alignments

## Overview

This activity is intended to guide you through simple examples of draft genome alignment using two alignment packages: **Mauve** and **MUMmer**. You should work with **MUMmer**  at the command line, and the **Mauve** GUI, following the examples as indicated.

The purpose of this activity is

1. to demonstrate and give practice in use of the **Mauve** GUI, and **MUMmer**/**NUCmer** at the command-line.
2. to demonstrate and give practice in **MUMMer**'s dotplot visual output.
3. to demonstrate and give practice in **Mauve**'s visual output.
4. to demonstrate and give practice in the use of Sanger's **ACT** package to visualise draft genome alignments.

When you have completed this activity you should be able 

1. to use **Mauve** and **MUMmer**/**NUCmer** to produce a draft genome alignment.
2. to visualise and navigate around draft genome alignments using **ACT**.

# B. Alignment of a draft and reference bacterial genome

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

## 1. Run **MUMmer**

Aligning a draft genome to a reference genome in **MUMmer** is almost identical to aligning two complete genomes. Our draft assembly takes the form of 41 unordered contig sequences from a *Dickeya solani* isolate, in FASTA format in the file `MK10_draft.fasta`. We will align this against the reference *Dickeya dadantii* genome sequence in `NC_014500.fna`.

To align and visualise the alignment of the draft assembly contigs in `MK10_draft.fasta` to the reference in `NC_014500.fna`, run the following commands:

```
$ nucmer --prefix=wga_output/MK10_vs_NC_014500 ../data/NC_014500.fna ../data/MK10_draft.fasta

$ mummerplot --png wga_output/MK10_vs_NC_014500.delta -R ../data/NC_014500.fna -Q ../data/MK10_draft.fasta --filter --layout --prefix=wga_output/MK10_vs_NC_014500
```

![**NUCmer** alignment of draft genome](./images/MK10_vs_NC_014500.png?raw=True =400x)

The resulting output shows the draft assembly contigs on the Y-axis, ordered to produce an optimal alignment against the reference. The quality of the overall alignment is clear from the strong red diagonal indicating aligned regions. **NUCmer**'s optimal ordering of contigs can be read from the Y-axis. 

A human-readable co-ordinate table can be produced as before, with

```
$ delta-filter -q wga_output/MK10_vs_NC_014500.delta > wga_output/MK10_vs_NC_014500.filter

$ show-coords -rcl wga_output/MK10_vs_NC_014500.filter > wga_output/MK10_vs_NC_014500_filtered.coords
```

and the very last column of `MK10_vs_NC_014500_filtered.coords` gives the ID of each aligned contig. This output can also be processed to recover **NUCmer**'s optimal ordering of contigs.


## 2. Run **Mauve**

**Mauve** is an alignment package produced by the Genome Evolution Laboratory at the University of Wisconsin-Madison, and can be obtained at [http://gel.ahabs.wisc.edu/mauve/](http://gel.ahabs.wisc.edu/mauve/). There are command-line and GUI versions of **Mauve**, but this activity will only use the GUI version.

The algorithm that **Mauve** uses is equally applicable to complete genome sequences, and can be run in single, and *progressive* (i.e. iterated) modes. In this activity, you will use the progressive mode to align the draft `MK10_draft.fasta` assembly against the `NC_014500.fna` reference.

Start **Mauve** at the command-line with the command:

```
$ Mauve
```

You should be greeted by the splash screen, and then a window as below:

![Mauve start window](images/mauve1.png?raw=True =300x)

### Progressive alignment

From the `File` menu item, select `Align with progressiveMauve`. You should then see a file selection dialogue box:

![Mauve `File` options](images/mauve2.png?raw=True =200x)

![Mauve file selection dialogue](images/mauve3.png?raw=True =200x)

There are tabs available for you to change parameters related to the alignment - you should ignore these for now:

![Mauve progressive alignment parameter options](images/mauve4.png?raw=True =200x)

![Mauve progressive alignment scoring options](images/mauve5.png?raw=True =200x)

Once the appropriate files have been selected, and an output file location has been chosen, click on the `Align…` button:

![Mauve progressive alignment file choices](images/mauve6.png?raw=True =200x)

A window will pop up, informing you of progress. This can be saved as a log of the alignment.

![Mauve progressive alignment progress window](images/mauve7.png?raw=True =400x)

When it is complete, **Mauve** will present you with a visualisation of the draft genome alignment:

![Mauve alignment](images/mauve8.png?raw=True =400x)

This view shows the (complete) `NC_014500` *Dickeya dadantii* genome in the top row, and the draft *Dickeya solani* genome in the bottom row. The major LCBs are indicated in colour blocks (there are 18 indicated, with minimum weight 2260).

On the bottom row, the draft genome contigs are indicated, with contig boundaries marked as red lines (i.e. the regions between consecutive red lines indicate individual contigs). By dragging the mouse/cursor along this bottom row, you should see the currently active contig IDs change in the status bar at the bottom (you will note that all the contigs remain in numerical order, and have not been moved). Blocks above the line are in the forward strand with respect to the input data, and those below the line are reversed with respect to the input.

The thin lines linking LCBs are guides to the eye, so that rearrangements can be 

This alignment produces four output files: 

```
$ tree wga_output/Alignment/
wga_output/Alignment/
├── MK10_vs_NC_014500.mauve
├── MK10_vs_NC_014500.mauve.backbone
├── MK10_vs_NC_014500.mauve.bbcols
└── MK10_vs_NC_014500.mauve.guide_tree
```

Other than the `.guide_tree` file, these outputs are all in **Mauve**-specific formats, and are not immediately helpful.

The important thing to note is that, while **Mauve** has identified LCBs and it is clear that there are some contigs that could be moved to improve the overall alignment, it has not actually performed these movements.

### Mauve Contig Mover (MCM)

The **Mauve Contig Mover** is more helpful, if we want an output that we can process further, more easily. To use it, select `Tools -> Move Contigs` from the menu:

![MCM menu](images/mcm1.png?raw=True =200x)

and select a directory in which to keep the resulting progressive Mauve output (I used `wga_output/Reordering`):

![MCM directory selection](images/mcm2.png?raw=True =200x)

An informative message will appear

![MCM information](images/mcm3.png?raw=True =400x)

and then you should choose input files in the dialogue box.

**NOTE:** This process requires that your **reference** sequence is named first, so choose `NC_014500.fna` as the first file.

![MCM file choice](images/mcm4.png?raw=True =200x)

Clicking on `start` will bring up another log window that will inform you about progress.

Unlike **Mauve** proper, **MCM** will place intermediate alignments and output in subdirectories `alignmentN` of your chosen directory. This allows you to inspect individual alignments, and to recover the final reordering of your draft contigs, conveniently. In this case the number of alignments required is 4, and you should take `alignment4` to be your final alignment.

This order of contigs on the bottom row is not necessarily the same order as originally given to progressive Mauve (indeed, this is what we want!). It should resemble the reordering as given by **MUMmer** above.

![MCM final alignment](images/mcm5.png?raw=True =400x)

Interpretation of the visualisation is the same as for progressive Mauve: coloured blocks represent LCBs, and contig boundaries are represented as red lines. The output files delivered by **MCM** are more useful, too:

```
$ tree wga_output/Reordering/alignment4/
wga_output/Reordering/alignment4/
├── MK10_draft.fasta
├── MK10_draft.fasta.sslist
├── MK10_draft_contigs.tab
├── NC_014500.fna
├── NC_014500.fna.sslist
├── alignment4
├── alignment4.backbone
├── alignment4.bbcols
└── alignment4.guide_tree
```

In practical terms, **MCM** has the advantage of producing output with explicit ordering information, in the `MK10_draft_contigs.tab` file, and the final `MK10_draft.fasta` file has the input contigs reordered into the best alignment, as determined by **Mauve**:

```
$ grep '>' wga_output/Reordering/alignment4/MK10_draft.fasta | head
>1
>2
>33
>3
>4
>5
>6
>7
>8
>9
```

**Mauve** is a very useful package and, as a visualisation tool, will let you resolve LCB boundaries interactively, down to single-base level:

![MCM boundary](images/mcm6.png?raw=True =400x)


## 3. Visualising reordered (or not) fragments and alignments in **ACT**

If you attempt to align a draft genome in multiple fragments, reordered or not, and visualise it in **ACT** without any further process, you will probably be disappointed with the results. For example, running the alignment:

```
$ blastn -query wga_output/Reordering/alignment4/MK10_draft.fasta -subject ../data/NC_014500.fna -outfmt 6 -out wga_output/disappointing.crunch -task megablast
```

and attempting to visualise the results will give you the following impossible to interpret image:

![ACT disappointment](images/disappointing.png?raw=True =400x)

The issue here is that although **ACT** runs together the different fragments in the lower genome (indicated by the brown markers), and the `.crunch` file contains enough information for an alignment viewer to render the alignment correctly by offsetting alignment co-ordinates appropriately, **ACT** does not currently do this. 

One solution is to stitch the contigs together, and we can do this using the `stitch_six_frame_stops.py` script in the `scripts` directory. This will join the contigs in the order presented in the input FASTA file, connecting them with the sequence

```
NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN
```

which contains start and stop codons in all frames. You can do this by running:

```
$ ../scripts/stitch_six_frame_stops.py -i wga_output/Reordering/alignment4/MK10_draft.fasta -o wga_output/MK10_draft_stitched.fasta --id=MK10_stitched
```

The script also generates a corresponding `.gff3` file for convenience, which details the locations of the original contigs on the output stitched sequence.

Now we can run the alignment as:

```
$ blastn -query wga_output/MK10_draft_stitched.fasta -subject ../data/NC_014500.fna -outfmt 6 -out wga_output/MK10_vs_NC_014500.crunch -task megablast
```

and the resulting **ACT** visualisation is a bit more sensible:

![ACT success](images/success.png?raw=True =400x)