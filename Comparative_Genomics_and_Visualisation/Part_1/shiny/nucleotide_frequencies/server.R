# nucleotide_frequencies/ui.R
#
# This R script is part of an activity module designed to illustrate
# differences between nucleotide frequencies of bacterial genomes,
# for a University of Dundee bioinformatics module on 
# Comparative Genomics and Visualisation.
#
# The UI is displayed interactively in the user's browser, using Shiny
# (http://www.rstudio.com/shiny/), using 
# runApp("shiny/nucleotide_frequencies") in RStudio's interactive
# window. It displays some activity description text, and an input 
# selection dropdown allowing for a choice of four bacterial 
# species, and some collections of species.
#
# The UI displays boxplots of mono-, di-, tri- and tetranucleotide
# frequencies for the chosen set of bacterial sequences.

# We need Shiny to use Shiny.
library(shiny);
# GGPlot2 is used to render pretty boxplots, and gridExtra is needed
# to show four of them in an organised way
library(ggplot2);
library(gridExtra);
# Biostrings is a component of BioConductor, and makes handling 
# biological sequences easier. In particular, it provides our
# k-mer frequency counting function
library(Biostrings);
# We need to reshape (melt) our frequency data for plotting
library(reshape2);
# We use a hash (similar to a Python dictionary) to manage
# data in response to the selected dataset.
library(hash);

# Define server logic required to plot nucleotide frequencies
shinyServer(function(input, output) {
  
  # Expression to generate a plot of nucleotide frequencies.
  # This is wrapped in a call to renderPlot so that the output
  # type is a plot, and it is reactive to input choice
  output$freqPlot = renderPlot({
  
    # Define a hash of genome filenames, for each organism name.
    orghash = hash(keys=c("Escherichia coli", "Mycoplasma pneumoniae",
                          "Mycoplasma genitalium", "Mycobacterium tuberculosis",
                          "Mycoplasma spp.", "All", "Unknown"),
                   values = list(list("NC_000913.fna", "NC_002695.fna",
                                      "NC_004431.fna", "NC_010468.fna"), 
                                 list("NC_000912.fna", "NC_016807.fna", 
                                      "NC_017504.fna", "NC_020076.fna"),
                                 list("NC_018495.fna", "NC_018496.fna",
                                      "NC_018497.fna", "NC_018498.fna"),
                                 list("NC_016934.fna", "NC_017523.fna",
                                      "NC_022350.fna", "NC_000962.fna"),
                                 list("NC_000912.fna", "NC_016807.fna", 
                                      "NC_017504.fna", "NC_020076.fna",
                                      "NC_018495.fna", "NC_018496.fna",
                                      "NC_018497.fna", "NC_018498.fna"),
                                 list("NC_000913.fna", "NC_002695.fna",
                                      "NC_004431.fna", "NC_010468.fna", 
                                      "NC_000912.fna", "NC_016807.fna", 
                                      "NC_017504.fna", "NC_020076.fna",
                                      "NC_018495.fna", "NC_018496.fna",
                                      "NC_018497.fna", "NC_018498.fna",
                                      "NC_016934.fna", "NC_017523.fna",
                                      "NC_022350.fna", "NC_000962.fna"),
                                 list("NC_002695.fna")
                                 ));
        
    datadir = "../../data";  # Hardcoded wrt Shiny location for demo
    filelist = orghash[[input$orgname]];  # Choose files for dataset

    # Create dataframes to collect data
    mono <- data.frame(org=character(), ID=character(), GC=list());
    di <- data.frame(org=character(), ID=character(), GC=list());
    tri <- data.frame(org=character(), ID=character(), GC=list());
    tetra <- data.frame(org=character(), ID=character(), GC=list());
    
    # Loop over sequence files for the chosen organism, and calculate frequencies
    for(filename in filelist) {
      # Read files in one-by-one
      contigs <- readDNAStringSet(file.path(datadir, filename));
      
      # Get base composition data, and compile into a single dataframe for each composition type.
      monodata = data.frame(org=input$orgname, ID=names(contigs), length=width(contigs), 
                            f=(oligonucleotideFrequency(contigs, 1)/width(contigs)));
      mono <- rbind(mono, monodata);
      didata = data.frame(org=input$orgname, ID=names(contigs), length=width(contigs), 
                            f=(oligonucleotideFrequency(contigs, 2)/width(contigs)));
      di <- rbind(di, didata);
      if(input$all_plots == TRUE) {
        tridata = data.frame(org=input$orgname, ID=names(contigs), length=width(contigs), 
                            f=(oligonucleotideFrequency(contigs, 3)/width(contigs)));
        tri <- rbind(tri, tridata);
        tetradata = data.frame(org=input$orgname, ID=names(contigs), length=width(contigs), 
                            f=(oligonucleotideFrequency(contigs, 4)/width(contigs)));
        tetra <- rbind(tetra, tetradata);
      }
    }

    # Plot oligonucleotide frequencies as boxplots
    # Reorganise (melt) data for plotting
    m.mono = melt(mono, id.vars=c("org"), measure.vars=4:7); 
    m.di <- melt(di, id.vars=c("org"), measure.vars=4:19);
    # Generate boxplots
    p1 = ggplot(m.mono, aes(x=variable, y=value, fill=variable)) + geom_boxplot(outlier.colour="red") + 
      guides(fill=FALSE) + scale_fill_brewer() + 
      labs(title="Chromosome base frequencies", x='Base', y='% content') + facet_wrap(~org, ncol=5);
    p2 <- ggplot(m.di, aes(x=variable, y=value, fill=variable)) + geom_boxplot(outlier.colour="red") + 
      guides(fill=FALSE) + 
      labs(title="Chromosome dinucleotide frequencies", x='Dinucleotide', y='% content') + 
      facet_wrap(~org);

    if(input$all_plots == TRUE) {
      m.tri <- melt(tri, id.vars=c("org"), measure.vars=4:67);
      m.tetra <- melt(tetra, id.vars=c("org"), measure.vars=4:259);    
      p3 <- ggplot(m.tri, aes(x=variable, y=value, fill=variable)) + geom_boxplot(outlier.colour="red") + 
        guides(fill=FALSE) + 
        labs(title="Chromosome trinucleotide frequencies", x='Trinucleotide', y='% content') + 
        facet_wrap(~org, ncol=1) + theme(axis.text.x = element_blank());
      p4 <- ggplot(m.tetra, aes(x=variable, y=value, fill=variable)) + geom_boxplot(outlier.colour="red") + 
        guides(fill=FALSE) + 
        labs(title="Chromosome tetranucleotide frequencies", x='Tetranucleotide', y='% content') + 
        facet_wrap(~org, ncol=1) + theme(axis.text.x = element_blank());
    }
    
    # Return the gridded ggplots to the UI
    if(input$all_plots == TRUE) {
      print(grid.arrange(p1, p2, p3, p4, ncol=1));
    } else {
      print(grid.arrange(p1, p2, ncol=1));      
    }
    
  }, height=800)
  
});