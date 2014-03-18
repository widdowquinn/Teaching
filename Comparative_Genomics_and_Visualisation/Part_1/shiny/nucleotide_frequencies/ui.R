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

# We need the Shiny library to use Shiny.
library(shiny);

# Define UI for dataset choice and nucleotide frequency views
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Nucleotide Frequency Plots"),
  
  # Sidebar with activity instructions, and controls to select dataset
  sidebarPanel(
    
    # Activity instructions as Paragraph elements
    p("This interactive example allows you to visualise and compare 
      mono- and dinucleotide plots for a small set of bacterial genomes."),
    p("Use the drop-down box below to select one of the example organisms, 
      all Mycoplasma species, or the whole example set."),
    p("Q1: Which genomes are GC-rich, and which are AT-rich?"),
    p("Q2: Compare the plots for individual Mycoplasma species to those for both 
      Mycoplasma species. What are the differences?"),
    p("Q3: Compare the plots for all genome sets, to the individual species plots.
      What are the differences?"),
    p("Q4: Compare the plot for the 'Unknown' genome to those of the other 
      genomes. To which species does it most likely belong?"),
    
    # Dataset chosen by a dropdown selection box
    selectInput("orgname", "Choose a dataset to load:",
                choices = c("Escherichia coli", "Mycoplasma pneumoniae",
                            "Mycoplasma genitalium", "Mycobacterium tuberculosis",
                            "Mycoplasma spp.", "All", "Unknown")),
    
    # Draw all plots (on slow machines, probably don't want to...)
    checkboxInput(inputId = "all_plots",
                  label = strong("Draw 3-mer and 4-mer plots (slower...)"),
                  value = FALSE)
    
  ),
  
  # Show plot of the mono-, di- and tri-nucleotide frequencies
  # of the selected dataset in the main panel.
  mainPanel(
    plotOutput("freqPlot")
    )
  
));