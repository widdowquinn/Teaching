# R script to install R packages used in the Comparative Genomics and 
# Visualisation course/workshop

# Set mirror
local({r = getOption("repos")
       r["CRAN"] = "http://cran.r-project.org"
       options(repos=r)
})
install.packages("shiny")
install.packages("ggplot2")
install.packages("hash")
install.packages("gridExtra")
install.packages("reshape2")
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("Biostrings")