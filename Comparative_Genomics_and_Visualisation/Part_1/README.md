#README.md - Part_1

## Overview

This repository contains teaching materials used in delivering the first part of the 3rd year bioinformatics course lectures on "Comparative Genome Analysis and Visualisation".

## Getting Started

You can grab a local copy of all the files for this workshop/lecture using `git clone`:

```
$ git clone https://github.com/widdowquinn/Teaching.git
```

The lecture is presented in a tutorial/workshop format, with interactive examples and exercises using [iPython notebooks](http://ipython.org/notebook.html). 

* Slides for the talk are provided in the `presentation` subdirectory, in Powerpoint `.pptx`, `.ppt` and `.pdf` format. 
* Example data can be found in the `data` directory.
* iPython notebooks are found in the directory containing this `README.md` file

### Prerequisites

The exercises and examples have been written, and are known to run, with the following software, but they may also run happily on other versions:

* **Python 2.7**
* **iPython (with pylab and matplotlib) 1.1.0**
* **Biopython 1.63+**
* **Pandas 0.13**

## Executing iPython notebooks

To start iPython in a suitable form in your browser, execute

```
$ ipython notebook --pylab inline
```

at the command-line. It is possible to start iPython in other ways (e.g. in the terminal window, without inline plots, or using the Qt console), but this is the way the course examples/exercises were intended to be run.
