# README.md - find_rbbh

## Overview

This directory contains files for in-class activity centred around reciprocal best BLAST hits (RBBH). The activities and explanation are in the iPython notebook `find_rbbh.ipynb`, which should be started using:

```
ipython notebook --pylab inline
```

and selecting `find_rbbh` in the browser window that appears.

The shell script `run_rbbh_blast.sh` should already have been run to generate data in the `rbbh_output` directory. If this is missing, then issuing

```
./run_rbbh_blast.sh
```

should regenerate the output.

## Prerequisites

The notebook was written and tested with the following software versions, though other versions may also work.

* **BLAST+** 2.2.29+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>
* **Python** 2.7 <http://www.python.org>
* **iPython** 1.1.0 (with Matplotlib and SciPy/Pylab) <http://ipython.org/>
* **Biopython** 1.63 <http://www.biopython.org>
* **Pandas** 0.13.0 <http://pandas.pydata.org/>