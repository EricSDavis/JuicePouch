# JuicePouch

## Overview

JuicePouch is a collection of scripts that extend the functionality of Juicebox, straw, and Juicer Tools.


## Tools (At a glance)

### extractCounts

This script extracts reads from a list of `.hic` files corresponding to a bedpe file using the straw api. Counts for each `.hic` file are appended as columns to the bedpe file and sent to a new output file.

#### Usage

1. Load python modules
    ```bash
    module load python/3.6.6
    ```

2. Install required python libraries (i.e. straw)
    ```bash
    python3 -m pip install hic-straw
    ```
3. Launch `extractCounts.py` with a sbatch wrapper
    ```bash
    sbatch -p general -t 1440 --wrap='python3 path/to/extractCounts.py -l <LOOPFILE> -f <HICFILES> -o <OUTFILE>'
    ```
    
    `<LOOPFILE>` is a BEDPE file (i.e. `.txt` or `.bedpe`)  
    `<HICFILES>` is a `.hic` file or list of `.hic` files  
    `<OUTFILE>` is a user-specified output file (default: `loopsOut.txt`)


#### Dependencies
* BEDPE formatted file
* .hic file(s)


### extractNorm

This script will extract the normalization factors for a Hi-C file for all chromosomes and resolutions. It uses juicer tools (juicer dump norm) in a parallelized fashion for speed.

#### Dependencies
* [Java](https://www.java.com/en/)
* [Juicer Tools](https://github.com/aidenlab/juicer/wiki/Download)
* .hic file