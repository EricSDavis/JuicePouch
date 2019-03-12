# JuicePouch

## Overview
**********************

JuicePouch is a collection of scripts that extend the functionality of Juicebox, straw, and Juicer Tools.


## Tools (At a glance)
**********************

### extractNorm

This script will extract the normalization factors for a Hi-C file for all chromosomes and resolutions. It uses juicer tools (juicer dump norm) in a parallelized fashion for speed.

### Dependencies
* [Java](https://www.java.com/en/)
* [Juicer Tools](https://github.com/aidenlab/juicer/wiki/Download)
* .hic file


### extractCounts

This script extracts reads from a .hic file corresponding to a bedpe file using the straw api. Reads are appended to the end of the input file as a new output file.

### Dependencies
* [Straw API](https://github.com/aidenlab/straw/wiki)
* .hic file
*BEDPE formatted file
