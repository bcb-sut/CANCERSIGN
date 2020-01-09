# CANCERSIGN manual
M. Bayati, H.R. Rabiee, et al., and H. Alinejad-Rokny, “**CANCERSIGN: a user-friendly and robust tool for identification and classification of mutational signatures and patterns in cancer genomes**”, preparing for submission.


Prerequisites
============

Operating System
-------------------------
The recommended operating systems for using this tool are **Linux** and **MacOSX**.

R
-------------------------
**R** is required to be installed and the command “**Rscript**” must be available. In addition, the following packages must be installed for **R**:
* BSgenome.Hsapiens.UCSC.hg19
* data.table
* doParallel
* ggplot2
* configr

Initial setup
============
Download CANCERSIGN.gz and uncompress it:
```
$ gunzip </path/to/CANCERSIGN.gz>
```
Inside CANCERSIGN directory, run the following commands to activate ```cancersign```:
```
$ cd src
$ chmod +x cancersign
```
And then add the path to the CANCERSIGN **src** directory to the PATH variable: 
```
$ export PATH=$PATH:</path/to/CANCERSIGN/src>
```
You can edit ```.bashrc``` or any other appropriate shell script to set the above configuration persistently.

Run CANCERSIGN
============

Before running CANCERSIGN, the user has to provide a simple text file which contains the custom configurations including the path to the input data, the path to the output directory where the user wants the results to be stored in and the analysis parameters (guidelines for input data format and writing the configuration file are provided in the subsequent sections). Assuming that the configurations are written in a file named config.txt, CANCERSIGN starts the analyses with the following command:
```
$ cancersign --config </path/to/config.txt>
```
where </path/to/config.txt> is the path to the configuration file.



Input data for CANCERSIGN
============
The input data must be a tab-delimited file with the following fields:

1.  **sample_id** - The ID of the sample.
2.  **chromosome** - The name of the chromosome (chr1, chr2, chr3, …, chrX, chrY or chrM).
3.  **position** - The position of the mutation in the chromosome.
4.  **reference** - The nucleotide at the corresponding location on the reference genome.
5.  **mutated_to** - The mutated nucleotide at the corresponding location in the sample.
