



# CANCERSIGN manual
M. Bayati, H.R. Rabiee, et al., and H. Alinejad-Rokny, “**CANCERSIGN: a user-friendly and robust tool for identification and classification of mutational signatures and patterns in cancer genomes**”, preparing for submission.

1: Prerequisites and notes
============
> **_NOTE:_**  This tool performs all analyses based on hg19 genome build.

Operating system
-------------------------
The recommended operating systems for using this tool are **Linux** and **MacOSX**.

R
-------------------------
This pipeline requires **R** version 3.6 or higher and needs the following packages to be installed:
* BSgenome.Hsapiens.UCSC.hg19
* data.table
* doParallel
* ggplot2
* configr

2: Initial setup
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

3: Run CANCERSIGN
============

Before running CANCERSIGN, the user has to provide a simple text file which contains the custom configurations including the path to the input data, the path to the output directory where the user wants the results to be stored in and the analysis parameters (guidelines for input data format and writing the configuration file are provided in the subsequent sections). Assuming that the configurations are written in a file named config.txt, CANCERSIGN starts the analyses with the following command:
```
$ cancersign --config </path/to/config.txt>
```
where `</path/to/config.txt>` is the path to the configuration file.



4: Input data for CANCERSIGN
============
The input data must be a tab-delimited file with the following fields:

1.  **sample_id** - The ID of the sample.
2.  **chromosome** - The name of the chromosome (chr1, chr2, chr3, …, chrX, chrY or chrM).
3.  **position** - The position of the mutation in the chromosome.
4.  **reference** - The nucleotide at the corresponding location on the reference genome.
5.  **mutated_to** - The mutated nucleotide at the corresponding location in the sample.


5: Creating configuration file
============
The configuration file is a simple text file with lines defining the settings required for the analyses. The first two lines must define the path to the input file and the path to the output directory where the user wants the results of analyses to be stored in. These two lines are written as follows:

```
input_file = </path/to/input_data_file>
output_dir = </path/to/output_dir>
```
where `</path/to/input_file>` and `</path/to/output_dir>` are paths to the input data file and to the output directory respectively.

The rest of the lines in the configuration file determine the desired types of analyses as well as the parameters for those analyses. In short, one line is written to “enable” each desired analysis and it is followed by other lines which set the corresponding parameters for that analysis. The following sections explain the configurations for each analysis type.

5.1: Infer 3-mer mutational signatures
-------------------------
To enable this analysis, write the following line in the configuration file:
```
infer_3mer_signatures = yes
```
All parameters for this analysis have initial **default** values unless the user specifies them with the following lines:

`N_min_3mer = <a number, default: 1>` 
*The minimum number for testing the number of signatures*

`N_max_3mer = <a number, default: 10>`
*The maximum number for testing the number of signatures*

`CPU_3mer = <a number, default: 30>`
*The number of allocated CPU cores for this analysis*

`nmf_iters_3mer = <a number, default: 1e4>`
*Number of iterations for NMF algorithm in each epoch*

`nmf_conv_3mer = <a number, default: 1e-5>`
*Convergence threshold for stopping NMF iterations*

`nmf_max_3mer = <a number, default: 5e5>`
*Maximum number of iterations for NMF algorithm*

`boot_iters_3mer = <a number, default: 30 or number of available CPUs>`
*Number of bootstrap iterations performed in each epoch*

`boot_conv_3mer = <a number, default: 0.01>`
*The convergence threshold for stopping bootstrap iterations*

`boot_max_3mer = <a number, default: 600>`
*Maximum number of bootstrap iterations*

The results of this analysis are stored in the output directory in a folder named **infered_3mer_signatures**.


5.2: Infer 5-mer mutational signatures
-------------------------
To enable this analysis, write the following line in the configuration file:
```
infer_5mer_signatures = yes
```
The parameter that the user must provide for this analysis is the desired set of 3-mer motifs. Based on this parameter, the tool expands the specified 3-mer motifs to all possible 5-mer motifs (which contain the specified 3-mer motifs) and then infers the 5-mer mutational signatures corresponding to them. This parameter is specified with the following line:

`selected_3mer_motifs_for_5mer_signatures = "motif-1", . . ., "motif-n"`
where `"motif-i"` is a selected 3mer motif in a standard format. An example of this standard format is `G[C>T]A` which means `C>T` mutation with `G` as left flanking nucleotide and `A` as the right flanking nucleotide.

All other parameters for this analysis have initial **default** values unless the user specifies them with the following lines:

`N_min_5mer = <a number, default: 1>` 
*The minimum number for testing the number of signatures*

`N_max_5mer = <a number, default: 10>`
*The maximum number for testing the number of signatures*

`CPU_5mer = <a number, default: 30>`
*The number of allocated CPU cores for this analysis*

`nmf_iters_5mer = <a number, default: 1e4>`
*Number of iterations for NMF algorithm in each epoch*

`nmf_conv_5mer = <a number, default: 1e-5>`
*Convergence threshold for stopping NMF iterations*

`nmf_max_5mer = <a number, default: 5e5>`
*Maximum number of iterations for NMF algorithm*

`boot_iters_5mer = <a number, default: 30 or number of available CPUs>`
*Number of bootstrap iterations performed in each epoch*

`boot_conv_5mer = <a number, default: 0.01>`
*The convergence threshold for stopping bootstrap iterations*

`boot_max_5mer = <a number, default: 600>`
*Maximum number of bootstrap iterations*

The results of this analysis are stored in the output directory in a folder named **infered_5mer_signatures**.



