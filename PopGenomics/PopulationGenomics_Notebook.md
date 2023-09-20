# Population Genomics

## Author: Carolyn Hanrahan

### Affiliation: University of Vermont, RSENR

### E-mail contact: [carolyn.hanrahan\@uvm.edu](mailto:carolyn.hanrahan@uvm.edu)

### Start Date: 09/11/2023

### End Date: TBD

### Project Descriptions: This notebook markdown file will document my workflow of the bioinformatics of the Population Genomics section of Ecological Genomics, fall 2023.

# Table of Contents:

-   [Entry 1: 2023-09-11](#id-section1)
-   [Entry 2: 2023-09-13](#id-section2)
-   [Entry 3: 2023-09-18](#id-section3)
-   [Entry 4: 2023-09-20](#id-section4)



<div id='id-section1'/>

### Entry 1: 2023-09-11.

- We reviewed the red spruce study system and exome cature data 
- We discussed the structure of the fastqc files (DNA sequence plus Qscores)
- Using the program FastQC

### Entry 2: 2023-09-13. 

- After discussing the FastQC results, we saw good quality sequence data for most of the read length. 
- The initial 5 bp or so had more variable base frequencies, and the very end of the readshad slightly lower Q-scores. 
- Based on this, we set up an analysis to trim the reads using the fastp program. 
- We ran the bash script `fastp.sh` for this. 
- We looked at the html files produced by `fastp` 

<div id='id-section3'/>

### Entry 3: 2023-09-18.
- We set up this lab notebook
- We ran the mapping.sh file over the weekend for our red spruce populations. 
- Visualize sequence alignment files (*.sam) from our mapping
- Process the files to binary and started running `bam_stats.sh` as well as `process_bam.sh`
- Talked to Nicole about paper selection for class

### Entry 4: 2023-09-20
- We finished processing `bam_stats.sh` by altering the vim file and running `bash bam_stats.sh`. Following this, we examined the read data by using the head command. 
- Created an ANGSD directory and a script for running the ANGSD program to determine genotype likelihood.
- Began running this script using `tmux` in the background. 
