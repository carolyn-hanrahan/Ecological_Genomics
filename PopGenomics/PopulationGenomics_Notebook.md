# [Population Genomics]{style="color:purple;"}

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
-   [Entry 5: 2023-09-25](#id-section5)
-   [Entry 6: 2023-09-27](#id-section6)
-   [Entry 7: 2023-10-02](#id-section7)
-   [Entry 8: 2023-10-04](#id-section8)

<div id='id-section1'/>

### Entry 1: 2023-09-11.

-   We reviewed the red spruce study system and exome cature data
-   We discussed the structure of the fastqc files (DNA sequence plus Qscores)
-   Using the program FastQC

### Entry 2: 2023-09-13.

-   After discussing the FastQC results, we saw good quality sequence data for most of the read length.
-   The initial 5 bp or so had more variable base frequencies, and the very end of the readshad slightly lower Q-scores.
-   Based on this, we set up an analysis to trim the reads using the fastp program.
-   We ran the bash script `fastp.sh` for this.
-   We looked at the html files produced by `fastp`

<div id='id-section3'/>

### Entry 3: 2023-09-18.

-   We set up this lab notebook
-   We ran the mapping.sh file over the weekend for our red spruce populations.
-   Visualize sequence alignment files (\*.sam) from our mapping
-   Process the files to binary and started running `bam_stats.sh` as well as `process_bam.sh`
-   Talked to Nicole about paper selection for class

### Entry 4: 2023-09-20

-   We finished processing `bam_stats.sh` by altering the vim file and running `bash bam_stats.sh`. Following this, we examined the read data by using the head command.
-   Created an ANGSD directory and a script for running the ANGSD program to determine genotype likelihood.
-   Began running this script using `tmux` in the background.

### Entry 5: 2023-09-25

- Today's learning objectives: 
  - Calculate diversity stats for our focal pops (SFS, theta-W, theta-Pi, Tajima’s D)
  - Summarize the results in R and share to google doc
  - Introduce Fst in ANGSD using genotype probabilities

### Entry 6: 2023-09-27

- Today's learning objectives: 
  - Review the diversity stats for our focal pops on the google doc
  - Estimate genetic differentiation (Fst) in ANGSD between our focal red spruce pops and black spruce
  - Visualize population structure using PCA (principle component analysis) and Admixture
  
### Entry 7: 2023-10-02

- Today is a coding and bioinformatics day. Goals include: 
  - Review the population structure results
  - Perform a scan for selection and identify contigs with outlier loci
  - Identify and visualize outlier loci
  
## Entry 8: 2023-10-04

- Today's goals include:
  - Revisit pcANGSD selection scan to output genetic PCs and outlier loci list
  - Extract climate data for red spruce localities and summarize with PCA
  - Run genotype-environment association (GEA) analysis running
  
- Notes: 
  - Genotype data (SNPs) + environmental data (climate) --> genome wide average 
  - Env (y) = Genotype (locus specific) + covariates (genetic PCs)

- Steps for today: 
  1. in R: revisit selection outliers --> outlier list (PC1)
  2. genetic PCA --> genetic PC1 + PC2 as covariates
  3. Get bioclim environment variables 
  4. Transfer files to server -- run GEA. 
  
  ## Entry 9: 2023-10-09
  
- Transcriptomics unit; reading paper on Colorado Potato Beetle and insecticide resistance. 
  
- Coding goals for the day: 
  - Review Acartia hudsonica ecology and biogeography and the experimental evolution/transcriptomics experimental design.

  - Develop questions that can be addressed and hypotheses that can be tested with the A. hudsonica experiment.

  - Understand the general work flow or “pipeline” for processing and analyzing RNAseq data.

  - Visualize and interpret the quality of our Illumina data.

  - Assess our previously assembled de novo transcriptome assembly using Trinity.
Start mapping reads and quantifying abundance simultaneously using Salmon.