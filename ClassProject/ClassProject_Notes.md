# [Class Project]{style="color:purple;"}

## Author: Carolyn Hanrahan

### Affiliation: University of Vermont, RSENR

### E-mail contact: [carolyn.hanrahan@uvm.edu](mailto:carolyn.hanrahan@uvm.edu)

### Start Date: 11/13/2023

### End Date: End of semester!

### Project Descriptions: A place for notes on my group project.

###Entry 1:

Goals:

-   Give opportunity to dig deeper into an area of interest
-   Pose new questions
-   Explore different analyses
-   Integrate across different datasets
    -   how can you integrate our datasets with datasets online?
-   Environmental data (Bioclim), phenotype data, etc....

Expectations:

-   Working collaboratively in teams
-   work with one or more of the class datasets (red spruce, copepods, or purple sea urchins)
-   Present initial plan for feedback Wednesday (second half of class)
-   continue brainstorming on Wednesday --> share pitch to class during second half of class on Wednesday - Devise and implement an analysis strategy ("divide and conquer")
-   Post-Thanksgiving: 2 weeks of open lab session. Come in and work together on your analyses.

Review of the datasets:

1.  Red Spruce dataset:

-   Exome-capture data for 12 populations (N=95 individuals across the range).
-   Approximately 16 Black Spruce individuals
-   Asking questions about population structure on the landscape and their associations with climate gradients. Also discussed hybridization.
-   Future directions...:
    -   introgression mapping between Black and Red Spruce populations.
    -   deeper dive into selection outliers --> what is the ancestral source?
    -   Integrating the Black Spruce dataset more to determine how selection is acting.
    -   Mechanisms for rate of change between different geographic regions (rates of adaptation vs. gene flow in southern fragmented populations vs. more continuous northern populations)
- Other data to incorporate: 
  - Survival data, diameter, height growth (measures of fitness), measures of cold injury (how resilient the foliage is to freezing mid-winter)
  - RNASeq data: growth chamber experiment with control, heat, and drought + warming treatment group. 73 total replicates. (you could do GWAS analysis to link genotype to phenotype for this data)
  - Reference genome for Black Spruce (much closer reference)
  
2. Copepod dataset: 

- Acartia Hudsonica, cold-water adapted copepod species. Experimentally evolved for 12 generations (F0-F11)
- 4 treatment groups (AM, OA, OW, OWA)
- Future directions..
  - Exploring plasticity - loss/gain? 
  - You would expect to see plasticity in the F0 generations. Whatever is happening in the first generation is responding to their conditions. Looking for loss of plasticity across time... 
  - Examining assimilation (if gene expression remains consistent across time) vs. compensation (when gene expression changes across time)
  - Allele frequency of RNASeq data
  
  
Additional datasets: 
  - A. tonsa RNASeq data F15 reciprocal transplant; transplanted OWA animals to ambient conditions and vice versa, and tracked gene expression data for three generations 
  - capture-seq data (sequenced from genomic DNA) - gDNA, all treatment groups at F25. Also an initial sample at F0. 
  
3. Purple Urchin dataset 

- 20 individuals from 7 different populations; total of 140 sequenced urchins in total from the West Coast. High gene flow between populations. Good model for studying local adaptation! 
- Sequenced whole genome (low coverage, WGS)
- Sequence means mapped to the reference genome, bcf files created. 
- You can look at structural variation, or you can look at just SNP data. 

Future directions: 
- take bcf file and filter for SNPs. Examine different environmental variables to determine where SNPs correlate. This has already been done for pH variability. Examine other environmental variables (temperature, salinity, etc.)
  - which alleles correlate with these environmental variables 
- Examine linkage (R package called `plink`)
  - Can examine linkage and structural variants. 
- Look at all the variation in the bcf
  - phenotype data
- Continue structural variation analysis: look at bcf file; examine insertions, distribution across genome (look for hotspots in structural variation)
- Local PCA, but for all chromosomes. Collapse bcf files into one bcf file. MDS plot would compare windows between different chromosomes. 
- Alter parameters. Define window size by base pairs instead of SNPs. 

Other datasets: 
- Differential gene expresssion 
- Single generation selection experiments. 
- Environmental data (found online --> temperature  gradient, bioclim variables for ocean data) 
  