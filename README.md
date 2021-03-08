# Merge-tool
Merge phenotype database (ClinVar) with genotype databases (gnomAD, ALFA)

# Overview
The purpose of the merge tool is to obtain estimates of allele frequencies for the pathogenic and likely pathogenic variants that will later be used as input parameters for a decision model. It pulls phenotype data from ClinVar and genotype data from gnomAD or ALFA to generate a dataset of variants for the model.

This tool uses python to pull ClinVar and ALFA databases, gsutils to pull gnomAD databases, and BCFtools to format and merge the databases. This can all be completed using a linux-based platform. 

This tool is fast and eliminates the need to manually curate this list of variants for models.

# General pipeline
1. Download databases of interest to linux-based platform.
      - ClinVar: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
      - gnomAD: https://gnomad.broadinstitute.org/downloads
         OR
      - ALFA: https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/
2. Tabix and reformat files. Additional steps needed here for ALFA (chromosome name notation needs to be changed).
3. Annotate ClinVar file with gnomAD or ALFA. 
4. Filter out necessary criteria. Additional steps needed here for ALFA (computational step to separate out the AC:AN column). There are different scripts depending on what variables you want outputted (i.e. race-specific AF).
5. Export to R and filter further, as needed.

# Useful links:
1. Differences in gnomAD versions: https://gnomad.broadinstitute.org/faq
      - For merge tool projects, we use gnomad.exome database. 
      - For AF look-up projects, we use gnomad.genome database. 
2. ALFA overview: https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/
3. BCFtools documentation: http://samtools.github.io/bcftools/bcftools.html

