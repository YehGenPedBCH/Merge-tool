# AF look-up tool
Merge phenotype database (ClinVar) with genotype databases (gnomAD, ALFA)

# Overview
The purpose of the AF look-up tool is to obtain estimates of allele frequencies and other genotype information of a list of variants from gnomAD and/or ALFA databases. 

This tool uses python to pull ALFA databases, gsutils to pull gnomAD databases, and BCFtools to format and query AF from databases. This can all be completed using a linux-based platform. 

Add more. (what command line code is for it after make changes).

# General pipeline
1. Download databases of interest to linux-based platform.
      - gnomAD: https://gnomad.broadinstitute.org/downloads
         OR
      - ALFA: https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/
2. Add more/
