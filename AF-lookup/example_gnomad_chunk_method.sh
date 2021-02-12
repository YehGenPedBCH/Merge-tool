#!/bin/bash

# This is the AF look-up script for gnomAD chunking. Only download gnomAD databases of chromosomes that you're interested in. 
## This example uses chrom 1, 6, 10 for gnomAD 2.1.1. Alter as needed. 

#SBATCH -c 1                    # Request one core
#SBATCH -N 1                    # Request one node

#SBATCH -t 0-05:00              # Run time D-HH:MM format
#SBATCH -p short                # Partition to run in
#SBATCH --mem=200               # Memory total in MB (in all cores)

# Loading required modules
module load gcc/6.2.0 samtools
module load vcftools
module load htslib/1.10.2
module load bcftools
module load python

# 1. Pull databases.
    ## See gsutils.txt for more details.

#2. Edit variant_list.txt to include your variants of interest.

# 3.VCF with only variants in list
bcftools view --include ID==@variant_list.txt gnomad.genomes.r2.1.1.sites.1.vcf.bgz  --output-file variant_gnomad_1.vcf
bcftools view --include ID==@variant_list.txt gnomad.genomes.r2.1.1.sites.6.vcf.bgz  --output-file variant_gnomad_6.vcf
bcftools view --include ID==@variant_list.txt gnomad.genomes.r2.1.1.sites.10.vcf.bgz  --output-file variant_gnomad_10.vcf

# 2. Filter out variables of interest into tsv
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' variant_gnomad_1.vcf  > variant_gnomad_1.tsv
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' variant_gnomad_6.vcf  > variant_gnomad_6.tsv
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' variant_gnomad_10.vcf  > variant_gnomad_10.tsv

# 3. Add header in, concatenate all chromosome
        # Copy and paste variables from above query for column names
echo -e "%ID\t%CHROM\t%POS\t%REF\t%INFO/AC\t%INFO/AN\t%INFO/AF\n" | cat - variant_gnomad_1.tsv variant_gnomad_6.tsv variant_gnomad_10.vcf > variant_gnomad_final.tsv
