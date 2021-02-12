#!/bin/bash

# This is the AF look-up script. This will pull allele frequency info from gnomAD or ALFA of variants you are interested in (variant_list.txt) and generate a flat file.

## Depending on how much memory you have available, you may not be able to download the entire gnomad.genomes database. 
## An easy workaround is to download only the chromosomes you're interested in, and then concatenate the results together. See example_gnomad_chunk_method.sh example for more details. 

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
    ## If want gnomAD database, use gsutil.txt. See gsutils.txt for more details.
    ## Python pull for ALFA.
    ## gsutil.txt and pull_database.py are both in parent Merge-tool folder.
python pull_database.py

#2. Edit variant_list.txt to include your variants of interest.

#3. VCF with only variants in list
        # Insert variant list you want
bcftools view --include ID==@variant_list.txt genotype_database.vcf.gz  --output-file variant_1.vcf

# 2. Filter out variables of interest into tsv file
cut -f 1,2,3,4,5,21 variant_1.vcf | awk 'NR > 9  {print}'  > variant_2.tsv

# 3. ALFA ONLY: Separate AC and AN. Divide to create new AF variable. Note: Need to remove commas in AC and AN beforehand, if any.
awk '{sub(/\:/, " ", $6)};1' variant_2.tsv | awk '{print $3, $1, $2, $4, $5, $6, $7, $7/$6}'  > variant_3.tsv

# 4. Add header in. Change input to variant_3.tsv if ran previous step for ALFA. 
echo -e 'ID\tCHROM\tPOS\tREF\tALT\tAC\tAN\tAF\n'  | cat - variant_2.tsv > variant_final.tsv

#gnomAD chunk method example:
## This example uses chrom 1, 6, 10 for gnomAD 2.1.1. Alter as needed. 

# 1.VCF with only variants in list
bcftools view --include ID==@variant_list.txt gnomad.genomes.r2.1.1.sites.1.vcf.bgz  --output-file variant_gnomad_1.vcf
bcftools view --include ID==@variant_list.txt gnomad.genomes.r2.1.1.sites.6.vcf.bgz  --output-file variant_gnomad_6.vcf
bcftools view --include ID==@variant_list.txt gnomad.genomes.r2.1.1.sites.10.vcf.bgz  --output-file variant_gnomad_10.vcf

# 2. Filter out variables of interest into tsv
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' variant_gnomad_1.vcf  > variant_gnomad_1.tsv
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' variant_gnomad_6.vcf  > variant_gnomad_6.tsv
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' variant_gnomad_10.vcf  > variant_gnomad_10.tsv

# 3. Add header in, concatenate all chromosome
        # Copy and paste variables from above query for column names
echo -e "%ID\t%CHROM\t%POS\t%REF\t%INFO/AC\t%INFO/AN\t%INFO/AF\n" | cat - variant_gnomad_1.tsv variant_gnomad_6.tsv > variant_gnomad_final.tsv

