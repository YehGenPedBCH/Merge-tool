#!/bin/bash

#SBATCH -c 1                    # Request one core
#SBATCH -N 1                    # Request one node

#SBATCH -t 0-12:00              # Run time D-HH:MM format
#SBATCH -p short                # Partition to run in
#SBATCH --mem=2000              # Memory total in MB (in all cores)
#SBATCH -e /home/gao4/clinvar-master/error/hostname_%j.err  # File to which STDERR will be written, including jobID      # Change/Delete this
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN, END, FAIL, ALL
#SBATCH --mail-user=grace.obrien@childrens.harvard.edu  # Email where notifications send to                              # Change/Delete this


# Loading required module
#cd /home/gao4/clinvar-master/src
module load gcc/6.2.0 samtools/1.10
module load vcftools
module load bcftools
module load htslib/1.10.2

# 1. Tabix vcf file so can use BCFtool commands

#tabix -p vcf $1


# 2. Reformat and name output, depending on input genotype database

# gnomad exome
if [[ $1  =~ ^gnomad.exome ]]; then
        # Pull only variants in txt file from database. This is the slow step- can take over 1 hr

bcftools view --include ID==@$2 $1  --output-file gnomad_exome_AF_1.vcf

        # Filter out variables of interest into tsv

bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' gnomad_exome_AF_1.vcf  > gnomad_exome_AF_1.tsv

        # Add header in
        # Copy and paste variables from above query for column names
echo -e "%ID\t%CHROM\t%POS\t%REF\t%INFO/AC\t%INFO/AN\t%INFO/AF\n" | cat - gnomad_exome_AF_1.tsv > gnomad_exome_AF_$(date +%F).tsv

fi


# gnomad genome
if [[ $1  =~ ^gnomad.genome ]]; then

         # Pull only variants in txt file from database. This is the slow step- can take over 1 hr

bcftools view --include ID==@$2 $1  --output-file gnomad_genome_AF_1.vcf

        # Filter out variables of interest into tsv
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' gnomad_genome_AF_1.vcf  > gnomad_genome_AF_1.tsv

        # Add header in
        # Copy and paste variables from above query for column names
echo -e "%ID\t%CHROM\t%POS\t%REF\t%INFO/AC\t%INFO/AN\t%INFO/AF\n" | cat - gnomad_genome_AF_1.tsv > gnomad_genome_AF_$(date +%F).tsv

fi


# ALFA
if [[ $1  =~ ^ALFA ]]; then
         # Pull only variants in txt file from database. This is the slow step- can take over 1 hr

#bcftools view --include ID==@$2 $1  --output-file ALFA_AF_1.vcf

        # Filter out variables of interest into tsv file, change cut numbers if want more columns

cut -f 1,2,3,4,5,21 ALFA_AF_1.vcf | awk '/#CHROM/,/end/  {print}'  > ALFA_AF_1.tsv

        # Remove "," in AC:AN column, if any
awk '{gsub(/,0/, "", $6)}1' ALFA_AF_1.tsv > ALFA_AF_2.tsv
awk '{gsub(/0,/, "", $6)}1' ALFA_AF_2.tsv > ALFA_AF_3.tsv

        # Separate AC and AN. Divide to create new AF variable.

awk NR\>1 ALFA_AF_3.tsv | awk '{sub(/\:/, " ", $6)};1' | awk '{print $3, $1, $2, $4, $5, $6, $7, $7/$6}'  > ALFA_AF_4.tsv

        # Add header in

echo -e 'ID\tCHROM\tPOS\tREF\tALT\tAC\tAN\tAF\n'  | cat - ALFA_AF_4.tsv >  ALFA_AF_$(date +%F).tsv

fi
