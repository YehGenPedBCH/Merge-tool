#!/bin/bash

#SBATCH -c 1                    # Request one core
#SBATCH -N 1                    # Request one node

#SBATCH -t 0-12:00              # Run time D-HH:MM format
#SBATCH -p short                # Partition to run in
#SBATCH --mem=2000              # Memory total in MB (in all cores)
#SBATCH -e /home/gao4/clinvar-master/src/error/hostname_%j.err  # File to which STDERR will be written, including jobID
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN, END, FAIL, ALL
#SBATCH --mail-user=grace.obrien@childrens.harvard.edu  # Email where notifications send to

# Loading required module
module load gcc/6.2.0 samtools/1.10
module load vcftools
module load htslib/1.10.2
module load python
pip install urllib


# 1. Run pull script. This pulls most updated ALFA database. It is GR38
#python ALFA_pull.py


# Add in date
#mv ALFA_freq.vcf.gz ALFA_freq_$(date +%F).vcf.gz

# 2.  Reformatting the chromosome variable
#tabix -p vcf $1

# 1.VCF with only variants in list

#bcftools view --include ID==@$2 $1  --output-file variant_1.vcf


#bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%INFO/SAMN10492695\n' variant_1.vcf  > variant__1.tsv

        # Filter out variables of interest into tsv file

cut -f 1,2,3,4,5,21 variant_1.vcf | awk 'NR > 9  {print}'  > variant_1.tsv

        # Separate AC and AN. Divide to create new AF variable. Note: Need to remove commas in AC and AN beforehand.

#awk '{sub(/\:/, " ", $6)};1' variant_1.tsv | awk '{print $3, $1, $2, $4, $5, $6, $7, $7/$6}'  > variant_2.tsv

        # Add header in
#echo -e 'ID\tCHROM\tPOS\tREF\tALT\tAC\tAN\tAF\n'  | cat - variant_2.tsv >  ALFA_$(date +%F).tsv

