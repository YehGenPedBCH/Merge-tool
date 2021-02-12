#!/bin/bash

# This script will 

#SBATCH -c 1                    # Request one core
#SBATCH -N 1                    # Request one node

#SBATCH -t 0-12:00              # Run time D-HH:MM format
#SBATCH -p short                # Partition to run in
#SBATCH --mem=2000              # Memory total in MB (in all cores)

# Loading required modules
module load gcc/6.2.0 samtools/1.10
module load vcftools
module load bcftools
module load htslib/1.10.2
module load python 

# 1. Run pull_database.py script. This pulls both ClinVar and/or ALFA databases.
python pull_database.py

# 2. OPTIONAL: Subset ClinVar vcf file to only the chromosome(s) of interest. This step will make the merge/annotation faster.
bcftools view --regions 17 clinvar.vcf.gz --output clinvar_17.vcf.gz

# 3. ALFA only: Additional steps are needed here. 
Add here.

# 4. Tabix both files to provide index file for the merge
tabix -p vcf gnomad.exomes.r2.1.1.sites.17.vcf.bgz
tabix -p vcf clinvar_17.vcf.gz

# 5. Merge: This steps ANNOTATES the ClinVar vcf with the genotype data of interest from gnomAD (i.e AF, AC, AN etc.).
bcftools annotate -a gnomad.exomes.r2.1.1.sites.17.vcf.bgz -c CHROM,POS,REF,ALT,INFO/AF clinvar_GR37.vcf.gz > annotate_1.vcf

# 6. Extracting dataset that we want. 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ALLELEID\t%INFO/AF\t%INFO/CLNDN\t%INFO/CLNREVSTAT\t%INFO/CLNSIG\t%INFO/CLNVC\t%INFO/GENEINFO\t%INFO/ORIGIN\n' annotate_1.vcf > annotate_2.tsv

# 7. Add header in. Copy and paste names from above. Note if want to use R to filter- stop here and upload annotate_2_header.tsv to R.
echo -e "CHROM\tPOS\tREF\tALT\tALLELEID\tAF\tCLNDN\tCLNREVSTAT\tCLNSIG\tCLNVC\tGENEINFO\tORIGIN\n" | cat - annotate_2.tsv > annotate_2_header.tsv

# 8. OPTIONAL: Filter by chromosome, gene, pathogenicity, gold stars, and germline status. Alter as needed.
awk 'NR==1 || ($1~/17/) && !($6 == ".") && !($6 == "0") && (($8=="criteria_provided,_multiple_submitters,_no_conflicts") || ($8 == "reviewed_by_expert_panel")) && ($9 ~/athogenic/) && ($11 ~/TP53/) && (($12==1) || ($12==3))' annotate_2_header.tsv > annotate_TP53.tsv





