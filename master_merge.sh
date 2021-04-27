#!/bin/bash

#SBATCH -c 1                    # Request one core
#SBATCH -N 1                    # Request one node

#SBATCH -t 0-12:00              # Run time D-HH:MM format
#SBATCH -p short                # Partition to run in
#SBATCH --mem=2000              # Memory total in MB (in all cores)
#SBATCH -e /home/gao4/clinvar-master/error/hostname_%j.err  # File to which STDERR will be written, including jobID      # Change/Delete this
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN, END, FAIL, ALL
#SBATCH --mail-user=grace.obrien@childrens.harvard.edu  # Email where notifications send to                              # Change/Delete this

# Loading required modules
module load gcc/6.2.0
module load bcftools
module load htslib/1.10.2

# 1. Tabix both files to provide sorted and indexed file for the merge. Output is tabixed VCF files (.tbi).

tabix -p vcf $1
tabix -p vcf $2

# 2. Merge genotype and phenotype files by "annotating" ClinVar  with the variables from gnomAD that want. Output is VCF file.

bcftools annotate -a $1  -c CHROM,POS,REF,ALT,INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt,INFO/AC_popmax,INFO/AN_popmax,INFO/AF_popmax,INFO/AC_afr,INFO/AN_afr,INFO/AF_afr $2  > temp_merge_1.vcf


# 3. Query variables that you want. Add additional variables here. Output is TSV file.

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ALLELEID\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/nhomalt\t%INFO/CLNDN\t%INFO/CLNREVSTAT\t%INFO/CLNSIG\t%INFO/CLNVC\t%INFO/GENEINFO\t%INFO/ORIGIN\n' temp_merg$


# 4. Add header in. Make sure column names align with the order that they were queried. Output is TSV file.

echo -e "CHROM\tPOS\tREF\tALT\tALLELEID\tAC\tAN\tAF\tnhomalt\tCLNDN\tCLNREVSTAT\tCLNSIG\tCLNVC\tGENEINFO\tORIGIN\n" | cat - temp_merge_2.tsv > temp_merge_3.tsv


# 5.A Filter file to:
                #a. Exclude variants that are missing genotype info (chose AC to filter here)
                #b. Include clinical review status of CLNREVSTAT= criteria provided, multiple submitters, no conficts OR reviewed by expert panel
                #c. Include clinical significance of CLNSIG= pathogenic, likely pathogenic OR likley pathogenic/pathogenic (this code grabs any CLNSIG value that conatins "athogenic")

awk 'NR==1 || !($6 == ".") && !($6 == "0") && (($11=="criteria_provided,_multiple_submitters,_no_conflicts") || ($11== "reviewed_by_expert_panel")) && ($12 ~/athogenic/)' temp_merge_3.tsv > temp_merge_4.tsv

# OPTIONAL 5.B  Filter by origin- currently commented out because easier to check for germline manually in final output. Can use this filtration step instead.
                # This code will include only the following ORIGIN values:
                                #a. germline=1
                                #b. germline,somatic=1+2=3
                                #c. germline,somatic,maternal=1+2+16=19

#awk 'NR==1 || (($15==1) || ($15==3) || ($15==16))' temp_merge_4.tsv > temp_merge_4.tsv
# 6. Generate the final output file.
        #a. Include gnomAD database type (exome vs genome) in final output filename
        #b. Filter by GENEINFO provided in command line. Include GENEINFO name in final outuput filename
        #c. Filter by CLNDN (clinical disease name) if provided in command line. Include CLNDN name in final output filename

if [[ $1  =~ ^gnomad.exome ]]; then

for i in $3;
do
awk  'NR==1 || ($14 ~/'$i'/)' temp_merge_4.tsv  > clinvar_gnomad_exome_${i}_$(date +%F).tsv;
done

if (( $4 = 1 )); then
for j in $4;
do
awk 'NR==1 || ($10 ~/'"$j"'/)' clinvar_gnomad_exome_${i}.tsv > clinvar_gnomad_exome_${i}_${j}_$(date +%F).tsv;
done

fi
fi


if [[ $1 =~ ^gnomad.genome ]]; then

for i in $3;
do
awk  'NR==1 || ($14 ~/'$i'/)' temp_merge_4.tsv  > clinvar_gnomad_genome_${i}_$(date +%F).tsv;

if (( $4 = 1 )); then
for j in $4;
do
awk 'NR==1 || ($10 ~/'$j'/)' clinvar_gnomad_genome_${i}.tsv > clinvar_gnomad_genome_${i}_${j}_$(date +%F).tsv;
done

fi
done
fi


