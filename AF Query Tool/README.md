# Allele Frequency Query Tool Documentation 04/23/2021

## Purpose:
The purpose of this tool is to extract allele frequencies (AF) and other genotype information from genotype databases (gnomAD, ALFA) to generate a flat file of variants of interest and their genotype information.

## Overview:
The allele frequency query tool was developed to pull AF of specific variants of interest from different genotype databases. This tool eliminates the need to manually look up the allele frequency and genotype details for a long list of SNPs. 
The tool uses python to pull ALFA database, gsutils to pull gnomAD databases, and bash script  and BCFtools to format and pull info from the databases. This can all be completed using a Linux/Unix-based platform.
Database Background:
GnomAD (Genomic Aggregation Database) is a public archive of exome and whole-genome sequences from disease-specific and population genetic studies [1]. It contains genotype information, including allele frequency estimates of the overall population and by specific subpopulations. There are two gnomAD databases that are available: genome and exome. The most recent version available for it is v3.1.1, updated November 2021. The Broad Institute recommends using the exome database for coding region analyses and the genome database for non-coding region analyses (https://gnomad.broadinstitute.org/faq).  

The Allele Frequency Aggregator (ALFA) is anopen access database that contains allele frequencies for variants in dbGaP [2]. It contains population-specific genotype data for over one million dbGaP subjects.
The gnomAD and ALFA databases can be downloaded as VCF files. The AF query tool uses BCFtools, a program for querying, sorting, and manipulating VCF files, to perform operations on these genomic databases [3]. The primary commands we use in the AF query tool are view and query. For further information on BCFtools, visit the source code and documentation sites (https://www.htslib.org/, https://github.com/samtools/bcftools).  

## Setup:
To use this tool, you need access to a Linux/Unix-based computing platform. Download the github repository and install the required programs (https://github.com/graceannobrien/Merge-tool/AF_query). See General Pipeline for more details on customizing these functions.

## Input formats:
This tool works for genotype databases that are VCF only. 

## Output formats:
The output is a flat file of AF estimates for the gene of interest. 

## Default Settings:
The default settings for the allele frequency query tool output files contain the following variables:
![image](https://user-images.githubusercontent.com/67425562/116264537-5c1fbd80-a748-11eb-8cc3-6ab16d5f197c.png)

*Note: In ALFA database, SAMN10492705 variable is the total across all population’s AC:AN ratio. The AF query tool manipulates the data to output AC and AN separately and computes AF.

If you wish to alter these, see the Data Dictionaries for each database, located in the git repository, to see the variable options. Note: the allele frequency query tool code will need to be altered to include those variables as additional columns. 

## General Pipeline:
1.	Curate a list of SNPs that they are on that you are interested in.
2.	Download genotype database(s) (gnomAD, ALFA). Depending on your platform’s storage level, you may be able to download the entire gnomAD database. Since Harvard’s O2 does not allow that large of database to be downloaded, we needed to take additional steps to only download the chromosomes we were interested in. Use ALFA_pull.sh script to pull ALFA and gsutils to pull gnomAD.
3.	Perform AF query step. In the command line, you need to specify the genotype database and the text file of variants. Depending on the number of variants in your list, it will take approximately 1-2 hours to run the query step. 
4.	Note: there are default conditions set. See Default Settings and how to change if needed.
5.	Export flat file to desktop.

## Example command line code:
In this example, we want to produce a file with the allele frequencies of a list of variants associated with cardiomyopathy in childhood cancer survivors. 

1.	Create text file with the following SNPs of interest. Name cardio_variants_2021-04-23.txt.
\
rs1056892
\
rs1786814
\
rs2232228
\
rs2232228
\
rs2229774
\
rs17863783
\
rs7853758
\
rs4982753
\
rs4149178

2.	Download databases
\
ClinVar: Using FTP python code
\
[o2_username@login03 ~] sbatch ALFA_pull.sh
\
gnomAD: Using gsutils in transfer directory. Repeat for chromosomes 6, 9, 12, 14, 16, 18, 21.
\
[o2_username@transfer01 ~]$ gsutil cp gs://gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr1.vcf.bgz /home/o2_username/AF query tool location/gnomad.genomes.v3.1.1.sites.chr2_$(date +%F).vcf.bgz

3.	Run AF_query with ALFA database
\
[o2_username@login03 ~]  sbatch master_AF_query.sh ALFA_freq_2021-04-23.vcf.gz cardio_variants_2021-04-23.txt 
\
This will output the flat TSV file:
\
ALFA_AF_cardio_variants_2021-04-23.tsv

4.	Run AF_query with gnomAD genome database for chromosome 1.
\
[o2_username@login03 ~]  sbatch master_AF_query.sh gnomad.genomes.v3.1.1.sites.chr2_2021-04-23.vcf.bgz cardio_variants_2021-04-23.txt
\
This will output the flat TSV file:
\
gnomad.genomes.v3.1.1.sites.chr2_AF_cardio_variants_2021-04-23.tsv
\
Repeat this step for all other chromosomes that the variants of interest are on.

5.	Export flat files using FileZilla or other file transfer software.


## Citations:
1.	Karczewski, K.J., et al., The mutational constraint spectrum quantified from variation in 141,456 humans. Nature, 2020. 581(7809): p. 434-443.
2.	L. Phan, Y.J., H. Zhang, W. Qiang, E. Shekhtman, D. Shao, D. Revoe, R. Villamarin, E. Ivanchenko, M. Kimura, Z. Y. Wang, L. Hao, N. Sharopova, M. Bihan, A. Sturcke, M. Lee, N. Popova, W. Wu, C. Bastiani, M. Ward, J. B. Holmes, V. Lyoshin, K. Kaur, E. Moyer, M. Feolo, and B. L. Kattman, "ALFA: Allele Frequency Aggregator.". National Center for Biotechnology Information, U.S. National Library of Medicine, 10 Mar. 2020.
3.	Danecek, P., et al., Twelve years of SAMtools and BCFtools. Gigascience, 2021. 10(2).

