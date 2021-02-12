# This script is used to pull ClinVar and/or ALFA databases directly onto your linux-based platform.
## Note: Check the links in the README file to ensure that the names/locations of these files have not changed. 
import urllib

# ClinVar database- GR37 or GR38. 
urllib.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz', 'clinvar_GR37.vcf.gz')

# ALFA database
urllib.urlretrieve('ftp://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/freq.vcf.gz', 'ALFA_freq.vcf.gz')

