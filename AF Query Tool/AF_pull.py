import urllib

# Pull GR38 newest VCF ALFA file
urllib.urlretrieve('ftp://ftp.ncbi.nih.gov/snp/population_frequency/latest_release/freq.vcf.gz', 'ALFA_freq.vcf.gz')
