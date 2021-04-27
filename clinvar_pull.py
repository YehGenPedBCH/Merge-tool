import urllib

# Pull GR37 clinvar dataset
urllib.urlretrieve('ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz', 'clinvar_GR37.vcf.gz')
