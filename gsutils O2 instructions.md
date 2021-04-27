# Gsutils instructions for O2 04/23/2020

## Downlaod gnomAD exome file:
\
The easiest way to download this directly to O2 is via gsutils. For this example, we will show how to download Chromosome 14. See O2 internal wiki page  for more information (https://wiki.rc.hms.harvard.edu/display/O2/File+Transfer).
\
\
Command line: This will prompt you to provide password and DUO authentication to switch over to transfer directory.
\
[o2_username@login03 ~] ssh o2_username@transfer.rc.hms.harvard.edu
\
\
Command line: Once successfully in the transfer directory, copy gnomAD exomes Chr14 into your merge tool directory location. 
\
[o2_username@transfer01 ~]$ gsutil cp gs://gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.14.vcf.bgz /home/o2_username/merge tool location/gnomad.exomes.r2.1.1.sites.14_$(date +%F).vcf.bgz
\
\
*Note: if you want to explore gnomAD options, in command line:
\
[o2_username@transfer01  ~]$ gsutil ls gs://gcp-public-data--gnomad/
\
\
Once transfer is complete, switch back over to regular O2 by either opening up new O2 session or submit the command below. This will prompt you to provide password and DUO authentication to switch over to default directory. 
\
[o2_username@transfer01 ~]$ ssh o2_username@login.rc.hms.harvard.edu
