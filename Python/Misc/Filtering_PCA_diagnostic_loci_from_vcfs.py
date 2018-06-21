
# coding: utf-8

# ### Filtering Ho_filtered vcfs for the Diagnostic SNP markers ID'd using PC loadings

# In[13]:

import sys

wd = sys.argv[1]

whitelist = open(str(wd+"/Diagnostic_snps.txt"), 'r').readlines()
original_vcf = open(str(wd+"/Ho_filtered_altered.vcf"), 'r').readlines()
original_sumstats = open(str(wd+"/batch_1.sumstats.tsv"), 'r').readlines()
Diagnostic_SNPS = open(str(wd+"/Diagnostic_SNPS.vcf"), 'w')
Diagnostic_SNPS_sumstats = open(str(wd+"/Diagnostic_SNPS_sumstats.tsv"), 'w')

whitelist = [i.strip() for i in whitelist]
for line in original_vcf:
    if line.startswith("#"):
        Diagnostic_SNPS.write(line)
    elif line.split()[2] in whitelist:
        Diagnostic_SNPS.write(line)

sumstats_lines = [] ## for the loci that have more than one snp - put the loc ID in here and if there is more than one snp it will only keep the first - i.e. the same snp as in the vcf
for line in original_sumstats:
    if line.startswith("#"):
        Diagnostic_SNPS_sumstats.write(line)
    elif line.split()[1] in whitelist and line.split()[1] not in sumstats_lines:
        sumstats_lines.append(line.split()[1])
        Diagnostic_SNPS_sumstats.write(line)
        
Diagnostic_SNPS_sumstats.close()
Diagnostic_SNPS.close()

