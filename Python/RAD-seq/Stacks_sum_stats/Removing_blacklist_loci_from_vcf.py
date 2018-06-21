import sys
# coding: utf-8

# In[ ]:

working_dir = sys.argv[1]

blacklist = open(working_dir+"/blacklist.txt", 'r').readlines()
vcf = open(working_dir+"/batch_1.vcf", 'r').readlines()
filtered_vcf = open(working_dir+"/Ho_filtered.vcf", 'w')

blacklist_stripped = [i.strip() for i in blacklist]

for line in vcf:
    if line.startswith("#"):
        filtered_vcf.write(line)

for snp in vcf:
    if '#' not in snp and snp.strip().split()[2] not in blacklist_stripped:
        filtered_vcf.write(snp)


filtered_vcf.close()

