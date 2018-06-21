import sys
# coding: utf-8

# In[ ]:

vcf_path= sys.argv[1]
blacklist = sys.argv[2]


blacklist = open(blacklist, 'r').readlines()
vcf = open(vcf_path, 'r').readlines()
filtered_vcf = open(vcf_path.rpartition(".")[0]+"_pure_cru_blacklisted.vcf", 'w')

blacklist_stripped = [i.strip() for i in blacklist]

for line in vcf:
    if line.startswith("#"):
        filtered_vcf.write(line)

for snp in vcf:
    if '#' not in snp and snp.strip().split()[2] not in blacklist_stripped:
        filtered_vcf.write(snp)


filtered_vcf.close()

