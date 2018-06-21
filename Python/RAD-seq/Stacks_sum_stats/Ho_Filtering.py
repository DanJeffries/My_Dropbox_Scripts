
# coding: utf-8

# In[31]:

import sys

working_dir=sys.argv[1]
sumstats = open(working_dir+"batch_1.sumstats.tsv", 'r').readlines()
populations = open(working_dir+"uniq_pops.txt", 'r').readlines()



# In[2]:

print "populations", len(populations)


# In[33]:

## make a dictionary of keys being the pop ID and the values being the locus ID of loci with higher than 0.6 Ho
black_list = {}
for mykey in range(1,len(populations)+1):
    black_list[mykey] = []
    line_number = 1
    for line in sumstats:
        if '#' not in line and line.split()[5] == str(mykey) and float(line.split("\t")[10]) > 0.5:
            #print line_number,line.split()[1],mykey, line.split()[5], line.split("\t")[10]
            black_list[mykey].append(line.split()[1])
        line_number +=1


# In[34]:

## ok, now for all loci - if the loci is present in lists of more than one pop, then get rid of it!

## make lists of uniq locus names (so you dont count loci with more than one SNP twice)

def f1(seq):
    set = {}
    map(set.__setitem__, seq, [])
    return set.keys()

locus_list = []
for i,v in black_list.items():
    for loc in f1(v):
        locus_list.append(loc)

## Count the number of populations a locus comes up in        

from collections import Counter

locus_counts = Counter(locus_list)
#print locus_counts

## Write any loci that have Ho of higher than 0.5 in more than 2 populations to a file

blacklist_file = open(working_dir+"blacklist.txt", 'w')
whitelist_file = open(working_dir+"whitelist.txt", 'w')

loc_count = 0
for i,v in locus_counts.items():
    if locus_counts[i] > 1:
        loc_count += 1
        blacklist_file.write(i+"\n")
    elif locus_counts[i] == 0:
        whitelist_file.write(i+"\n")
        

blacklist_file.close()
whitelist_file.close()
print "Number of blacklisted loci = ", loc_count



# ###Note that the sumstats file contains all snps, where as the vcf contains only one per tag if i use the --write_single_snp flag. So the number of snps removed from the vcf is very likely to be lower than the number in the blacklist!
