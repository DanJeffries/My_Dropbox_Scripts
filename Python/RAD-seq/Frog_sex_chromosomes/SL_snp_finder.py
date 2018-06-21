
# coding: utf-8

# In[ ]:

from __future__ import division
import vcf
import matplotlib
#matplotlib.use('Agg') ## this allows the drawing of plots in matplotlib on the cluster, which doesn't use the X-server backend. This has something to do with display (but I don't know what)
import matplotlib.pyplot as plt
import numpy as np
import os.path
import sys
import time
import gzip


# ## This script identifies sex linked markers from a VCF file using the criteria of female-male allele freq > 0.4.
# Written by D. Jeffries 09/2015 

# #### Workflow:
#     1. Filters loci that are present in the user-specified number of samples
#     2. Calculates the allele frequencies for males and females separately
#     3. Subtracts male from female frequencies and filter loci that show signs of X or Z linkage
#     4. Outputs all male and female frequencies and female-male outputs to a single file called "yourinput.vcf.all_frequencies.tsv" (where yourinput = the name and path of your vcf file). Loci identified as X or Z linked are labelled as such in this file.
#     6. Outputs all putative X or Z linked markers to separate fasta files if any are identified.
#     7. Outputs a histogram of the distribution of female-male frequencies called "yourinput.vcf.fem-male_freqs.pdf"
#     8. All suplus information is recorded to a log file, with a summary at the end of this file. 
# 

# In[35]:

## Get arguments from the command line

Usage_message = "\n## This script identifies sex linked markers from a VCF file. Written by D. Jeffries 09/2015\n\n## USAGE (args in this order):\n\n  SL_snp_finder.py <vcf> <pop_map.txt> <catalog.tags.tsv_file> <female-male_cutoff> <sample_presence_cutoff> <coverage_per_sample_cutoff> <maf_cutoff>\n\n  <vcf> -- Path to your VCF file\n  <pop_map.txt> -- Path to the pop_map.txt file that you used in Stacks\n  <catalog.tags.tsv_file> -- Path to the catalog tags file for getting sequences for SL tags\n  <female-male_cutoff> -- The threshold value for the (female frequency - Male frequency) calulations\n  <sample_presence_cutoff> -- The minimum percent (rounded up to nearest whole individual) of samples a locus has to be present in\n  <coverage_per_sample_cutoff> -- minimum coverage for a locus\n  <maf_cutoff> -- minimum minor allele frequnce\n\n  ## Note all paths should be absolute, not relative\n\n## Workflow:\n   1. Filters loci in the vcf based on user-specified sample presence, maf and coverage thresholds\n   2. Calculates the allele frequencies for males and females separately\n   3. Subtracts male from female frequencies to identify signal of sex linkage\n   4. SL loci are identified as having allele freq of > 0.95 in all of the homogametic sex, and fem-male freq score of approx 0.5 (Xlinked) or -0.5 (Zlinked).\n   ## Note that a maximum threshold for the <female-male_cutoff> is generated to create a symmetrical window around 0.5 or -0.5,\n    eg. a <female-male_cutoff> = 0.4 will allow for scores between 0.4-0.6, whereas <female-male_cutoff> = 0.45 will allow for scores between 0.45-0.55 \n\n## Outputs:\n   1. All male and female frequencies and female-male scores are outputed to a single file called 'yourinput.vcf.all_frequencies.tsv'. Loci identified as X or Z linked are labelled as such in this file.\n   2. All putative X or Z linked markers to separate fasta files if any are identified.\n   3. Histogram of the distribution of female-male frequencies called 'yourinput.vcf.fem-male_freqs.pdf'\n   4. A barplot of the average tag coverage at all loci and at all loci not filtered called 'yourinput.vcf.coverage_by_sample.pdf\n   5. A barplot of the amount of missing data in each sample called 'missing_data_by_sample.pdf'\n   6. All suplus information is recorded to a log file, with a summary at the end of this file.\n\n"

if len(sys.argv) == 1:
    sys.exit(Usage_message)

elif len(sys.argv) < 8: ## If not enough args are supplied print error message
    sys.exit("\n##Error, not enough arguments\n"+Usage_message)

elif len(sys.argv) == 8:
    myvcfpath = sys.argv[1]
    alteredvcfpath = "%s%s" %(myvcfpath, ".altered")
    popmappath = sys.argv[2]
    catalog_tags_file = sys.argv[3]
    X_or_Z_freq_threshold = float(sys.argv[4])
    sample_presence_cutoff = float(sys.argv[5])
    coverage_threshold = int(sys.argv[6])
    maf_threshold = float(sys.argv[7])
    
else:
    sys.exit("Unknown Error\n"+Usage_message)



## set the window around the freq threshold. The window automatically tightens and relaxes around 0.5 or -0.5 

lower_thresh = X_or_Z_freq_threshold
upper_thresh = 0.5 + (0.5 - X_or_Z_freq_threshold)


# First thing to do is alter the metadata in the vcf outputted by stacks 1.30. I am not sure if it is stacks or pyvcf that is wrong, but stacks encodes the data in the allele depth field as an interger, while pyvcf expects a float. Changing the metadata line in the vcf header to contain "Number=." instead of "Number=1" fixes the issue.

myvcf = open(myvcfpath, 'r').readlines()
alteredvcf = open(alteredvcfpath, 'w')

for line in myvcf:
    if "Allele Depth" not in line:
        alteredvcf.write(line)
    elif "Allele Depth" in line:
        line = '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele Depth">\n'
        alteredvcf.write(line)
alteredvcf.close()


# ### Now calculate allele frequencies for males and females at each SNP
# ####Requires:
#     1. pyvcf module installed (can use pip, remember to add to python path. This is on the cluster!
#     2. altered vcf file from above
#     3. pop_map.txt file. Same format as used for stacks. Sample names must be the same. And males and females must be denoted by M or F (case sensitive) respectively. Must be the same file as used in populations to creat the VCF. If there are additional samples in this file the allele frequencies will be wrong!

vcf_reader = vcf.Reader(open(alteredvcfpath, 'r')) ## load in altered vcf file


## Open file for all frequency and locus info for all kept genotypes
all_frequencies = open("%s%s" %(myvcfpath, ".all_frequencies.tsv"), 'w')
all_frequencies.write("#Locus_ID\tPOS\tN_MaleChroms\tN_FemChroms\tMalefreq_REF\tMalefreq_ALT\tFemalefreq_REF\tFemalefreq_ALT\tFemREF-MaleREF\tSex_linked\n")

## Open files to output X or Z linked loci to and add headers

Putative_Xlinked_makers = []

Putative_Zlinked_makers = []

## Make a log file for all auxillary info
freq_ratios_log = open("%s%s" %(myvcfpath, ".freq_ratios.log"), 'w')
freq_ratios_log.write("Script run on %s\n " % (time.strftime("%c")))

## Make list to keep frequency information in for histogram
freq_ratio_data = []

## get male/female info for samples for working out frequencies below38      3460    50      26      0.940   0.060   0.846   0.154   -0.094  NotSexLinked

pop_map = open(popmappath, 'r').readlines()

sample_dict = {}
sample_dict["females"] = []
sample_dict["males"] = []
fem_samples = 0
male_samples = 0
sample_counter = 0
low_maf_counter = 0
kept_loci = 0
sample_missing_dict = {}
sample_cov_dict = {}
sample_cov_kept_dict = {}

for sample in pop_map: ## gets list of males and females from pop_map.txt file. At the moment there are no verbose checks for samples present in pop_map but not in vcf and vice versa
    name = sample.split()[0]
    sex = sample.strip().split()[1]
    sample_missing_dict[name] = 0

    #sample_counter += 1
    if sex == "F" or sex == "f":
        sample_dict["females"].append(name)
        fem_samples += 1
    elif sex == "M" or sex == "m":
        sample_dict["males"].append(name)
        male_samples += 1

#print sample_dict        
locus_dict = {}


## Find number of loci in input file
locus_counter = 0


for record in vcf_reader:
    locus_counter += 1
print "Number of loci = %s" % (locus_counter)
    
## Write some general stats and input options to the log file

freq_ratios_log.write("\n## User specified options:\n")
freq_ratios_log.write("Input vcf: %s\n" % (myvcfpath))
freq_ratios_log.write("pop_map file used: %s\n" % (popmappath))
freq_ratios_log.write("Catalog file used: %s\n" % (catalog_tags_file))
freq_ratios_log.write("Specified female-male threshold: %s\n" % (X_or_Z_freq_threshold))
freq_ratios_log.write("Female-male thresh range: (+/-) %s-%s\n" % (lower_thresh, upper_thresh))
freq_ratios_log.write("Min percentage samples present: %s\n" % (sample_presence_cutoff))
freq_ratios_log.write("Min coverage per genotype: %s\n" % (coverage_threshold))
freq_ratios_log.write("Min maf per locus: %s\n" % (maf_threshold))
freq_ratios_log.write("Number of female samples = %s\n" % (fem_samples))
freq_ratios_log.write("Number of male samples = %s\n" % (male_samples))
freq_ratios_log.write("Number of loci = %s\n" % (locus_counter))
                      

## Some quick counters

numb_putative_Xlinked = 0
numb_putative_Zlinked = 0
low_data_loci = 0

vcf_reader = vcf.Reader(open(alteredvcfpath, 'r'))

info_rec = vcf_reader.next()

for sample in info_rec:
    name = sample.sample    
    sample_cov_dict[name] = []
    sample_cov_kept_dict[name] = []


for record in vcf_reader:
    femREF_count = 0    ## set the counters for the reference and alternative allele (encoded as 0 in the vcf)
    femALT_count = 0
    malREF_count = 0
    malALT_count = 0   
    fem_none_count = 0
    male_none_count = 0
    low_cov_samples = 0
    n_genotypes = 0
    male_genotypes = 0
    fem_genotypes = 0
    number_of_samples = len(record.samples)
    
    
    for sample in record.samples:
        
        name = sample.sample  
        sample_cov_dict[name].append(sample['DP'])
    
    
    if record.aaf[0] < maf_threshold: ## if locus has minor allele freq lower than specified threshold then skip it
        low_maf_counter += 1
        freq_ratios_log.write("\n#LOCUS_ID: %s, Locus_POS: %s\n\n" %(record.ID, record.POS)) 
        freq_ratios_log.write("Minor allele frequence of locus is lower than specified cutoff (%s)\n" % (record.aaf[0]))
        pass  
    
    elif record.aaf >= maf_threshold:
    
        for sample in record.samples: 
            
            name = sample.sample
            
            ### For each sample, if the coverage is too low, remove the genotype for that individual.
            
            
            if sample['DP'] < coverage_threshold:
                genotype = None
                low_cov_samples += 1
                #freq_ratios_log.write("Sample %s thrown out due to low coverage (%s)\n" % (name, sample['DP']))
                
            elif sample['DP'] >= coverage_threshold:
                genotype = sample['GT']
                sample_cov_kept_dict[name].append(sample['DP'])
            ## Now calculate the female and male frequencies separately
            
            if name in sample_dict["females"]: 
                #print "Locus=", record.ID, "Female=", sample.sample, "depth=", sample['DP'], "Orig_GT=", sample['GT'], "assignedGT=", genotype                    
                #print "FEMALE Locus=", record.ID, "n_genotypes", n_genotypes
                if genotype == None: ## if no genotype exists
                    sample_missing_dict[name] += 1
                    pass
                elif genotype == "0/0":
                    femREF_count += 2
                    n_genotypes +=1
                    fem_genotypes +=1
                elif genotype == "0/1":
                    femREF_count += 1
                    femALT_count += 1
                    n_genotypes +=1
                    fem_genotypes +=1
                elif genotype == "1/0":
                    femREF_count += 1
                    femALT_count += 1
                    n_genotypes +=1
                    fem_genotypes +=1
                elif genotype == "1/1":
                    femALT_count += 2
                    n_genotypes +=1
                    fem_genotypes +=1
                #print "\tN_REF=", femREF_count, "N_ALT=", femALT_count, "Nnone=", fem_none_count
            elif name in sample_dict["males"]:
                #print "Locus=", record.ID, "Male=", sample.sample, "depth=", sample['DP'], "Orig_GT=", sample['GT'], "assignedGT=", genotype                    
                #print "MALE Locus=", record.ID, "n_genotypes", n_genotypes
                if genotype == None: ## if no genotype exists
                    sample_missing_dict[name] += 1
                    pass
                elif genotype == "0/0":
                    malREF_count += 2
                    n_genotypes +=1
                    male_genotypes +=1
                elif genotype == "0/1":
                    malREF_count += 1
                    malALT_count += 1
                    n_genotypes +=1
                    male_genotypes +=1
                elif genotype == "1/0":
                    malREF_count += 1
                    malALT_count += 1
                    n_genotypes +=1
                    male_genotypes +=1
                elif genotype == "1/1":
                    malALT_count += 2
                    n_genotypes +=1
                    male_genotypes +=1
            #else:
                #print "\n##SAMPLE NAME NOT IN POP_MAP.TXT: Sample = %s" % (name)
                #print "\tN_REF=", malREF_count, "N_ALT=", malALT_count, "Nnone=", male_none_count
                
        ## Filter loci that have too many missing samples, including samples thrown out due to low coverage!
        samples_at_locus = n_genotypes
        chromosomes_at_locus = n_genotypes*2
        percent_samples_present = n_genotypes/number_of_samples
        #print "Locus", record.ID
        #print samples_at_locus
        #print percent_samples_present
        
        if percent_samples_present >= sample_presence_cutoff:
            kept_loci += 1
            ## Calculate frequencies
            
            femREF_freq = femREF_count/(fem_genotypes*2)
            femALT_freq = femALT_count/(fem_genotypes*2)
            
            maleREF_freq = malREF_count/(male_genotypes*2)
            maleALT_freq = malALT_count/(male_genotypes*2)
                                
                        
            ## Output female stats
            freq_ratios_log.write("\n#LOCUS_ID: %s\n\n" %(record.ID))   
            freq_ratios_log.write("Number of female genotypes for this locus = %s\n" %(fem_genotypes))
            freq_ratios_log.write("Female reference count = %s\n" % (femREF_count))
            freq_ratios_log.write("Female alternative count = %s\n" % (femALT_count))
            freq_ratios_log.write("Female reference frequency = %.3f\n" % (femREF_freq))
            freq_ratios_log.write("Female alternative frequency = %.3f\n" % (femALT_freq))
                                  
            ## check fem freqs
            if not (femREF_freq) + (femALT_freq) == 1:
                freq_ratios_log.write("\n******ERROR, summed frequencies do not add up to 1******\n")
            elif (femREF_freq) + (femALT_freq) == 1:
                freq_ratios_log.write("Summed female ref and alt frequencies OK! (= %.3f)\n" % (femREF_freq + femALT_freq))
        
            ## Output male stats
            freq_ratios_log.write("Number of male genotypes for this locus = %s\n" %(male_genotypes))
            freq_ratios_log.write("Male reference count = %s\n" % (malREF_count))
            freq_ratios_log.write("Male alternative count = %s\n" % (malALT_count))
            freq_ratios_log.write("Male reference frequency = %.3f\n" % (maleREF_freq))
            freq_ratios_log.write("Male alternative frequency = %.3f\n" % (maleALT_freq))
                                  
            ## check male freqs
            if not (maleREF_freq) + (maleALT_freq) == 1:
                freq_ratios_log.write("\n******ERROR, summed frequencies do not add up to 1******\n")
            elif (maleREF_freq) + (maleALT_freq) == 1:
                freq_ratios_log.write("Summed female ref and alt frequencies OK! (= %.3f)\n" % (maleREF_freq + maleALT_freq))
        
            locus_dict[record.ID] = {}
            locus_dict[record.ID]["female_freqs"] = ["%.3f" % (femREF_freq), "%.3f" % (femALT_freq)]
            locus_dict[record.ID]["male_freqs"] = ["%.3f" % (maleREF_freq), "%.3f" % (maleALT_freq)]
            

            ####### ==============================================================================================
            ### So now I have the allele frequencies for males and females, I can subtract them and see if the distribution fits the expectition for X linked or Z linked!
            
            freq_ratio = femREF_freq - maleREF_freq
            #print freq_ratio
            freq_ratio_data.append(freq_ratio)
            

            ## Write files for X or Z linked loci
            
            if freq_ratio >= lower_thresh and freq_ratio <= upper_thresh and femREF_freq >= 0.95:  ## for X linked
                linked_status = "Xlinked"
                freq_ratios_log.write("Locus %s DOES FIT X linked criteria <------------------------\n" % (record.ID))
                Putative_Xlinked_makers.append("%s" % (record.ID))
                numb_putative_Xlinked += 1
            
            elif freq_ratio >= -upper_thresh and freq_ratio <= -lower_thresh and maleREF_freq >= 0.95:  ## for Z linked
                linked_status = "Zlinked"
                freq_ratios_log.write("Locus %s DOES FIT Z linked criteria <------------------------\n" % (record.ID))
                Putative_Zlinked_makers.append("%s" % (record.ID))
                numb_putative_Zlinked += 1
            else:
                freq_ratios_log.write("Locus %s does not fit X or Z linked criteria\n" % (record.ID))
                linked_status ="NotSexLinked"
        
            ## Write the main info file for male and female frequencies, ratios etc.
            
            all_frequencies.write("%s\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n" % (record.ID, record.POS, male_genotypes*2, fem_genotypes*2, maleREF_freq, maleALT_freq, femREF_freq, femALT_freq, freq_ratio, linked_status))
            
        
        elif percent_samples_present < sample_presence_cutoff:  ## If not enough samples at a locus then log it but don't use for female-male calculations
            freq_ratios_log.write("\n#LOCUS_ID: %s, Locus_POS: %s\n\n" %(record.ID, record.POS)) 
            freq_ratios_log.write("Number of samples at locus is lower than sample presence cutoff (%s)\n" % (samples_at_locus))
            low_data_loci += 1
    
## Look at coverage...


    
    
    
    
    
    
## Print some quick summary stats - also written at the end of the log file
print "Number of samples =", number_of_samples
print "Number of loci with too few samples = %s" % (low_data_loci)
print "Number of loci with low MAF = %s" % (low_maf_counter)
print "Number of loci with enough data = %s" % (kept_loci)
print "Number of putative X linked snps = %s" % (numb_putative_Xlinked)
print "Number of putative X linked tags = %s" % (len(set(Putative_Xlinked_makers)))
print "Number of putative Z linked markers = %s" % (numb_putative_Zlinked)
print "Number of putative Z linked tags = %s" % (len(set(Putative_Zlinked_makers)))

freq_ratios_log.write("\nSUMMARY....\n\n")
freq_ratios_log.write("Number of loci with too few samples = %s\n" % (low_data_loci))
freq_ratios_log.write("Number of loci with enough data = %s\n" % (kept_loci))
freq_ratios_log.write("Number of loci with low MAF = %s\n" % (low_maf_counter))
freq_ratios_log.write("Number of putative X linked snps = %s\n" % (numb_putative_Xlinked))
freq_ratios_log.write("Number of putative X linked tags = %s\n" % (len(set(Putative_Xlinked_makers))))
freq_ratios_log.write("Number of putative Z linked snps = %s\n" % (numb_putative_Zlinked))
freq_ratios_log.write("Number of putative Z linked tags = %s\n" % (len(set(Putative_Zlinked_makers))))


## plot histogram of frequency ratios

plt.hist(freq_ratio_data, bins = 100, color = '0.5')
plt.xlabel('Female_REFfreq - Male_REFfreq', fontsize=8)
plt.ylabel("Number of SNPs", fontsize=8)
plt.axvline(x = np.mean(freq_ratio_data), color='r', linestyle='dashed', linewidth = 2)
plt.axvline(x = X_or_Z_freq_threshold, color='b', linestyle='dashed', linewidth = 1)
plt.axvline(x = -X_or_Z_freq_threshold, color='b', linestyle='dashed', linewidth = 1)
plt.savefig("%s%s" %(myvcfpath, ".fem-male_freqs.pdf"), format = 'pdf')
plt.close()
#plt.show()


## plot histogram of missing data per sample

plt.bar(range(len(sample_missing_dict)), sample_missing_dict.values(), align='center', color = '0.5')
plt.xticks(range(len(sample_missing_dict)), sample_missing_dict.keys(), rotation = 90, fontsize=8)
plt.ylabel("Number of loci with missing data", fontsize=8)
plt.savefig("%s%s" %(myvcfpath, ".missing_data_by_sample.pdf"), format = 'pdf')
plt.close()
#plt.show()

## Plot coverage per sample

n_groups = len(sample_cov_dict)

Xlabs = []
all_means = []
all_std = []

kept_means = []
kept_std = []

for sample in sorted(sample_cov_dict.keys()):
    Xlabs.append(sample)
    all_means.append(np.mean(sample_cov_dict[sample]))
    all_std.append(np.std(sample_cov_dict[sample]))
for sample in sample_cov_kept_dict.keys():
    kept_means.append(np.mean(sample_cov_kept_dict[sample]))
    kept_std.append(np.std(sample_cov_kept_dict[sample]))

fig, ax = plt.subplots()

index = np.arange(n_groups)
bar_width = 0.35

opacity = 0.4
error_config = {'ecolor': '0.3'}

rects1 = plt.bar(index, all_means, bar_width,
                 alpha=opacity,
                 color='0.5',
                 #yerr=all_std,
                 error_kw=error_config,
                 label='All')

rects2 = plt.bar(index + bar_width, kept_means, bar_width,
                 alpha=opacity,
                 color='0',
                 #yerr=kept_std,
                 error_kw=error_config,
                 label='Kept')

plt.xlabel('Sample')
plt.ylabel('Mean coverage')
plt.title('Average coverage per sample')
plt.xticks(index + bar_width, (Xlabs), rotation = 90, size = 8)
plt.legend()

plt.tight_layout()
plt.savefig("%s%s" %(myvcfpath, ".coverage_by_sample.pdf"), format = 'pdf')




## Write fasta files of putative X or Z linked loci if there are any

Putative_Xlinked_makers = set(Putative_Xlinked_makers)
Putative_Zlinked_makers = set(Putative_Zlinked_makers)

if catalog_tags_file.endswith("gz"):
    catalog = gzip.open(catalog_tags_file, 'r').readlines()
else:
    catalog = open(catalog_tags_file, 'r').readlines()

if numb_putative_Xlinked > 0:
    Putative_Xlinked_makers_file = open("%s%s" %(myvcfpath, ".Putative_Xlinked_makers.fa"), 'w')
    
    for locus in Putative_Xlinked_makers:
        for tag in catalog:
            if locus == tag.split()[2]:
                Putative_Xlinked_makers_file.write(">X_linkedLocusID_%s\n" % (locus))
                Putative_Xlinked_makers_file.write("%s\n" % (tag.split()[8]))
    Putative_Xlinked_makers_file.close()

if numb_putative_Zlinked > 0:
    Putative_Zlinked_makers_file = open("%s%s" %(myvcfpath, ".Putative_Zlinked_makers.fa"), 'w')
    
    for locus in Putative_Zlinked_makers:
        for tag in catalog:
            if locus == tag.split()[2]:
                Putative_Zlinked_makers_file.write(">Z_linked|LocusID_%s\n" % (locus))
                Putative_Zlinked_makers_file.write("%s\n" % (tag.split()[8]))
    Putative_Zlinked_makers_file.close()


## close all unclosed files 
freq_ratios_log.close()
all_frequencies.close()

            
           

print "\n***DONE!***\n"





