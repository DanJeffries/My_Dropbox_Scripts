
# coding: utf-8

# ### This script finds Y, or W specific tags that are present only in the heterogametic sex
# 
# #### Inputs:
#     1. catalog.tags file path
#     2. Ustacks outs file path
#     3. pop_map file path (This file should include sex information for each sample).

# In[26]:

import sys
import os
import gzip
import linecache
import pprint as pp
import math


# In[25]:

# Error and usage

Usage_message = "\n##USAGE (args in this order):\nSL_snp_finder.py <path/to/catalog.tags.tsv> <path/to/Ustacks_outputs> <path/to/pop_map.txt> <per_sex_presence_cutoff>\n\n##All paths should be absolute, not relative\n\n"

if len(sys.argv) == 1:
    sys.exit(Usage_message)

elif len(sys.argv) < 5: ## If not enough args are supplied print error message
    sys.exit("\n##Error, not enough arguments\n"+Usage_message)


# ### Cline args for troubleshooting
# 
# catalog_tags_file = "/home/djeffrie/Data/Pperezi/Stacks_outs/batch_1.catalog.tags.tsv.gz"
# U_outs_path = "/home/djeffrie/Data/Pperezi/Stacks_outs/"
# popmappath = "/home/djeffrie/Data/Pperezi/Stacks_outs/pop_map_kept.txt"
# sex_presence_thresh = 0.5
# 
# log = open("%s/%s" % (catalog_tags_file.rpartition("/")[0], "Presence_absence.log"), 'w')     
# log.write("Input parameters:\n\nCatalog file: %s\nUstacks ouptuts: %s\nPop_map_file: %s\nMinimum percentage of M of F present: %1.2f\n" % (catalog_tags_file, U_outs_path, popmappath, sex_presence_thresh))
# 

# In[47]:

catalog_tags_file = sys.argv[1]
U_outs_path = sys.argv[2]
popmappath = sys.argv[3]
sex_presence_thresh = float(sys.argv[4])

## open log file and write params
log = open("%s/%s" % (catalog_tags_file.rpartition("/")[0], "Presence_absence.log"), 'w')     
log.write("Input parameters:\n\nCatalog file: %s\nUstacks ouptuts: %s\nPop_map_file: %s\nMinimum percentage of M of F present: %1.2f\n\n" % (catalog_tags_file, U_outs_path, popmappath, sex_presence_thresh))


# In[72]:

## First get sample ID and sex info for each sample
## Include only samples in the pop_map file

sex_file = open(popmappath, 'r').readlines()

kept_sample_names = []

for line in sex_file:
    name = line.split()[0]
    kept_sample_names.append(name)

sample_dict = {}
sex_dict = {}
sex_dict["Males"] = []
sex_dict["Females"] = []

male_count = 0
female_count = 0

for root, dirs, files in os.walk(U_outs_path):
    for fil in files:
        sex = None
        if "tags.tsv" in fil and "catalog" not in fil:
            sample_name = fil.split(".tags")[0]
            
            if sample_name in kept_sample_names:
                
                sample_dict[sample_name] = []
            
                if fil.endswith("gz"):
                    tagsfile = gzip.open("%s/%s" % (root, fil), 'r')
                    tagsfile.readline() ## pass over first line
                    ID = tagsfile.readline().split()[1]
                    sample_dict[sample_name].append(ID)
                else:
                    tagsfile = open("%s/%s" % (root, fil), 'r')
                    tagsfile.readline() ## pass over first line
                    ID = tagsfile.readline().split()[1]
                    sample_dict[sample_name].append(ID)
                    
                for sample in sex_file:
                    if sample_name == sample.split()[0]:
                        sex = sample.split()[1]
                
                if sex == "M":
                    sample_dict[sample_name].append(sex)
                    sex_dict["Males"].append(ID)
                    male_count += 1
                elif sex == "F":
                    sample_dict[sample_name].append(sex)
                    sex_dict["Females"].append(ID)
                    female_count += 1
            else:
                log.write("Sample %s in stacks outputs but not in pop_map file\n" % (sample_name))
               
log.write("\n\n############################\n\nLOG:\n\n# 'BAD' tag status means that one or more samples had two or more tags at the same \ncatalog locus, thus, this locus was not used in analyses\n\n")
log.write("Tag_status\tTag_ID\tN_males\tN_females\tY_linked\tW_linked\n") ## headers for log info


# ####Now look through each tag in catalog file:
#     find tags which are present in sex_presence_thresh males and no females (Y linked)
#     find tags which are present in sex_presence_thresh females and no males (W linked)

# In[73]:

if catalog_tags_file.endswith("gz"):
    cat_file = gzip.open(catalog_tags_file, 'r').readlines()
else:
    cat_file = open(catalog_tags_file, 'r').readlines()

    


# In[74]:

n_males_required = math.ceil(male_count*sex_presence_thresh) ## math.ceil rounds up to nearest whole individual
n_females_required = math.ceil(female_count*sex_presence_thresh)

putative_Ylinked_tags = []
putative_Wlinked_tags = []


for tag in cat_file[1:]:
    
    ## First get all the info I need
    
    N_males_at_locus = 0 
    N_females_at_locus = 0
    duplicates_present = "no" ## Checking for replicate samples (Often a tag can be present in the same individual twice, need to remove these)
    samples_at_locus = [] ## for checking for replicate samples
    
    tag_ID =  tag.split()[2]
    samples_field = tag.split()[7] ## the field in the catalog file that contains the samples present at that locus
    
    Y_linked_str = "NO" ## defaults for log file
    W_linked_str = "NO"
    tag_status = "BAD"
    
    for sample in samples_field.split(","):
        
        sample_ID = sample.split("_")[0]
        
        if sample_ID not in samples_at_locus: ## check that there are no replicates
            if sample_ID in sex_dict["Males"]:
                N_males_at_locus += 1
                samples_at_locus.append(sample_ID)
            elif sample_ID in sex_dict["Females"]:
                N_females_at_locus += 1
                samples_at_locus.append(sample_ID)
        
        elif sample_ID in samples_at_locus:
            duplicates_present = "yes"
            
    if duplicates_present == "no": ## If there are no sample duplicates at the tag
        
        tag_status = "OK"
        
        ## Look for Y linked tags
        
        if N_females_at_locus == 0 and N_males_at_locus > n_males_required:
            seq = tag.split()[8]
            putative_Ylinked_tags.append((tag_ID, N_females_at_locus, N_males_at_locus, seq))
            Y_linked_str = "YES"
            #print "Y-LINKED: Tag_ID=", tag_ID, "N_males=", N_males_at_locus, "N_females=", N_females_at_locus, seq
        
        ## look for W linked tags
        
        elif N_males_at_locus == 0 and N_females_at_locus > n_females_required:
            seq = tag.split()[8]
            putative_Wlinked_tags.append((tag_ID, N_females_at_locus, N_males_at_locus, seq))
            W_linked_str = "YES"
            #print "W-LINKED: Tag_ID=", tag_ID, "N_males=", N_males_at_locus, "N_females=", N_females_at_locus, seq
            
    log.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (tag_status, tag_ID, N_males_at_locus, N_females_at_locus, Y_linked_str, W_linked_str))
    
              
if len(putative_Ylinked_tags) > 0:
    Yfa = open("%s/%s" % (catalog_tags_file.rpartition("/")[0], "Y_linked_tags.fa"), 'w')
    for line in putative_Ylinked_tags:
        Yfa.write(">Y_linked|LocusID_%s|N_females_%s|N_males_%s\n%s\n" % (line[0], line[1], line[2], line[3]))
    Yfa.close()
              
if len(putative_Wlinked_tags) > 0:
    Wfa = open("%s/%s" % (catalog_tags_file.rpartition("/")[0], "W_linked_tags.fa"), 'w')
    for line in putative_Wlinked_tags:
        Wfa.write(">W_linked|LocusID_%s|N_females_%s|N_males_%s\n%s\n" % (line[0], line[1], line[2], line[3]))
    Wfa.close()
    
## Write log file summary

log.write("\nSUMMARY:\nNumber of females: %s\n" % (female_count))
log.write("Number of males: %s\n" % (male_count))
log.write("Number of Putative Y linked tags: %s\n" % (len(putative_Ylinked_tags)))
log.write("Number of Putative W linked tags: %s\n" % (len(putative_Wlinked_tags)))

## Print summary to STDOUT

print "\nSUMMARY:\nNumber of females: %s" % (female_count)
print "Number of males: %s" % (male_count)
print "Number of Putative Y linked tags: %s" % (len(putative_Ylinked_tags))
print "Number of Putative W linked tags: %s" % (len(putative_Wlinked_tags))
              
log.close()
print "\n ### DONE! ###\n"
    


# In[ ]:



