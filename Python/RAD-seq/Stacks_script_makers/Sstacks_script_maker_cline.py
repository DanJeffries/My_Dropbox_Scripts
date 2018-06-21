
# coding: utf-8

# In[2]:

import sys
import os

### Check and get arguments
if len(sys.argv) == 1:
    sys.exit("\nUsage:\nSstacks_scriptmaker.py <inputdir> <catalog_path> <batch_id> <threads>\n")  

if len(sys.argv) < 5:
    sys.exit("\nNOT ENOUGH ARGUMENTS\nUsage:\nSstacks_scriptmaker.py <inputdir> <catalog_path> <batch_id> <threads>\n")  

input_dir = sys.argv[1]
catalog = sys.argv[2]
batch_id = sys.argv[3]
threads = sys.argv[4]

## Make script
sstacks_script = open("%s/%s" % (input_dir, "Sstacks_scripts.sh"), 'w')
sstacks_script.write("#!/bin/bash/\n")

for root, dirs, files in os.walk(input_dir):
    for fil in files:
        if "tags.tsv" in fil and "catalog" not in fil:
            infile = "%s/%s" % (root, fil.partition(".")[0])
#	    infile = "%s/%s.fil"  % (root, fil.partition(".")[0]) ## for the crucian file names	

            command = "sstacks -b %s -c %s -s %s -o %s -p %s;\n" % (batch_id, catalog, infile, input_dir, threads)
            sstacks_script.write(command)
            print command.strip()
sstacks_script.close()

