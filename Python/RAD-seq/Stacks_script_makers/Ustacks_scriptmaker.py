
# coding: utf-8

# ###Make scripts for ustacks.  
# ###To use, put all fastq files in the same directory. 
# #### Usage: Ustacks_scriptmaker.py \<inputdir>  \<outputdir>  \<M value> \<m value>"
# ####\<inputdir> = path to fastq files
# ####\<outputdir> = path you want to write Ustacks outputs to
# #### You can specify M and m values. All other params are default. To alter, manually edit line 23 of this script.

# In[5]:

import sys
import os

if len(sys.argv) < 5:
    print "\n***********************\nNOT ENOUGH ARGUMENTS\nUsage:\nUstacks_scriptmaker.py <inputdir>  <outputdir>  <M value> <m value>\n********************\n"

indir='%s/' %(sys.argv[1])
outdir='%s/' %(sys.argv[2])
Mval = '%s' %(sys.argv[3])
mval = '%s' %(sys.argv[4])

def Ustacks_script_maker(parent_dir, output_dir):  ## Parent_dir = directory containing all files for Ustacking
    samples = []
    Ustacks_commands = []
    for root, dirs, files in os.walk(parent_dir):
        for fil in files:
            if fil.endswith("fq") or fil.endswith("fq_1") or fil.endswith("fq_2") or fil.endswith("fq.gz") or fil.endswith("fq_1.gz") or fil.endswith("fq_2.gz") and 'discards' not in fil:
                samples.append(root+fil)
                samples = sorted(samples)
    
    ID = 1
    for sample in samples:
        Ustacks_commands.append("ustacks -t fastq -f "+sample+" -d -r -i "+str(ID)+" -M "+Mval+" -m "+mval+" -o "+output_dir+" -p 7")
        ID += 1
        
    f = open(output_dir+'Ustacks_scripts.sh', 'w')
    f.write('#!/bin/bash\n')
    for command in Ustacks_commands:
        f.write(command+'\n')
    f.close()    
        
    return Ustacks_commands
                
U_commands = Ustacks_script_maker(indir, outdir)

