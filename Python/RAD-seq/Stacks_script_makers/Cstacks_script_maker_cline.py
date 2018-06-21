
# coding: utf-8

# ##Make scripts for Cstacks.  
# ###To use, put all fastq files in the same directory. 
# ### Usage: Ustacks_scriptmaker.py \<inputdir>  \<outputdir>  \<n value> \<p value>"
# ###\<inputdir> = path to fastq files
# ###\<outputdir> = path you want to write Ustacks outputs to
# ### You can specify N and p values. All other params are default. To alter, manually edit script
# ### Note, this script needs a pop_codes file like the one needed in populaitons, format:
# 
#     sample_1    pop_1
#     sample_2    pop_1
#     sample_3    pop_2
#     sample_4    pop_2
#     ....

# ## Make scripts

# In[8]:

import os
import re
import sys

## Check usage 

if len(sys.argv) < 6:
    sys.exit("\nNOT ENOUGH ARGUMENTS\nUsage<\nCstacks_scriptmaker.py <inputdir>  <outputdir> <pop_file_path> <n value> <p value>\n")

indir='%s/' %(sys.argv[1])
outdir='%s/' %(sys.argv[2])
pop_file='%s' %(sys.argv[3])
nval = '%s' %(sys.argv[4])
pval = '%s' %(sys.argv[5])

## define function

def Cstacks_script_maker(input_dir, output_dir, pop_file_path, n_value, p_value):
    samples = []
    pop_list = []
    samples_ready = []
    cwd = input_dir
    
    pops = open(pop_file_path, 'r').readlines()
    for pop in pops:
        pop_list.append(pop.strip().split()[1])
    print 'Populations for catalog = ', set(pop_list) , 'N populations = '+str(len(set(pop_list)))
    
    for root, dirs, files in os.walk(input_dir):
        for fil in files:
            #print fil
            for population in pop_list: 
                if 'tags.tsv' in fil:
                    if fil.endswith('tsv'):
                        if population in fil:
                           samples.append(fil.rpartition('.')[0].rpartition('.')[0])
                    elif fil.endswith('gz'):
                           samples.append(fil.rpartition('.')[0].rpartition('.')[0].rpartition('.')[0]) 
    samples = sorted(set(samples))
    
    for sample in samples:
        samples_ready.append(' -s '+input_dir+sample)
    samples = ''.join(samples_ready)
    
    Ccommand = 'cstacks -b 1 -n'+n_value+' -p '+p_value+' '+samples+' -o '+output_dir
    
    f = open(output_dir+'Cstacks_scripts.sh','w')
    f.write('#!/bin/bash\n##working directory = '+str(cwd)+'\n')
    for command in Ccommand:
        f.write(command)
    f.close()
    
    return Ccommand

## run function using command line args

Cstacks_script_maker(indir, outdir, pop_file, nval, pval)
  


# In[ ]:



