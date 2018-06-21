
# coding: utf-8

# In[11]:

import sys

## Usage Pairwise_pop_map_maker [master_pop_map_path] [output path to put new pop_map] [pop1] [pop2]

# master_pop_map = open("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/pop_codes.txt", 'r').readlines() ## For testing

master_pop_map = open(sys.argv[1], 'r').readlines()
output_dir= sys.argv[2]
pop1 = "BF" #sys.argv[3]
pop2 = "V" #sys.argv[4]

PW_pop_map = open(output_dir+"/"+"pop_codes.txt",'w')

for line in master_pop_map:
    if pop1 == line.strip().split()[1]:
        PW_pop_map.write(line)
    elif pop2 == line.strip().split()[1]:
        PW_pop_map.write(line)


PW_pop_map.close()

