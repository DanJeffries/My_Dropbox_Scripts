
# coding: utf-8

# In[ ]:

## label changer
import sys

data_file = open(sys.argv[1], 'r').readlines()
outfile = open(sys.argv[2], 'w')
label_keys = open(sys.argv[3], 'r').readlines()

for datline in data_file:
	for line in label_keys:
		old = line.split()[0]
		new = line.strip().split()[1]
		if old in datline:
	            datline = datline.replace(old, new)
	        else:
			pass
	outfile.write(datline)
outfile.close()



