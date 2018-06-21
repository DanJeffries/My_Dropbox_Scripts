import vcf
import numpy as np
import sys

def sample_Avcov_from_vcf(vcf_path, output_file):

    import vcf
    import numpy as np
    
    myvcf = vcf.Reader(open(vcf_path, 'r'))
    
    sample_cov = {}
    
    for record in myvcf:
	for sample in record.samples:
		if not sample['GT'] == "./.":  ## disregard missing loci!!
        		if sample.sample not in sample_cov.keys():
                	    sample_cov[sample.sample] = []
	                else:
        	            sample_cov[sample.sample].append(sample['DP'])
		else:
	                pass     
    
    
    sample_avgs = {}
    for sample in sample_cov.keys():
        sample_avgs[sample] = np.round(np.mean(sample_cov[sample]), 2)
    
    outfile = open(output_file, 'w')
        
    for sample in sample_avgs:
        outfile.write("%s\t%s\n" %(sample.rpartition("_")[0], sample_avgs[sample]))
                      
    outfile.close()                 



sample_Avcov_from_vcf(sys.argv[1],sys.argv[2])

