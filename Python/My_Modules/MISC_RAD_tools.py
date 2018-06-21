from __future__ import division

def fasta_maka(whitey, cat, out = None):

    """
    whitey = whitelist (either a python list or a file path) containing locus IDs in the form of "<Tag_ID>_<Position>"
    cat    = path to the catalog file to get sequences from

    """

    import sys
    import gzip
    
    if isinstance(whitey, str):
        loci = open(whitey, 'r').readlines()
    elif isinstance(whitey, (list, set)):
        loci = whitey
    else:
        sys.exit("Unknown whitelist format - expected a python list or a file path")
        
    if cat.endswith("gz"):
        tags = gzip.open(cat, 'r').readlines()
    else:
        tags = open(cat, 'r').readlines()

    ## Pull out the locus ID's from the whitelist

    Loc_IDs = []
    for locus in loci:
        if locus.startswith("compli"):
            Loc_id = locus.split("_")[1]
        else:
            Loc_id = locus.split("_")[0]
        Loc_IDs.append(Loc_id.strip())
    
    print "Number of tags in whitelist:",len(Loc_IDs)

    ## Write the fasta
    
    if not out == None:
        fasta = open(out, 'w')
        outpath = out
    else:
        fasta = open("%s/%s" % (cat.rpartition('/')[0], 'Whitelist_tags.fa'), 'w')
        outpath = "%s/%s" % (cat.rpartition('/')[0], 'Whitelist_tags.fa')

    count = 0
    for line in tags:
        if 'consensus' in line:
            Tag_ID = line.split()[2]
            if Tag_ID in Loc_IDs:
                count+=1
                fasta.write('>'+ Tag_ID +'\n'+line.split()[8]+'\n')
                
    print count, "sequences written to", outpath

    fasta.close()

def filter_vcf_tag_ID_only(vcfpath, tags_to_keep, outfile_name):

    """
    Filters a vcf ("vcfpath"), keeping only the loci that are in the <list> "tags_to_keep"
    """

    import vcf

    XY_set = set(tags_to_keep)

    out_vcf = open("%s/%s" % (vcfpath.rpartition("/")[0], outfile_name), 'w')

    vcf_parsed = vcf.Reader(open(vcfpath, 'r'))

    XY_ids = []

    for line in open(vcfpath,'r').readlines():
        if line.startswith("#"):
            out_vcf.write(line)
        else:
            loc_id = "%s" % (line.split()[2])
            if loc_id in XY_set:
                out_vcf.write(line)
    out_vcf.close()



def filter_vcf(vcfpath, tags_to_keep, outfile_name):
    
    """
    Filters a vcf ("vcfpath"), keeping only the loci that are in the <list> "tags_to_keep"
    """

    import vcf

    XY_set = set(tags_to_keep)
        
    out_vcf = open("%s/%s" % (vcfpath.rpartition("/")[0], outfile_name), 'w')
    
    vcf_parsed = vcf.Reader(open(vcfpath, 'r'))

    XY_ids = []

    for line in open(vcfpath,'r').readlines():
        if line.startswith("#"):
            out_vcf.write(line)
        else:
	    loc_id = "%s_%s" % (line.split()[2], line.split()[1])
            if loc_id in XY_set:
                out_vcf.write(line)
    out_vcf.close()



def missing_data_finder(input_file, output_file_path, blacklisted_tags_cutoff=100):
    
    '''
    
    This function will count the number of tags that droupout in all samples
    in the specified VCF. It will output the data and simple bar plots of both 
    the number of dropouts per sample, and the number of dropouts per tag. It will also 
    output a file called blacklist.txt containing locus IDs which can be used to remove 
    these loci from subsequent runs of the stacks populations module.
    
        input_file              -  full path to vcf
        output_file_path        -  path to output data, plots and blacklist to
        blacklisted_tags_cutoff -  number of 'worst' tags to blacklist
        
    '''
    
    
    import sys
    import vcf 
    import numpy as np
    import matplotlib.pyplot as plt 
    from collections import Counter
    from operator import itemgetter


    ## Alter the allele depth header line - which doesn't seem to be compatible between Stacks and PyVCF. 
    
    myvcf = open(input_file, 'r').readlines()
    alteredvcfpath = "%s%s" %(input_file, ".altered")
    alteredvcf = open(alteredvcfpath, 'w')
    
    for line in myvcf:
        if "Allele Depth" not in line:
            alteredvcf.write(line)
        elif "Allele Depth" in line:
            line = '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele Depth">\n'
            alteredvcf.write(line)
    alteredvcf.close()
    
    ## Now open altered file and continue
    
    alteredvcf = open(alteredvcfpath, 'r')
    
    vcf_reader = vcf.Reader(alteredvcf) 
    record = next(vcf_reader)
    
    samples = vcf_reader.samples ## get all the samples in the file

    ## Count the missing data --------------------------------------------------------------------------------------------------------
    
    missing_data_count = {}
    missing_data_tags = {}
    
    for record in vcf.Reader(open(alteredvcfpath, 'r')): 
        for sample in samples: ## for each sample
            if sample not in missing_data_count:
                missing_data_count[sample] = 0
            
            if record.ID not in missing_data_tags:
                missing_data_tags[record.ID] = 0
            
            if record.genotype(sample)['GT'] == None: ## If the current samples genotype for this record is missing ('./.'), count it with no_genotype_counter
                missing_data_count[sample] += 1
                missing_data_tags[record.ID] += 1

    
    ######### TAG DROPOUTS PER SAMPLE ######## -------------------------------------------------------------------------------------------
    
    sorted_tags_per_sample = sorted(missing_data_count.items(), key=itemgetter(1), reverse=True) ## write file from this sorted tags list
    
    ## Print tag dropouts per sample to file
    f = open(output_file_path+'/Tag_dropouts_per_sample.txt', 'w')
    f.write('sample\t#_missing_loci\n')
    for sample in sorted_tags_per_sample:
        f.write("%s\t%s\n" % (sorted_tags_per_sample[0], sorted_tags_per_sample[1]))
    f.close()
    

    
    ######### TAG DROPOUTS PER SAMPLE ######## -------------------------------------------------------------------------------------------
    
    x = []
    y = []
    xlabs = []
    counter = 1

    plt.figure(figsize=(70,50))

    for i, j in sorted_tags_per_sample:
        xlabs.append(i)        
        x.append(counter)
        y.append(int(j))
        counter += 1

    ## Plot the per-tag droupouts
    plt.bar(x,y, edgecolor = "none",align="center")
    plt.yticks(range(0,max(y),100), fontsize = 30)
    plt.xticks(x, xlabs, fontsize = 30, rotation = "vertical")
    plt.ylabel("Number of dropouts", fontsize = 30)
    plt.xlabel("Sample", fontsize = 30)
    plt.axhline(np.mean(y), linewidth=4, color='r')
    plt.title("Tag dropouts per sample", fontsize = 40)

    plt.tight_layout()
    plt.savefig("%s/Tag_dropouts_per_sample.pdf" % output_file_path)
    plt.show()
    plt.close()
    
    
    ######### TAG DROPOUTS PER TAG ######## -------------------------------------------------------------------------------------------
    
    sorted_tags = sorted(missing_data_tags.items(), key=itemgetter(1), reverse=True)
    
    x = []
    y = []
    xlabs = []
    blacklisted = []

    cut_off = blacklisted_tags_cutoff
    counter = 1
    plt.figure(figsize=(70,50))

    for i, j in sorted_tags:
        xlabs.append(i)
        x.append(counter)
        y.append(int(j))
        if counter <= blacklisted_tags_cutoff:
            blacklisted.append(j) ## write this to file
        counter += 1

    xarray = np.asarray(x)
    yarray = np.asarray(y)

    mask1 = xarray < cut_off
    mask2 = xarray >= cut_off

    ## write per-tag droupout data to file

    tag_drop_outs = open("%s/Tag_dropouts_per_tag.txt" % output_file_path, 'w')
    tag_drop_outs.write("TagID\tNumber_of_dropouts\n")

    for tag in sorted_tags:
        tag_drop_outs.write("%s\t%s\n" % (tag[0], tag[1]))
    tag_drop_outs.close()
    
    ## Write the blacklist file
    
    blacklist = open(output_file_path+'/blacklist.txt', 'w')
    for tag_ID in set(blacklisted):
        blacklist.write("%s\n" % (tag_ID))
    blacklist.close()

    ## Plot the per-tag droupouts
    plt.bar(xarray[mask1],yarray[mask1], edgecolor = "none")
    plt.bar(xarray[mask2],yarray[mask2], color = "grey", edgecolor = "none")
    plt.yticks(range(0,max(y),50), fontsize = 30)
    plt.xticks(range(0,max(x),100), fontsize = 30, rotation = "vertical")
    plt.ylabel("Number of dropouts", fontsize = 30)
    plt.xlabel("Tag", fontsize = 30)
    plt.axhline(np.mean(yarray), linewidth=4, color='r')
    plt.title("Tag dropouts per tag", fontsize = 40)

    #plt.tight_layout()
    
    plt.savefig("%s/Tag_dropouts_per_tag.pdf" % output_file_path)
    #plt.show()
    plt.close()
    
    print 'Outputs written to:'
    print "%s/Tag_dropouts_per_sample.txt" % output_file_path
    print "%s/Tag_dropouts_per_tag.txt" % output_file_path
    print "%s/Tag_dropouts_per_sample.pdf" % output_file_path
    print "%s/Tag_dropouts_per_tag.pdf" % output_file_path
    print "%s/blacklist.txt" % output_file_path
    
    return missing_data_count


def missing_data_per_pop(missing_data_per_sample_dict, pop_file_path, outpath):
    
    '''
    This function takes the outputs from the missing_data_finder function
    along with the pop_map file used in stacks and creates an additional 
    dataset and barplot for the tag dropouts per population.
    
        missing_data_per_sample_dict     - Outputs from missing_data_finder
        pop_file_path                    - Full path to the pop map text file used in stacks (popualtions)
        
    '''
    from matplotlib import pyplot as plt
    import numpy as np
    
    pop_map = open(pop_file_path, 'r').readlines()
    
    pop_missing_dict = {}
    
    for line in pop_map:
        sample_name = line.strip().split()[0]
        pop_name = line.strip().split()[1]
        
        if sample_name in missing_data_per_sample_dict:
            
            if pop_name not in pop_missing_dict:
                pop_missing_dict[pop_name] = [int(missing_data_per_sample_dict[sample_name])]
            else:
                pop_missing_dict[pop_name].append(int(missing_data_per_sample_dict[sample_name]))
                
    x = pop_missing_dict.keys()
    y = [sum (i) for i in pop_missing_dict.values()]
    z = zip(x,y)            
    ## write to file
    
    pop_dropouts = open("%s/Tag_dropouts_per_population.txt" % outpath, 'w')
    
    for pop in z:
        pop_dropouts.write("%s\t%s\n" % (pop[0], pop[1]))
    pop_dropouts.close()
    
    ## Plot population tag dropouts
    plt.figure(figsize=(70,50))
    plt.bar(range(len(x)), y, align='center', width = 0.6)
    plt.xticks(range(len(x)), x, fontsize = 30)
    plt.yticks(range(0,max(y),(np.round(max(y)/10))),range(0,max(y),(np.round(max(y)/10))) , fontsize = 30)
    plt.ylabel("Number of tag dropouts", fontsize = 30)
    plt.title("Tag dropouts per population", fontsize = 40)
    plt.savefig("%s/Tag_dropouts_per_population.pdf" % outpath)
    #plt.show()
    plt.close()


def Cstacks_script_maker(input_dir, output_dir, pop_file_path, nval, pval):
    
    '''
    
    Makes Cstacks scripts and saves them to .sh file for running.
    
    input_dir        -  full path to the ustacks output files
    output_dir       -  path to write the catalog to
    pop_file_path    -  full path to the pop_map file
    nval             -  value for n (mismatch)
    pval             -  value for p (threads)
    
    '''
    
    import os
    import re
    import sys
    
    samples = []
    pop_list = []
    samples_ready = []
    
    for root, dirs, files in os.walk(input_dir):
        for fil in files:
            #print fil
            if 'tags.tsv' in fil:
                if fil.endswith('tsv'):
                    samples.append(fil.rpartition('.')[0].rpartition('.')[0])
                elif fil.endswith('tsv.gz'):
                        samples.append(fil.rpartition('.')[0].rpartition('.')[0].rpartition('.')[0]) 
    samples = sorted(set(samples))
    
    for sample in samples:
        samples_ready.append(' -s '+input_dir+sample)
    samples = ''.join(samples_ready)
    
    Ccommand = 'cstacks -b 1 -n%s -p %s %s -o %s'  % (nval,pval,samples,output_dir)
    
    f = open("%s/Cstacks_scripts.sh" % output_dir,'w')
    f.write('#!/bin/bash\n\n')
    
    
    for command in Ccommand:
        f.write(command)
    f.close()
    
    return Ccommand


def Sstacks_script_maker(input_dir,threads, batch_id = 1):
    
    '''
    Makes Sstacks scripts
    input_dir
    threads
    '''
    
    import os
    
    sstacks_commands = []
    
    ## Make script
    sstacks_script = open("%s/%s" % (input_dir, "Sstacks_scripts.sh"), 'w')
    sstacks_script.write("#!/bin/bash/\n")

    ## First pass to get catalog file
    for root, dirs, files in os.walk(input_dir):
        for fil in files:
            if "tags.tsv" in fil and "catalog" in fil:
                catalog = "%s/%s" % (root, fil.split(".catalog")[0])
    
    ## Second pass to make script
    for root, dirs, files in os.walk(input_dir):
        for fil in files:
            if "tags.tsv" in fil and "catalog" not in fil:
                infile = "%s/%s" % (root, fil.partition(".")[0])

                command = "sstacks -b %s -c %s -s %s -o %s -p %s" % (batch_id, catalog, infile, input_dir, threads)
                script_command = "sstacks -b %s -c %s -s %s -o %s -p %s;\n" % (batch_id, catalog, infile, input_dir, threads)
                sstacks_script.write(script_command)
                sstacks_commands.append(command)
                
                
    sstacks_script.close()
    return sstacks_commands




def av_tag_cov(vcf_path, out_path):
    """
    Makes a histogram for tag coverage from a vcf.

    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    vcf = open(vcf_path, 'r')

    tag_means = []
    tag_counter = 0
    for line in vcf.readlines():
        tag_cov = []
        if not line.startswith("#"):
    
            tag_counter += 1
            for sample in line.split()[9:]:
                tag_cov.append(float(sample.split(":")[1]))
                
            tag_cov = np.asarray(tag_cov)
            tag_means.append(np.mean(tag_cov))
    
    ## Print some stats
    print "Number of tags = ", tag_counter
    print "Average tag coverage", np.mean(tag_means), "(+-", np.std(tag_means), ")"
    
     
    ## plot the histogram - can save this if I want to make it cline executable
    plt.hist(tag_means, bins = 100)
    plt.show()
    
    ## write to a file in the output directory specified
    nf = open(out_path, 'w')

    for tag in tag_means:
        if isinstance(tag, (int, float)):
            nf.write(str(tag)+"\n")
            
    nf.close()



def Super_av_tag_cov(vcf_path, whitelist = None, popmap = None):
    """
     
     
    This function alters behaviour depending on what you give it:
    
    If you give it just a VCF, it will make a simple histogram of 
    average tag coverage across all samples. 
    
    If you give it a whitelist, it will do the same, just for the tags IDs
    specified in the whitelist. This whitelist can be a file path (string) or 
    a list of tag IDs passed from within python. 
    
    If you give it a popmap (same format as Stacks), then it will give you
    the same histogram, as well as a box plot of average tag coverages for each 
    population. 

    """
    import numpy as np
    import matplotlib.pyplot as plt
    import vcf as VCF
    
    ## get whitelist if it exists
       
    if not whitelist == None:
        if isinstance(whitelist, str):
            whitelisted_loci = open(whitelist, 'r').readlines()
        
        elif isinstance(whitelist, list):
            whitelisted_loci = whitelist
    
    else:
        whitelisted_loci = []
    
    ## alter the header in the vcf if not done already
    
    if not vcf_path.endswith("altered"):
        myvcf = open(vcf_path, 'r').readlines()
        alteredvcfpath = "%s%s" %(vcf_path, ".altered")
        alteredvcf = open(alteredvcfpath, 'w')

        for line in myvcf:
            if "Allele Depth" not in line:
                alteredvcf.write(line)
            elif "Allele Depth" in line:
                line = '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele Depth">\n'
                alteredvcf.write(line)
        alteredvcf.close()
        vcf = VCF.Reader(open(alteredvcfpath, 'r'))
    else:
        vcf = VCF.Reader(open(vcf_path, 'r'))
    
    ## get popmap if it exists and make dict for the data

    if not popmap == None:
        
        pop_info = open(popmap, 'r').readlines()
        
        pop_sample_cov_dict = {}
        sample_pop_dict = {}
        
        for line in pop_info:
                    
            line = line.strip()

            sample_name = line.split()[0]
            pop_name = line.split()[1]
            
            if sample_name not in sample_pop_dict:  ## for easy access to the pop info for a sample when looping over vcf below
                sample_pop_dict[sample_name] = pop_name
            
            if pop_name not in pop_sample_cov_dict: ## all data goes here for popmap plots
                pop_sample_cov_dict[pop_name] = {}
            pop_sample_cov_dict[pop_name][sample_name] = []
        
    #print pop_sample_cov_dict
    
    
    
    ## Get coverages ---------------------------------------------------------------------------------------------------------------


    tag_means = []
    tag_counter = 0

    
    for record in vcf:
        
        if whitelist == None or record.ID in whitelisted_loci: ## process a locus if it is in the whitelist or if no whitelist was specified.
            
            tag_cov = []
            tag_counter += 1
            
            for sample in record.samples:
                tag_cov.append(sample["DP"]) ## for normal plot with no pop or sample info
                
                if not popmap == None: ## If a popmap was specified
                    pop_sample_cov_dict[sample_pop_dict[sample.sample]][sample.sample].append(sample["DP"])                            
        
        
            tag_means.append(np.mean(tag_cov))

    if not popmap == None:
        
        pop_sample_dict_means = {}

        for sex in pop_sample_cov_dict:
            if sex not in pop_sample_dict_means:
                pop_sample_dict_means[sex] = {}

                for sample in pop_sample_cov_dict[sex]:
                    if len(pop_sample_cov_dict[sex][sample]) > 0:
                        pop_sample_dict_means[sex][sample] = np.mean(pop_sample_cov_dict[sex][sample])
                    else:
                        pop_sample_dict_means[sex][sample] = 0
            
   
    
    ## Print some stats
    print "Number of tags = ", tag_counter
    print "Average tag coverage", np.mean(tag_means), "(+-", np.std(tag_means), ")"

    
    ## plot the histogram - can save this if I want to make it cline executable
    plt.figure(figsize = (10,10))
    plt.title("Histogram of average coverage per SNP locus")
    plt.xlabel("Coverage per SNP")
    plt.hist(tag_means, bins = 100, edgecolor = "white", color = "royalblue")
    plt.show()

    ### box plot for pop_map addition!!

    ## Great Box plot tutorial here for changing features of plot http://blog.bharatbhole.com/creating-boxplots-with-matplotlib/ 
    
    if not popmap == None:
        
        fig = plt.figure(figsize = (10,10))
        ax = plt.subplot(111)

        position = 1
        labels = []
        for sex in pop_sample_dict_means:
            bp = ax.boxplot([pop_sample_dict_means[sex][i] for i in pop_sample_dict_means[sex]], positions = [position], widths = 0.75, patch_artist=True)

            for box in bp['boxes']:
                # change outline color
                box.set( color='dodgerblue', linewidth=2)
                # change fill color
                box.set( facecolor = 'dodgerblue' )

            for median in bp['medians']:
                median.set(color='white', linewidth=2)

            for whisker in bp['whiskers']:
                whisker.set(color='dodgerblue', linewidth=2)

            for cap in bp['caps']:
                cap.set(color='dodgerblue', linewidth=2)

            position += 1
            labels.append(sex)

        ax.set_xlim(0, len(pop_sample_dict_means)+1)
        ax.set_xticks(range(1, len(pop_sample_dict_means)+1))
        ax.set_xticklabels(labels)
        plt.show()
        
        
    ## write to a file in the same dir as the vcf
    out_path = "%s/Av_tag_covs.txt" % vcf_path.rpartition("/")[0]
    
    nf = open(out_path, 'w')

    for tag in tag_means:
        if isinstance(tag, (int, float)):
            nf.write(str(tag)+"\n")

    nf.close()
    
    if not popmap == None:
        return tag_means, pop_sample_cov_dict
    else:
        return tag_means



def sample_Avcov_from_vcf(vcf_path):

    import vcf
    import numpy as np
    from matplotlib import pyplot as plt
    
    
    myvcf = open(vcf_path, 'r').readlines()
    alteredvcfpath = "%s%s" %(vcf_path, ".altered")
    alteredvcf = open(alteredvcfpath, 'w')
    
    for line in myvcf:
        if "Allele Depth" not in line:
            alteredvcf.write(line)
        elif "Allele Depth" in line:
            line = '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele Depth">\n'
            alteredvcf.write(line)
    alteredvcf.close()
    
    
    myvcf = vcf.Reader(open(alteredvcfpath, 'r'))
    
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
    
    fig = plt.figure(figsize = (15,10))
    
    
    sample_avgs = {}
    x = []
    xlabs = []
    y = []
    sample_count = 0
    for sample in sample_cov.keys():
        sample_count += 1
        average = np.round(np.mean(sample_cov[sample]), 2)
        x.append(sample_count)
        xlabs.append(sample)
        y.append(average)
        sample_avgs[sample] = average
        
    print "Mean sample coverage = %s (+/%s) " % (np.round(np.mean(sample_avgs.values()),2), np.round(np.std(sample_avgs.values()),2))
        
    plt.bar(x,y, color = "royalblue", edgecolor = "none")
    plt.ylabel("Snp Coverage")
    plt.xticks(x, xlabs, rotation = "vertical")
    plt.show()

    return sample_cov
                  


### Sex linked marker finding scripts =================================================================================================


def SL_snp_freq(myvcfpath, popmappath, catalog_tags_file, X_or_Z_freq_threshold = 0.4, sample_presence_cutoff = 0.75, coverage_threshold = 3, maf_threshold = 0.05, homgam_REF_freq = 0.95):    
   
    """
    This script identifies sex linked markers from a VCF file using the criteria of female-male allele freq 

    Arguments:
    
    myvcfpath                     - path to vcf file (note this will be altered to make header compatible with Pyvcf. New vcf will have same name 
                                    with ".altered" appended to the end)
    
    popmappath                    - path to population map file containing sex information. Same format as Stacks pop map file.
    
    catalog_tags_file             - The catalog tags file used to create the vcf
    
    X_or_Z_freq_threshold         - The lower threshold for the freq caluclation to find sex linked snps, e.g. for an XY system, a threshold
                                    of 0.4 means that f(F) - f(M) can be >= 0.4 and <= 0.6 (the upper threshold is automatically calculated to be
                                    the same distance above 0.5 as the lower threshold is below 0.5)
    
    sample_presence_cutoff        - a locus must be called in at least this proportion of all samples (not within populations) to be considered
    
    coverage_threshold            - a locus must have at least this threshold in a sample to be considered for that sample. Note that loci below this 
                                    threshold will be removed from a sample, and this can push the locus below the sample presence cut-off, which will
                                    then remove the locus.
    
    maf_threshold                 - minor allele frequency cutoff for a locus across all samples. 
    
    
    Workflow:
    1. Filters loci that are present in the user-specified number of samples
    2. Calculates the allele frequencies for males and females separately
    3. Subtracts male from female frequencies and filter loci that show signs of X or Z linkage
    4. Outputs all male and female frequencies and female-male outputs to a single file called 
       "yourinput.vcf.all_frequencies.tsv" (where yourinput = the name and path of your vcf file). Loci 
       identified as X or Z linked are labelled as such in this file.
    5. Outputs all putative X or Z linked markers to separate fasta files if any are identified.
    6. Outputs a histogram of the distribution of female-male frequencies called "yourinput.vcf.fem-male_freqs.pdf"
    7. All suplus information is recorded to a log file, with a summary at the end of this file.
    
    
    """
    
    import vcf
    import matplotlib
    #matplotlib.use('Agg') ## this allows the drawing of plots in matplotlib on the cluster, which doesn't use the X-server backend. This has something to do with display (but I don't know what)
    import matplotlib.pyplot as plt
    import numpy as np
    import os.path
    import sys
    import time
    import gzip
    
    
    print "\n##### Using SNP frequency approach #### \n"
    
    ## set the window around the freq threshold. The window automatically tightens and relaxes around 0.5 or -0.5 

    lower_thresh = X_or_Z_freq_threshold
    upper_thresh = 1 - X_or_Z_freq_threshold


    # First thing to do is alter the metadata in the vcf outputted by stacks 1.30. I am not sure if it is stacks or pyvcf that is wrong, but stacks encodes the data in the allele depth field as an interger, while pyvcf expects a float. Changing the metadata line in the vcf header to contain "Number=." instead of "Number=1" fixes the issue.

    myvcf = open(myvcfpath, 'r').readlines()
    
    alteredvcfpath = "%s%s" %(myvcfpath, ".altered")
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
    freq_ratios_log = []
    #open("%s%s" %(myvcfpath, ".freq_ratios.log"), 'w')
    freq_ratios_log.append("Script run on %s\n " % (time.strftime("%c")))

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

    for sample in pop_map:
        name = sample.split()[0]
        sex = sample.strip().split()[1]
        sample_missing_dict[name] = 0

        #sample_counter += 1
        if sex == "F" or sex == "f" or sex == "Female" or sex == "female" or sex == "Fem" or sex == "fem":
            sample_dict["females"].append(name)
            fem_samples += 1
        elif sex == "M" or sex == "m" or sex == "Male" or sex == "male" or sex == "Mal" or sex == "mal":
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

    freq_ratios_log.append("\n## User specified options:\n")
    freq_ratios_log.append("Input vcf: %s\n" % (myvcfpath))
    freq_ratios_log.append("pop_map file used: %s\n" % (popmappath))
    freq_ratios_log.append("Catalog file used: %s\n" % (catalog_tags_file))
    freq_ratios_log.append("Specified female-male threshold: %s\n" % (X_or_Z_freq_threshold))
    freq_ratios_log.append("Female-male thresh range: (+/-) %s-%s\n" % (lower_thresh, upper_thresh))
    freq_ratios_log.append("Min percentage samples present: %s\n" % (sample_presence_cutoff))
    freq_ratios_log.append("Min coverage per genotype: %s\n" % (coverage_threshold))
    freq_ratios_log.append("Min maf per locus: %s\n" % (maf_threshold))
    freq_ratios_log.append("Number of female samples = %s\n" % (fem_samples))
    freq_ratios_log.append("Number of male samples = %s\n" % (male_samples))
    freq_ratios_log.append("Number of loci = %s\n" % (locus_counter))


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
	loc_id = "%s_%s" % (record.ID, record.POS)

        for sample in record.samples:

            name = sample.sample  
            sample_cov_dict[name].append(sample['DP'])


        if record.aaf[0] < maf_threshold: ## if locus has minor allele freq lower than specified threshold then skip it
            low_maf_counter += 1
            freq_ratios_log.append("\n#LOCUS_ID: %s, Locus_POS: %s\n\n" %(record.ID, record.POS)) 
            freq_ratios_log.append("Minor allele frequence of locus is lower than specified cutoff (%s)\n" % (record.aaf[0]))
            pass  

        elif record.aaf >= maf_threshold:

            for sample in record.samples: 

                name = sample.sample

                ### For each sample, if the coverage is too low, remove the genotype for that individual.


                if sample['DP'] < coverage_threshold:
                    genotype = None
                    low_cov_samples += 1
                    

                elif sample['DP'] >= coverage_threshold:
                    genotype = sample['GT']
                    sample_cov_kept_dict[name].append(sample['DP'])
                ## Now calculate the female and male frequencies separately

                if name in sample_dict["females"]: 

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
                freq_ratios_log.append("\n#LOCUS_ID: %s\n\n" %(loc_id))   
                freq_ratios_log.append("Number of female genotypes for this locus = %s\n" %(fem_genotypes))
                freq_ratios_log.append("Female reference count = %s\n" % (femREF_count))
                freq_ratios_log.append("Female alternative count = %s\n" % (femALT_count))
                freq_ratios_log.append("Female reference frequency = %.3f\n" % (femREF_freq))
                freq_ratios_log.append("Female alternative frequency = %.3f\n" % (femALT_freq))

                ## check fem freqs
                if not (femREF_freq) + (femALT_freq) == 1:
                    freq_ratios_log.append("\n******ERROR, summed frequencies do not add up to 1******\n")
                elif (femREF_freq) + (femALT_freq) == 1:
                    freq_ratios_log.append("Summed female ref and alt frequencies OK! (= %.3f)\n" % (femREF_freq + femALT_freq))

                ## Output male stats
                freq_ratios_log.append("Number of male genotypes for this locus = %s\n" %(male_genotypes))
                freq_ratios_log.append("Male reference count = %s\n" % (malREF_count))
                freq_ratios_log.append("Male alternative count = %s\n" % (malALT_count))
                freq_ratios_log.append("Male reference frequency = %.3f\n" % (maleREF_freq))
                freq_ratios_log.append("Male alternative frequency = %.3f\n" % (maleALT_freq))

                ## check male freqs
                if not (maleREF_freq) + (maleALT_freq) == 1:
                    freq_ratios_log.append("\n******ERROR, summed frequencies do not add up to 1******\n")
                elif (maleREF_freq) + (maleALT_freq) == 1:
                    freq_ratios_log.append("Summed female ref and alt frequencies OK! (= %.3f)\n" % (maleREF_freq + maleALT_freq))

                locus_dict[record.ID] = {}
                locus_dict[record.ID]["female_freqs"] = ["%.3f" % (femREF_freq), "%.3f" % (femALT_freq)]
                locus_dict[record.ID]["male_freqs"] = ["%.3f" % (maleREF_freq), "%.3f" % (maleALT_freq)]


                ####### ==============================================================================================
                ### So now I have the allele frequencies for males and females, I can subtract them and see if the distribution fits the expectition for X linked or Z linked!

                freq_ratio = femREF_freq - maleREF_freq
                #print freq_ratio
                freq_ratio_data.append(freq_ratio)


                ## Write files for X or Z linked loci

                if freq_ratio >= lower_thresh and freq_ratio <= upper_thresh and femREF_freq >= homgam_REF_freq:  ## for X linked
                    linked_status = "Xlinked"
                    freq_ratios_log.append("Locus %s DOES FIT X linked criteria <------------------------\n" % (loc_id))
                    Putative_Xlinked_makers.append("%s" % (loc_id))
                    numb_putative_Xlinked += 1

                elif freq_ratio >= -upper_thresh and freq_ratio <= -lower_thresh and maleREF_freq >= homgam_REF_freq:  ## for Z linked
                    linked_status = "Zlinked"
                    freq_ratios_log.append("Locus %s DOES FIT Z linked criteria <------------------------\n" % (loc_id))
                    Putative_Zlinked_makers.append("%s" % (loc_id))
                    numb_putative_Zlinked += 1
                else:
                    freq_ratios_log.append("Locus %s does not fit X or Z linked criteria\n" % (record.ID))
                    linked_status ="NotSexLinked"

                ## Write the main info file for male and female frequencies, ratios etc.

                all_frequencies.write("%s\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n" % (record.ID, record.POS, male_genotypes*2, fem_genotypes*2, maleREF_freq, maleALT_freq, femREF_freq, femALT_freq, freq_ratio, linked_status))


            elif percent_samples_present < sample_presence_cutoff:  ## If not enough samples at a locus then log it but don't use for female-male calculations
                freq_ratios_log.append("\n#LOCUS_ID: %s, Locus_POS: %s\n\n" %(record.ID, record.POS)) 
                freq_ratios_log.append("Number of samples at locus is lower than sample presence cutoff (%s)\n" % (samples_at_locus))
                low_data_loci += 1


    ## Print some quick summary stats - also written at the end of the log file
    print "Number of samples =", number_of_samples
    print "Number of loci with too few samples = %s" % (low_data_loci)
    print "Number of loci with low MAF = %s" % (low_maf_counter)
    print "Number of loci with enough data = %s" % (kept_loci)
    print "Number of putative X linked snps = %s" % (numb_putative_Xlinked)
    print "Number of putative X linked tags = %s" % (len(set(Putative_Xlinked_makers)))
    print "Number of putative Z linked markers = %s" % (numb_putative_Zlinked)
    print "Number of putative Z linked tags = %s" % (len(set(Putative_Zlinked_makers)))

    freq_ratios_log.append("\nSUMMARY....\n\n")
    freq_ratios_log.append("Number of loci with too few samples = %s\n" % (low_data_loci))
    freq_ratios_log.append("Number of loci with enough data = %s\n" % (kept_loci))
    freq_ratios_log.append("Number of loci with low MAF = %s\n" % (low_maf_counter))
    freq_ratios_log.append("Number of putative X linked snps = %s\n" % (numb_putative_Xlinked))
    freq_ratios_log.append("Number of putative X linked tags = %s\n" % (len(set(Putative_Xlinked_makers))))
    freq_ratios_log.append("Number of putative Z linked snps = %s\n" % (numb_putative_Zlinked))
    freq_ratios_log.append("Number of putative Z linked tags = %s\n" % (len(set(Putative_Zlinked_makers))))


    ## plot histogram of frequency ratios

    plt.hist(freq_ratio_data, bins = 100, color = '0.5')
    plt.xlabel('Female_REFfreq - Male_REFfreq', fontsize=8)
    plt.ylabel("Number of SNPs", fontsize=8)
    plt.axvline(x = np.mean(freq_ratio_data), color='r', linestyle='dashed', linewidth = 2)
    plt.axvline(x = X_or_Z_freq_threshold, color='b', linestyle='dashed', linewidth = 1)
    plt.axvline(x = -X_or_Z_freq_threshold, color='b', linestyle='dashed', linewidth = 1)
    plt.savefig("%s%s" %(myvcfpath, ".fem-male_freqs.pdf"), format = 'pdf')
    plt.show()



    all_frequencies.close()
    

    print "\n***DONE!***\n"

    return set(Putative_Xlinked_makers), set(Putative_Zlinked_makers), freq_ratios_log



def SL_snp_het(myvcfpath, popmappath, catalog_tags_file, homogamtic_homozygosity_threshold = 0.9, heterogametic_heterozygosity_threshold = 0.5, sample_presence_cutoff = 0.75, coverage_threshold = 3, maf_threshold = 0.05):
    
   
    """
    This script identifies sex linked markers from a VCF file using the criteria of female-male allele freq 

    Arguments:
    
    myvcfpath                            - path to vcf file (note this will be altered to make header compatible with Pyvcf. 
                                           New vcf will have same name with ".altered" appended to the end)
    
    popmappath                           - path to population map file containing sex information. Same format as Stacks pop map file.
    
    catalog_tags_file                    - The catalog tags file used to create the vcf
    
    homogamtic_homozygosity_threshold    - The lower threshold for the proportion of homozygotes in the homogametic sex at a locus
    
    heterogamtic_heterozygosity_threshold- The lower threshold for the proportion of heterozygotes in the heterogametic sex at a locus
        
    sample_presence_cutoff               - a locus must be called in at least this proportion of all samples (not within populations) to be considered
    
    coverage_threshold                   - a locus must have at least this threshold in a sample to be considered for that sample. Note that loci below this 
                                           threshold will be removed from a sample, and this can push the locus below the sample presence cut-off, which will
                                           then remove the locus.
    
    maf_threshold                        - minor allele frequency cutoff for a locus across all samples. 
    
    
    Workflow: ### NOTE*** different to SL_snp_finder()
    
    1. Filters loci that are present in the user-specified number of samples
    2. Looks for loci that are heterozygous in all or most of the homogametic sex and homozygous in all or most of the heterogametic sex
    3. Outputs all putative X or Z linked markers to separate fasta files if any are identified.
    4. All suplus information is recorded to a log file, with a summary at the end of this file.
    
    
    """
    
    import vcf
    import matplotlib
    #matplotlib.use('Agg') ## this allows the drawing of plots in matplotlib on the cluster, which doesn't use the X-server backend. This has something to do with display (but I don't know what)
    import matplotlib.pyplot as plt
    import numpy as np
    import os.path
    import sys
    import time
    import gzip
    
    
    
    print "\n##### Using SNP heterozygosity approach #####\n "

    # First thing to do is alter the metadata in the vcf outputted by stacks 1.30. I am not sure if it is stacks or pyvcf that is wrong, but stacks encodes the data in the allele depth field as an interger, while pyvcf expects a float. Changing the metadata line in the vcf header to contain "Number=." instead of "Number=1" fixes the issue.  
    
    myvcf = open(myvcfpath, 'r').readlines()
    
    alteredvcfpath = "%s%s" %(myvcfpath, ".altered")
    alteredvcf = open(alteredvcfpath, 'w')

    for line in myvcf:
        if "Allele Depth" not in line:
            alteredvcf.write(line)
        elif "Allele Depth" in line:
            line = '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele Depth">\n'
            alteredvcf.write(line)
    alteredvcf.close()


    # ### Now calculate heterozygosity for males and females at each SNP

    vcf_reader = vcf.Reader(open(alteredvcfpath, 'r')) ## load in altered vcf file

    Putative_Xlinked_makers = []
    Putative_Zlinked_makers = []

    ## Make a log file for all auxillary info
    het_approach_log = []
    #open("%s%s" %(myvcfpath, ".het_approach.log"), 'w')
    het_approach_log.append("Script run on %s\n " % (time.strftime("%c")))

    ## get male/female info 
    
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

    for sample in pop_map:
        name = sample.split()[0]
        sex = sample.strip().split()[1]
        sample_missing_dict[name] = 0

        #sample_counter += 1
        if sex == "F" or sex == "f" or sex == "Female" or sex == "female" or sex == "Fem" or sex == "fem":
            sample_dict["females"].append(name)
            fem_samples += 1
        elif sex == "M" or sex == "m" or sex == "Male" or sex == "male" or sex == "Mal" or sex == "mal":
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

    het_approach_log.append("\n## User specified options:\n")
    het_approach_log.append("Input vcf: %s\n" % (myvcfpath))
    het_approach_log.append("pop_map file used: %s\n" % (popmappath))
    het_approach_log.append("Catalog file used: %s\n" % (catalog_tags_file))
    het_approach_log.append("homogamtic_homozygosity_threshold: %s\n" % (homogamtic_homozygosity_threshold))
    het_approach_log.append("heterogamtic_heterozygosity_threshold : %s\n" % (heterogametic_heterozygosity_threshold))
    het_approach_log.append("Min percentage samples present: %s\n" % (sample_presence_cutoff))
    het_approach_log.append("Min coverage per genotype: %s\n" % (coverage_threshold))
    het_approach_log.append("Min maf per locus: %s\n" % (maf_threshold))
    het_approach_log.append("Number of female samples = %s\n" % (fem_samples))
    het_approach_log.append("Number of male samples = %s\n" % (male_samples))
    het_approach_log.append("Number of loci = %s\n" % (locus_counter))


    ## Some quick counters

    numb_putative_Xlinked = 0
    numb_putative_Zlinked = 0
    low_data_loci = 0
    kept_loci = 0

    vcf_reader = vcf.Reader(open(alteredvcfpath, 'r'))

    info_rec = vcf_reader.next()

    for sample in info_rec:
        name = sample.sample
        sample_cov_dict[name] = []
        sample_cov_kept_dict[name] = []


    for record in vcf_reader:
        loc_id = "%s_%s" % (record.ID, record.POS)
        number_of_samples = len(record.samples) ## Number of samples (not genotypes)
        samples_present = record.num_called ## Number of samples called
        percent_samples_present = samples_present/number_of_samples ## percentage of samples genotyped at this locus
	low_cov_samples = 0        

        if percent_samples_present >= sample_presence_cutoff:
                
            if record.aaf[0] < maf_threshold: ## if locus has minor allele freq lower than specified threshold then skip it
                low_maf_counter += 1
                het_approach_log.append("\n#LOCUS_ID: %s, Locus_POS: %s\n\n" %(record.ID, record.POS)) 
                het_approach_log.append("MAF is lower than cutoff (%s), locus not used\n" % (record.aaf[0]))

            elif record.aaf >= maf_threshold:

                male_genotypes = 0
                fem_genotypes = 0

                het_males_count = 0 
                het_females_count = 0 
                hom_males_count = 0 
                hom_females_count = 0 

                for sample in record.samples:

                    if sample.called:  ## If the sample is called, keep going
                        name = sample.sample

                        if sample['DP'] < coverage_threshold: ## If the sample has enough coverage keep going
                            low_cov_samples += 1
                        elif sample['DP'] >= coverage_threshold:
                            sample_cov_kept_dict[name].append(sample['DP'])

                            if name in sample_dict["females"]: ## If the sample is female, adjust female counts accordingly
                                fem_genotypes += 1

                                if sample.is_het:
                                    het_females_count += 1
                                elif not sample.is_het:
                                    hom_females_count += 1

                            elif name in sample_dict['males']: ## or if the sample is male, adjust the male counts
                                male_genotypes += 1

                                if sample.is_het:
                                    het_males_count += 1
                                elif not sample.is_het:
                                    hom_males_count += 1
                    
		
		## Account for situations where there are no males or females
		
		if all([male_genotypes > 0, fem_genotypes > 0]):


	                ## Calculate frequencies per sex (these frequencies are based on the number of males or females called at each number, so accounts for missing data)

	                male_homozygosity = hom_males_count/male_genotypes
	                male_heterozygosity = het_males_count/male_genotypes

	                female_homozygosity = hom_females_count/fem_genotypes
	                female_heterozygosity = het_females_count/fem_genotypes


	                ## Output female stats
	                het_approach_log.append("\n#LOCUS_ID: %s\n\n" %(loc_id))   
	                het_approach_log.append("Number of female genotypes for this locus = %s\n" %(fem_genotypes))
	                het_approach_log.append("Female homozygosity = %.3f\n" % (female_homozygosity))
	                het_approach_log.append("Female heterozygosity = %.3f\n" % (female_heterozygosity))

	                ## Output male stats
	                het_approach_log.append("Number of male genotypes for this locus = %s\n" %(male_genotypes))
	                het_approach_log.append("Male homozygosity = %.3f\n" % (male_homozygosity))
	                het_approach_log.append("Male heterozygosity = %.3f\n" % (male_heterozygosity))


	                ## Find and write files for X or Z linked loci  

			N_genotypes_at_locus = male_genotypes + fem_genotypes
 
			if N_genotypes_at_locus/number_of_samples >= sample_presence_cutoff: ### Need to do a second check for enough samples at this locus after the filters!
				kept_loci += 1 ## This locus can be used!
	
		                if all([female_homozygosity >= homogamtic_homozygosity_threshold, male_heterozygosity >= heterogametic_heterozygosity_threshold]):
		                    linked_status = "Xlinked"
		                    het_approach_log.append("Locus %s DOES FIT X linked criteria <------------------------\n" % (loc_id))
		                    Putative_Xlinked_makers.append("%s" % (loc_id))
		                    numb_putative_Xlinked += 1
		
		                elif all([male_homozygosity >= homogamtic_homozygosity_threshold, female_heterozygosity >= heterogametic_heterozygosity_threshold]):  ## for Z linked
		                    linked_status = "Zlinked"
		                    het_approach_log.append("Locus %s DOES FIT Z linked criteria <------------------------\n" % (loc_id))
		                    Putative_Zlinked_makers.append("%s" % (loc_id))
		                    numb_putative_Zlinked += 1
		                else:
		                    het_approach_log.append("Locus %s does not fit X or Z linked criteria\n" % (loc_id))
		                    linked_status ="NotSexLinked"

			else:
				low_data_loci += 1
				het_approach_log.append("\n#LOCUS_ID: %s, Locus_POS: %s\n\n" %(record.ID, record.POS))
	                        het_approach_log.append("Too many sample genotypes filtered from locus, locus discarded\n")
				het_approach_log.append("%s genotypes / %s samples in total = %s, miniumum threshold = %s\n" % (N_genotypes_at_locus, number_of_samples, N_genotypes_at_locus/number_of_samples, sample_presence_cutoff))
        	else:
			het_approach_log.append("\n#LOCUS_ID: %s, Locus_POS: %s\n\n" %(record.ID, record.POS))
			het_approach_log.append("Locus is missing from one sex: Males = %s, Females = %s" %(male_genotypes, fem_genotypes))

        elif percent_samples_present < sample_presence_cutoff:  ## If not enough samples at a locus then log it but don't use for female-male calculations
            het_approach_log.append("\n#LOCUS_ID: %s, Locus_POS: %s\n\n" %(record.ID, record.POS)) 
            het_approach_log.append("Number of samples at locus (%s) is lower than sample presence cutoff (%s)\n" % (percent_samples_present, samples_present))
            low_data_loci += 1


    ## Print some quick summary stats - also written at the end of the log file
    print "Number of samples =", number_of_samples
    print "Number of loci with too few samples = %s" % (low_data_loci)
    print "Number of loci with low MAF = %s" % (low_maf_counter)
    print "Number of loci with enough data = %s" % (kept_loci)
    print "Number of putative X linked snps = %s" % (numb_putative_Xlinked)
    print "Number of putative X linked tags = %s" % (len(set(Putative_Xlinked_makers)))
    print "Number of putative Z linked markers = %s" % (numb_putative_Zlinked)
    print "Number of putative Z linked tags = %s" % (len(set(Putative_Zlinked_makers)))

    het_approach_log.append("\nSUMMARY....\n\n")
    het_approach_log.append("Number of loci with too few samples = %s\n" % (low_data_loci))
    het_approach_log.append("Number of loci with enough data = %s\n" % (kept_loci))
    het_approach_log.append("Number of loci with low MAF = %s\n" % (low_maf_counter))
    het_approach_log.append("Number of putative X linked snps = %s\n" % (numb_putative_Xlinked))
    het_approach_log.append("Number of putative X linked tags = %s\n" % (len(set(Putative_Xlinked_makers))))
    het_approach_log.append("Number of putative Z linked snps = %s\n" % (numb_putative_Zlinked))
    het_approach_log.append("Number of putative Z linked tags = %s\n" % (len(set(Putative_Zlinked_makers))))


    print "\n ### DONE! ### \n"
    
    return set(Putative_Xlinked_makers), set(Putative_Zlinked_makers), het_approach_log




def SL_tag_finder(catalog_tags_file, popmappath, sex_presence_thresh = 0.5):
    
    import gzip
    import math

    ## Note popmappath here should point to a file that contains usual two popmap columns plus a
    ## a column for the ID of each sample . . . 
    
    print "\n##### Using Sex specific tag approach ##### \n"
    
    ## open log file and write params
    log = []
    #open("%s/%s" % (catalog_tags_file.rpartition("/")[0], "Presence_absence.log"), 'w')     
    log.append("Input parameters:\n\nCatalog file: %s\nPop_map_file: %s\nMinimum percentage of M or F present: %1.2f\n\n" % (catalog_tags_file, popmappath, sex_presence_thresh))

    
    ## First get sample ID and sex info for each sample
    ## Include only samples in the pop_map file

    sex_file = open(popmappath, 'r').readlines()

    kept_sample_names = []
    
    sample_dict = {}
    sex_dict = {}
    sex_dict["Males"] = []
    sex_dict["Females"] = []
    male_count = 0
    female_count = 0

    ## Get sample IDs
    
    for line in sex_file:
        sample_name = line.split()[0]
        sex = line.split()[1]
        ID = line.split()[2]
        
        kept_sample_names.append(sample_name)
        sample_dict[sample_name] = []
        sample_dict[sample_name].append(ID)
        
        if sex == "M" or sex == "m" or sex == "Male" or sex == "male":
            sample_dict[sample_name].append(sex)
            sex_dict["Males"].append(ID)
            male_count += 1
        elif sex == "F" or sex == "f" or sex == "Female" or sex == "female":
            sample_dict[sample_name].append(sex)
            sex_dict["Females"].append(ID)
            female_count += 1

    log.append("\n\n############################\n\nLOG:\n\n# 'BAD' tag status means that one or more samples had two or more tags at the same \ncatalog locus, thus, this locus was not used in analyses\n\n")
    log.append("Tag_status\tTag_ID\tN_males\tN_females\tY_linked\tW_linked\n") ## headers for log info
    
    
    ## Now look through each tag in catalog file:
        ## find tags which are present in sex_presence_thresh males and no females (Y linked)
        ## find tags which are present in sex_presence_thresh females and no males (W linked)
    
    
    if catalog_tags_file.endswith("gz"):
        cat_file = gzip.open(catalog_tags_file, 'r').readlines()
    else:
        cat_file = open(catalog_tags_file, 'r').readlines()
    
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
                putative_Ylinked_tags.append(tag_ID)
                Y_linked_str = "YES"

            ## look for W linked tags

            elif N_males_at_locus == 0 and N_females_at_locus > n_females_required:
                putative_Wlinked_tags.append(tag_ID)
                W_linked_str = "YES"
        

        log.append("%s\t%s\t%s\t%s\t%s\t%s\n" % (tag_status, tag_ID, N_males_at_locus, N_females_at_locus, Y_linked_str, W_linked_str))


    ## Write log file summary

    log.append("\nSUMMARY:\nNumber of females: %s\n" % (female_count))
    log.append("Number of males: %s\n" % (male_count))
    log.append("Number of Putative Y linked tags: %s\n" % (len(putative_Ylinked_tags)))
    log.append("Number of Putative W linked tags: %s\n" % (len(putative_Wlinked_tags)))
    
    ## Print summary to STDOUT

    print "\nSUMMARY:\nNumber of females: %s" % (female_count)
    print "Number of males: %s" % (male_count)
    print "Number of Putative Y linked tags: %s" % (len(putative_Ylinked_tags))
    print "Number of Putative W linked tags: %s" % (len(putative_Wlinked_tags))
    print "\n ### DONE! ###\n"
    
    
    ## return list of sex linked tags
    
    return set(putative_Ylinked_tags), set(putative_Wlinked_tags), log

def SL_markers_from_sexed_fam(vcf_path, sex_info_file, catalog_path, offspring_presence_threshold, mendelian_threshold, het_sex_thresh, hom_sex_thresh):
    
    """
      
    Usage: SL_markers_from_sexed_fam  <vcf_path>  <sex_info_file>  <catalog_path> <offspring_presence_threshold>  <mendelian_threshold> <het_sex_thresh> <hom_sex_thresh>
    
    This script will identify sex linked markers using the genotype of parents and sexed offspring. Specifically, it will identify alleles specific 
    to one of the parents and then track these loci through the progeny. The script will filter loci that are not present in enough samples 
    (specified by the user) and also those that have segregation patterns indicative of incorrect genotype calles in either the parents or 
    offspring (See below)

    ### Informative locus filters:
    -Loci are used only if the number of offspring present is equal to or above the percentage specified by the user
    -Loci are used only if they are heterozygous in one sex and homozygous in another.
    
    ###Null allele filters:

    Null allele in Parents:
    -If an offspring is found to have an allele that isn't present in the parent, this locus is discarded

    Null allele in offspring:
    -If only 1 offspring is homozygous for the parental minor allele then that sample is assumed to possess a null allele and the data for that individual 
     is coded as missing. **Note, the sample presence/absence filter works after this step, so some loci may be pushed over the "missing data" limit by 
     this null allele filter. 

    -If two or more samples are found to be homozygous for the minor parental allele, it is more likely that this is due to the homozygous parent having a 
     null allele. In this case, the locus is again discarded.

    -Also, loci with non-mendelian properties indicative of allele dropout or similar are also filtered

	<> vcf_path				## The path to the vcf file 
	<> sex_info_file			## The path to the file containing info on males, females, parents and offspring (see below for example layout)
	<> catalog_path				## The path to the catalog tags file from stacks
	<> offspring_presence_threshold		## The minimum proportion of offspring that must have data for a locus for it to be considered
	<> mendelian_threshold			## The maximum proportion of heterozygotes or homozygotes in the offspring (default 0.75)
	<> het_sex_thresh			## The threshold number of the heterogametic sex that must be heterozygous at a sex linked marker
	<> hom_sex_thresh			## The threshold number of the homogametic sex that must be homozygous at a sex linked marker


	## Pop map file ##
	
	Sample_1	M
	Sample_1	F
	Sample_1	MO
	Sample_1	mo
	Sample_1	FO
	Sample_1	fo

	Male can be denoted by M or m
	Female can be denoted by F or f
	Male offspring can be denoted by MO or mo
	Female offspring can be denoted by FO or fo

    
    """
    #from __future__ import division
    import sys
    import gzip
    import vcf
    from collections import Counter
    import pprint as pp
    from pydoc import help
    

    
    ## Alter the header in the vcf file to be compatible with PyVCF ==========================
    
    alteredvcfpath = "%s%s" % (vcf_path, ".altered")

    oldvcf = open(vcf_path, 'r').readlines()
    alteredvcf = open(alteredvcfpath, 'w')    
    
    for line in oldvcf:
        if "Allele Depth" not in line:
            alteredvcf.write(line)
        elif "Allele Depth" in line:
            line = '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele Depth">\n'
            alteredvcf.write(line)
    alteredvcf.close()
    
    altered_vcf = open(alteredvcfpath, 'r')
    
    ## ========================================================================================
    
    
    #####
    #####
    #####
        
    
    ## Read in VCF, make counters and lists, open output files ================================
    
    myvcf = vcf.Reader(altered_vcf)
        
    Working_dir = "%s/" % (vcf_path.rpartition("/")[0])
    
    loci_used_for_male_map = 0
    loci_used_for_female_map = 0
    
    Loc_w_null_alleles = 0
    loc_w_excess_hom = 0
    
    offspring = []
    male_offspring = []
    female_offspring = []
    
    #ZW_tags = open("%sPutative_ZW_linked_tags.txt" % (MST_path), 'w')
    #XY_tags = open("%sPutative_XY_linked_tags.txt" % (MST_path), 'w')
    log_file = open("%sSL_markers_from_sexed_fam.log" % (Working_dir), 'w')
    
    ZW_tags = []
    XY_tags = []
    
    ## =========================================================================================================

       
    ## Get info for Parents and offspring ======================================================================
    
    sex_file = open(sex_info_file, 'r').readlines()
        
    for line in sex_file: ## Note, only samples in the sex_info file are used
        
        sample_name = line.split()[0]
        sex = line.split()[1]                    
        
        if sex == "F" or sex == "f":
            Mother = sample_name
            log_file.write("Mother = %s" % (sample_name))
        
        elif sex == "M" or sex == "m":
            Father = sample_name
            log_file.write("\nFather = %s\n" % (sample_name))        
    
        elif sex == "O" or sex == "o" or sex == "MO" or sex == "FO" or sex == "mo" or sex == "fo":
            offspring.append(sample_name)
            if sex == "MO" or  sex == "mo":
                male_offspring.append(sample_name)
            elif sex == "FO" or  sex == "fo":
                female_offspring.append(sample_name)
            
    offspring = sorted(offspring) ## sort the list of offspring - used to retain correct order throughout script
    male_offspring = sorted(male_offspring)
    female_offspring = sorted(female_offspring)
    
    
    
    ## ===========================================================================================================

    
    ## Looping over loci to get info =============================================================================
    ## First check loci for signs of wrong genotype calls, assign them to a list of "Good loci".
    ## And then check these good loci for XY or ZW signal    
    
    Good_loci = []
    
    Male_offs_GTs = {}
    Female_offs_GTs = {}
    N_uncalled_male_offs = 0
    N_uncalled_female_offs = 0
    
    
    N_ZW_testers = 0
    N_ZW_testers_with_enough_het_females = 0
    N_ZW_testers_with_enough_het_females_and_hom_males = 0
            
            
    N_XY_testers = 0
    N_XY_testers_with_enough_het_males = 0
    N_XY_testers_with_enough_het_males_and_hom_females = 0

    for record in myvcf:
        
        ## Make all variable defaults, lists and dicts.
        
        Mat_ref = None
        Mat_alt = None
        Pat_ref = None
        Pat_alt = None
        mat_het = None
        pat_het = None
        Mat_genotype = None
        Pat_genotype = None
        Mat_genotype_bases = None
        Pat_genotype_bases = None
        Non_par_allele = False
        N_offspring = 0
        offspr_w_allele_dropout = 0
        N_offsp_hom_for_minor_par_allele = 0
        
        Sample_headers = []
        sample_names = []
        Parental_alleles = []
        offspring_gt_types = []
        
        offspring_gt_types_dict = {}
        offspring_gt_bases_dict = {}
        allele_dropout_scores = {}
        
        ## lists for constructing the lines of the MST output file
        malemap_offspring_haploid_genotypes = []
        malemap_offspring_haploid_genotypes_compli = []
        femalemap_offspring_haploid_genotypes = []
        femalemap_offspring_haploid_genotypes_compli = []
    
    
        Loc_Id = "%s_%s" % (record.ID, record.POS)
        
        log_file.write("\n\nLocus: %s -------------\n" % (Loc_Id))
        
        for sample in record.samples:
            
            name = sample.sample
            genotype = sample['GT']
            allele_dropout_scores[record] = {}
            
            sample_names.append(name)
            
            
            ## Get Mother and Father information
            
            if name == Mother:
                
                if not sample.called == False:
                    Mat_ref = sample.gt_bases.split("/")[0]
                    Mat_alt = sample.gt_bases.split("/")[1]
                    Mat_genotype = genotype
                    Mat_genotype_bases = sample.gt_bases
                    
                    Parental_alleles.append(Mat_ref)
                    Parental_alleles.append(Mat_alt)
                    
                    if sample.gt_type == 0 or sample.gt_type == 2:
                        mat_het = False
                        log_file.write("\nMother is homozygous (%s)" % (sample.gt_bases))
                    elif sample.gt_type == 1:
                        mat_het = True   
                        log_file.write("\nMother is heterozygous (%s)" % (sample.gt_bases))
                elif sample.called == False:
                    log_file.write("\nLocus not called in mother")
                    
                
                    
            elif name == Father:
                
                if not sample.called == False:
                    Pat_ref = sample.gt_bases.split("/")[0]
                    Pat_alt = sample.gt_bases.split("/")[1]
                    Pat_genotype = genotype
                    Pat_genotype_bases = sample.gt_bases
                    
                    Parental_alleles.append(Pat_ref)
                    Parental_alleles.append(Pat_alt)
                    
                    if sample.gt_type == 0 or sample.gt_type == 2:
                        pat_het = False
                        log_file.write("\nFather is homozygous (%s)" % (sample.gt_bases))
                    elif sample.gt_type == 1:
                        pat_het = True
                        log_file.write("\nFather is heterozygous (%s)" % (sample.gt_bases))
                elif sample.called == False:
                    log_file.write("\nLocus not called in father")
            
        
        ## Get parental allele information and filter for unacceptable loci
        
        Par_allele_counts = Counter(Parental_alleles)
        
        
        if None in Par_allele_counts or Mat_genotype == None or Pat_genotype == None: ## If there is a missing call in the parents, discard locus
            log_file.write("\nMissing data in parents, locus discarded")
                        
        elif len(Par_allele_counts) > 2: ## If there are more than two alleles in the parents, discard locus
            log_file.write("\nparents contain more than 2 alleles, locus not used")
            
        elif pat_het == True and mat_het == True: ## If both parents are heterozygous, discard locus
            log_file.write("\nBoth parents heterozygous, locus not used")
            
        elif pat_het == False and mat_het == False: ## If both parents are homozygous, discard locus
            log_file.write("\nBoth parents homozygous, locus not used")
            

            
        elif len(Par_allele_counts) == 2 and not 2 in Par_allele_counts.values(): ## if both parents are called and both do not have the same genotype
            
            for allele in Par_allele_counts:
                if Par_allele_counts[allele] == 3:
                    Par_major_allele = allele
                elif Par_allele_counts[allele] == 1:
                    Par_minor_allele = allele
                    
        
            ## Get offspring information -----------------------------------------------------------------------
                 
            
            ## First - look for signs of allele dropout - using all offspring for this
            
            for sample in record.samples:            
                name = sample.sample
                genotype = sample['GT']            
                
                if name in offspring:  ## Only samples included in the pop_map file are used
                    N_offspring += 1


                    ### Filter sample genotypes for signs of allele droupout 

                    if genotype == None: ## If sample not called, record the lack of genotype
                        
                        offspring_gt_types.append(sample.gt_type) # record hom or het
                        offspring_gt_bases_dict[name] = (sample.gt_type, sample.gt_bases)
                        
                        
                    elif not genotype == None: ## If sample has been called

                        off_ref = sample.gt_bases.split("/")[0] ## get the samples ref allele
                        off_alt = sample.gt_bases.split("/")[1] ## get the samples alt allele
                        
                        if sample.gt_bases == "%s/%s" % (Par_minor_allele,Par_minor_allele): ## If sample homozygous for the parental minor allele, remove sample info, record the event for tallying later.
                            
                            log_file.write("\n%s is homozygous for minor parental allele (%s), likely allele dropout in offspring (if rare) or parent (if common)" % (name, sample.gt_bases))
                            offspr_w_allele_dropout += 1

                            offspring_gt_types.append(None) # Count as missing data
                            offspring_gt_bases_dict[name] = (None, None)
                            N_offsp_hom_for_minor_par_allele += 1

                        elif off_ref not in Parental_alleles or off_alt not in Parental_alleles: ## If sample has an allele that is not in parents, remove locus
                            
                            offspring_gt_types.append(sample.gt_type) # record hom or het
                            offspring_gt_bases_dict[name] = (sample.gt_type, sample.gt_bases)

                            log_file.write("\n%s contains non-parental allele in genotype (%s), likely allele dropout in a Parent" % (name, sample.gt_bases))
                            Non_par_allele = True


                        if off_ref in Parental_alleles and off_alt in Parental_alleles: ## And if the sample's alleles are both in the Parental genotypes then the locus is good
                            
                            offspring_gt_types.append(sample.gt_type) # record hom or het
                            offspring_gt_bases_dict[name] = (sample.gt_type, sample.gt_bases)

                            ## Record male and female genotypes for later use if the locus is ok.

                            if name in male_offspring:
                                if Loc_Id not in Male_offs_GTs:
                                    Male_offs_GTs[Loc_Id] = {}
                                Male_offs_GTs[Loc_Id][name] = (sample.gt_type, sample.gt_bases)

                                if sample.gt_type == None:
                                    N_uncalled_male_offs += 1

                            elif name in female_offspring:
                                if Loc_Id not in Female_offs_GTs:
                                    Female_offs_GTs[Loc_Id] = {}
                                Female_offs_GTs[Loc_Id][name] = (sample.gt_type, sample.gt_bases)

                                if sample.gt_type == None:
                                    N_uncalled_female_offs += 1

            
            counted = Counter(offspring_gt_types) ## Count the numbers of each genotype at this locus for looking at offspring presence threshold

            

            ## Caluculate the percentage of offspring missing at the locus

            if None in counted:
                perc_offspring_missing = counted[None]/N_offspring
            else:
                perc_offspring_missing = 0

            ## calculate the percentages of homozygotes/heterozygotes for mendelian frequency filters below

            perc_hom_0 = counted[0]/(counted[0]+counted[2]+counted[1])
            perc_hom_2 = counted[2]/(counted[0]+counted[2]+counted[1])
            perc_het = counted[1]/(counted[0]+counted[2]+counted[1])
            perc_off_w_all_dropout = offspr_w_allele_dropout/N_offspring                  

            ## Filter loci with non-mendelian genotype proportions and write data to repsective map files

            if perc_offspring_missing > 1-offspring_presence_threshold:
                log_file.write("\nToo many missing genotypes, locus not used")
                #print "Too many missing genotypes, locus not used"

            elif Non_par_allele == True:
                log_file.write("\nNon parental allele found, locus not used")
                #print "Non parental allele found, locus not used"

            elif N_offsp_hom_for_minor_par_allele > 1: ## if there is more than one offspring homozygous for the minor parental allele, then it is likely to be allele dropout in the parent, locus is discarded
                log_file.write("\nMore than one offspring homozygous for minor allele, locus not used")
                #print "More than one offspring homozygous for minor allele, locus not used"
            elif perc_hom_0 > mendelian_threshold or perc_hom_2 > mendelian_threshold:
                log_file.write("\nHomozygosity excess (>%s)" % (mendelian_threshold))
                log_file.write("\nLocus discarded, Perc_hom_0= %s, Perc_hom_2= %s" % (perc_hom_0, perc_hom_2))
                loc_w_excess_hom += 1
                #print "Homozygosity excess"
            else:
                #print "Locus passed null allele filters"
                Good_loci.append(Loc_Id)


            ## Now look through offspring again, to test for XY, look only at loci where Father is heterozygous and find those that 
            ## are heterozygous in all (or more than a threshold number) of males

            ##============================================================================================================================= 
            ## Now look specifically at loci that are het in only one parent. =============================================================
            ## ============================================================================================================================

            ## Get total number of CALLED females

            N_called_fem_offspring = len(female_offspring) - N_uncalled_female_offs
            N_called_male_offspring = len(male_offspring) - N_uncalled_male_offs
            
            

            

            if Loc_Id in Good_loci:            

                if pat_het == False and mat_het == True: ## If mother is het and father is hom, use to test for ZW
                    log_file.write("\nMother is heterozygous, Father is homozygous, locus used to test for ZW")
                    
                    N_ZW_testers += 1
                    
                    ## Now look through the male and female offspring - look for ZW patterns

                    N_het_females = 0
                    N_hom_males = 0

                    ## Find out the percentage of females called that were heterozygous
                    
                    for female in Female_offs_GTs[Loc_Id]:
                        if Female_offs_GTs[Loc_Id][female][0] == 1 and Mat_alt in Female_offs_GTs[Loc_Id][female][1].split("/"): ## If female is heterozygous count it. 
                            N_het_females += 1
                    
                    ## And then find out the percentage of males called that were homozygous
                    
                    for male in Male_offs_GTs[Loc_Id]:
                        if Male_offs_GTs[Loc_Id][male][0] in [0,2]: ## If male is homozygous count it. 
                            N_hom_males += 1
                        
                    

                    if N_het_females > 1:
                        perc_het_females = N_het_females/N_called_fem_offspring
                        log_file.write("\nN called female offspring = %s" % N_called_fem_offspring)
                        log_file.write("\nN heterozygous females = %s" % N_het_females)
                        log_file.write("\nPercent heterozygous female offspring = %s" % perc_het_females)
                        
                        N_ZW_testers_with_enough_het_females += 1
                        
                        if perc_het_females >= het_sex_thresh:
                            log_file.write("\nEnough heterozygous females!")

                            if N_hom_males > 1:
                                perc_hom_males = N_hom_males/N_called_male_offspring
                                log_file.write("\nPercent homozygous males = %s" % perc_hom_males)
                                
                                if perc_hom_males >= hom_sex_thresh:
                                    ZW_tags.append(Loc_Id)
                                    log_file.write("\nEnough homozygous males!")
                                    log_file.write("\n##Locus supports ZW system!")
                                    N_ZW_testers_with_enough_het_females_and_hom_males += 1
                                else:
                                    log_file.write("\nNot enough homozygous males")
                            else:
                                log_file.write("\nNo homozygous male offspring")
                        else:
                            log_file.write("\nNot enough heterozygous female offspring")
                    else:
                        log_file.write("\nNo heterozygous female offspring")

                elif pat_het == True and mat_het == False: ## If father is het and mother is hom, use to test for XY
                    log_file.write("\nFather is heterozygous, Mother is homozygous, locus used to test for XY")
                    N_XY_testers += 1
                    
                    
                    N_het_males = 0
                    N_hom_females = 0
                    
                    ## Find out the percentage of males called that were heterozygous
                    
                    for male in Male_offs_GTs[Loc_Id]:
                        if Male_offs_GTs[Loc_Id][male][0] == 1 and Pat_alt in Male_offs_GTs[Loc_Id][male][1].split("/"): ## If male is heterozygous count it. 
                            N_het_males += 1
                    
                    ## Find out the percentage of females called that were homozygous
                    
                    for female in Female_offs_GTs[Loc_Id]:
                        if Female_offs_GTs[Loc_Id][female][0] in [0,2]: ## If female is homozygous count it. 
                            N_hom_females += 1
                            
                    

                    if N_het_males > 1:
                        perc_het_males = N_het_males/N_called_male_offspring
                        
                        log_file.write("\nN called male offspring = %s" % N_called_male_offspring)
                        log_file.write("\nN heterozygous male offspring = %s" % N_het_males)
                        log_file.write("\nPercent heterozygous males = %s" % perc_het_males)
                        
                        if perc_het_males >= het_sex_thresh:
                            log_file.write("\nEnough heterozygous males!")
                            N_XY_testers_with_enough_het_males += 1
                            
                            if N_hom_females > 1:
                                perc_hom_females = N_hom_females/N_called_fem_offspring
                                log_file.write("\nPercent homozygous females = %s" % perc_hom_females)
                                
                                if perc_hom_females >= hom_sex_thresh:
                                    XY_tags.append(Loc_Id)
                                    log_file.write("\nEnough homozygous females!")
                                    log_file.write("\n##Locus supports XY system!")
                                    N_XY_testers_with_enough_het_males_and_hom_females += 1
                                else:
                                    log_file.write("\nNot enough homozygous females")

                            else:
                                log_file.write("\nNo homozygous female offspring")        
                        else:
                            log_file.write("\nNot enough heterozygous male offspring")           
                    else:
                        log_file.write("\nNo heterozygous male offspring")


    print "Number of good loci = %s" % len(Good_loci)

    print "\nN loci suitable for XY testing:", N_XY_testers
    print "N XY test loci with enough heterozygous males:", N_XY_testers_with_enough_het_males
    print "N loci that fit the specified XY criteria:", N_XY_testers_with_enough_het_males_and_hom_females
    
    print "\nN loci suitable for ZW testing:", N_ZW_testers
    print "N ZW test loci with enough heterozygous females:", N_ZW_testers_with_enough_het_females
    print "N loci that fit the specified ZW criteria:", N_ZW_testers_with_enough_het_females_and_hom_males
    

    
    ## Finally, make fasta files and vcf files of the sex linked markers
    
    XYset = set(XY_tags)
    ZWset = set(ZW_tags)
    
    if catalog_path.endswith("gz"):
        catalog = gzip.open(catalog_path, 'r').readlines()
    else:
        catalog = open(catalog_path, 'r').readlines()

        
    ## Outdir same as vcf dir:

    outdir = vcf_path.rpartition("/")[0]

    if len(XYset) > 0:
        Putative_XYlinked_makers_file = open("%s/%s" %(outdir, "Putative_XYlinked_makers.fa"), 'w')

        for locus in XYset:
            for tag in catalog:
                if locus.split("_")[0] == tag.split()[2]:
                    Putative_XYlinked_makers_file.write(">X_linkedLocusID_%s\n" % locus)
                    Putative_XYlinked_makers_file.write("%s\n" % (tag.split()[8]))
        Putative_XYlinked_makers_file.close()

    if len(ZWset) > 0:
        Putative_ZWlinked_makers_file = open("%s/%s" %(outdir, "Putative_ZWlinked_makers.fa"), 'w')

        for locus in ZWset:
            for tag in catalog:
                if locus.split("_")[0] == tag.split()[2]:
                    Putative_ZWlinked_makers_file.write(">Z_linked|LocusID_%s\n" % locus)
                    Putative_ZWlinked_makers_file.write("%s\n" % (tag.split()[8]))
        Putative_ZWlinked_makers_file.close()

    return XY_tags, ZW_tags





def Super_SLM_finder(Parameter_dict, execute_seq = "111", write_files = True):
    
    """
    
    Parameter_dict - contains all data and settings for all three approaches
    execute_seq -    (default = "111") tells the wrapper which approach to use, e.g. 111 means run all, 101 means run approach 1 and 3 but not 2.
    
    
    -----------------------------------------------------------------------------------------------------------------------------------------------------
    
    Approach 1. Identifies sex linked SNPs from a VCF file using the criteria of female-male allele freq 

    Workflow:
    1. Filters loci by presence (number of samples), coverage and maf.
    2. Calculates the allele frequencies for males and females separately
    3. Subtracts male from female frequencies and filter loci that show signs of X or Z linkage
    4. Outputs all male and female frequencies and female-male outputs to a single file called 
       "yourinput.vcf.all_frequencies.tsv" (where yourinput = the name and path of your vcf file). Loci 
       identified as X or Z linked are labelled as such in this file.
    5. Outputs all putative X or Z linked markers to separate fasta files if any are identified.
    6. Outputs a histogram of the distribution of female-male frequencies called "yourinput.vcf.fem-male_freqs.pdf"
    7. All suplus information is recorded to a log file, with a summary at the end of this file.
    
    -----------------------------------------------------------------------------------------------------------------------------------------------------
    
    Approach 2. Identifies sex linked SNPs from a VCF file using the heterozygosity difference between sexes 
    
    Workflow:
    1. Filters loci by presence (number of samples), coverage and maf.
    3. Looks for loci that are heterozygous in all or most of the homogametic sex and homozygous in all or most
       of the heterogametic sex
    5. Outputs all putative X or Z linked markers to separate fasta files if any are identified.
    6. All suplus information is recorded to a log file, with a summary at the end of this file.
    
    -----------------------------------------------------------------------------------------------------------------------------------------------------
    
    Approach 3. Identifies sex linked RAD tags from a catalog file using the presence and absence of tags in each sex 
    
    Workflow:
    # Using only samples specified in the Pop_map file:
    1. Finds tags which are present in <sex_presence_threshold> males and no females (Y linked)
    2. Find tags which are present in <sex_presence_threshold> females and no males (W linked)
    3. Outputs all putative X or W linked tags to a fasta file
    
    ====================================================================================================================================================
    
    ## Example parameter dict
    
    Parameter_dict = {}

    ##### Data ########################

    Parameter_dict['Catalog'] =  "/home/djeffrie/Data/RADseq/Rarvalis/Final_populations_outs/batch_1.catalog.tags.tsv.gz" ## Path to the catalog file - used by all approaches.
    Parameter_dict['VCF'] =  "/home/djeffrie/Data/RADseq/Rarvalis/Final_populations_outs/Guillaumes_vcf/batch_1.vcf" ## path to vcf file (note this will be altered to make header compatible with Pyvcf. New vcf will have same name with ".altered" appended to the end). Used by Approach i) and ii)
    Parameter_dict['Pop_map'] = "/home/djeffrie/Data/RADseq/Rarvalis/Final_populations_outs/sex_info_ID.txt" ## path to population map file containing sex information. Same format as Stacks pop map file. Used by all approaches.

    ###### threshold parameters #######

    # 1. Frequency approach
    Parameter_dict['X_or_Z_freq_threshold'] = 0.4  ## (Default = 0.4) The lower threshold for the freq caluclation to find sex linked snps, e.g. for an XY system, a threshold of 0.4 means that f(F) - f(M) can be >= 0.4 and <= 0.6 (the upper threshold is automatically calculated to be the same distance above 0.5 as the lower threshold is below 0.5) 
    Parameter_dict['sample_presence_cutoff1'] = 0.75 ## (Default = 0.75) a locus must be called in at least this proportion of all samples (not within populations) to be considered
    Parameter_dict['coverage_threshold1'] = 3 ## (Default = 3) a locus must have at least this threshold in a sample to be considered for that sample. Note that loci below this threshold will be removed from a sample, and this can push the locus below the sample presence cut-off, which will then remove the locus.
    Parameter_dict['maf_threshold1'] =  0.05 ## (Default = 0.05) minor allele frequency cutoff for a locus across all samples. 

    # 2. Heterozygosity approach
    Parameter_dict['homogamtic_homozygosity_threshold'] = 0.9 ## (Default = 0.9) The minimum number of the homogametic sex which must not have the tag for that tag to be considered linked to the sex-limited chromosome
    Parameter_dict['heterogamtic_heterozygosity_threshold'] = 0.5 ## (Default = 0.5) The lower threshold for the proportion of heterozygotes in the heterogametic sex at a locus 
    Parameter_dict['sample_presence_cutoff2'] = 0.75 ## (Default = 0.75) a locus must be called in at least this proportion of all samples (not within populations) to be considered
    Parameter_dict['coverage_threshold2'] = 3 ## (Default = 3) a locus must have at least this threshold in a sample to be considered for that sample. Note that loci bels this threshold will be removed from a sample, and this can push the locus below the sample presence cut-off, which will then remove the locus.
    Parameter_dict['maf_threshold2'] = 0.05 ## (Default = 0.05) minor allele frequency cutoff for a locus across all samples. 

    # 3. Sex specific presence or absence approach
    Parameter_dict['sex_presence_threshold'] =  0.5 ## (Default = 0.5) The minimum percenatage of the heterogametic sex that a tag must be present in.

    
    """

    from matplotlib import pyplot as plt
    from matplotlib_venn import venn2, venn3
    import gzip
    import time
 
    now = time.strftime("%c")
    
    Xlinked_freq = []
    Zlinked_freq = []
    freq_log = []
    Xlinked_het = []
    Zlinked_het = []
    het_log = []
    Ylinked_tags = []
    Wlinked_tags = []
    pres_abs_log = []
    
    
    ## freq approach ## -------------------------------------------------------------------
    
    if execute_seq[0] == 1 or execute_seq[0] == "1":
        Xlinked_freq, Zlinked_freq, freq_log = SL_snp_freq(Parameter_dict['VCF'], \
                                                 Parameter_dict['Pop_map'],\
                                                 Parameter_dict['Catalog'],\
                                                 Parameter_dict['X_or_Z_freq_threshold'],\
                                                 Parameter_dict['sample_presence_cutoff1'],\
                                                 Parameter_dict['coverage_threshold1'],\
                                                 Parameter_dict['maf_threshold1'],\
						 Parameter_dict['homogametic_REF_allele_freq'])
        
    ## Heterozygosity approach ## ---------------------------------------------------------
    
    if execute_seq[1] == 1 or execute_seq[1] == "1":
        Xlinked_het, Zlinked_het, het_log = SL_snp_het(Parameter_dict['VCF'],\
                   Parameter_dict['Pop_map'],\
                   Parameter_dict['Catalog'],\
                   Parameter_dict['homogamtic_homozygosity_threshold'],\
                   Parameter_dict['heterogamtic_heterozygosity_threshold'],\
                   Parameter_dict['sample_presence_cutoff2'],\
                   Parameter_dict['coverage_threshold2'],\
                   Parameter_dict['maf_threshold2'])
    
    
    ## Presence absence approach ## --------------------------------------------------------
    
    if execute_seq[2] == 1 or execute_seq[2] == "1":
        Ylinked_tags, Wlinked_tags, pres_abs_log = SL_tag_finder(Parameter_dict['Catalog'],\
                      Parameter_dict['Pop_map'],\
                      Parameter_dict['sex_presence_threshold'])

    
    ## Merge the tag lists and draw Venn diagrams for loci identified by each method. 
    
    if execute_seq == "111":
    
        ## For XY tags

        fig = plt.figure(figsize= (15,15))
        fig.add_subplot(1,2,1)
        a = venn3([Xlinked_freq,Xlinked_het, Ylinked_tags], ["Xlinked_freq","Xlinked_het", "Ylinked_tags"])
        plt.title("XY tags identified")

        XYset = Xlinked_freq.union(Xlinked_het,Ylinked_tags)        

        ## For ZW tags

        fig.add_subplot(1,2,2)
        b = venn3([Zlinked_freq, Zlinked_het, Wlinked_tags], ["Zlinked_freq","Zlinked_het", "Wlinked_tags"])
        plt.title("ZW tags identified")

        ZWset = Zlinked_freq.union(Zlinked_het,Wlinked_tags)

    elif execute_seq == "101":

        ## For XY tags

        fig = plt.figure(figsize= (15,15))
        fig.add_subplot(1,2,1)
        venn2([Xlinked_freq, Ylinked_tags], ["Xlinked_freq", "Ylinked_tags"])
        plt.title("XY tags identified")
        XYset = Xlinked_freq.union(Ylinked_tags)

        ## For ZW tags

        fig.add_subplot(1,2,2)
        venn3([Zlinked_freq, Wlinked_tags], ["Zlinked_freq", "Wlinked_tags"])
        ZWset = Zlinked_freq.union(Wlinked_tags)

    elif execute_seq == "011":

        ## For XY tags

        fig = plt.figure(figsize= (15,15))
        fig.add_subplot(1,2,1)
        venn2([Xlinked_het, Ylinked_tags], ["Xlinked_freq", "Ylinked_tags"])
        plt.title("XY tags identified")
        XYset = Xlinked_het.union(Ylinked_tags)

        ## For ZW tags

        fig.add_subplot(1,2,2)
        venn2([Zlinked_het, Wlinked_tags], ["Zlinked_freq", "Wlinked_tags"])
        ZWset = Zlinked_het.union(Wlinked_tags)
        plt.title("ZW tags identified")

    elif execute_seq == "110":

        ## For XY tags

        fig = plt.figure(figsize= (15,15))
        fig.add_subplot(1,2,1)
        venn2([Xlinked_freq,Xlinked_het], ["Xlinked_freq", "Ylinked_tags"])
        plt.title("XY tags identified")
        XYset = Xlinked_freq.union(Xlinked_het)

        ## For ZW tags

        fig.add_subplot(1,2,2)
        venn2([Zlinked_freq,Zlinked_het], ["Zlinked_freq", "Wlinked_tags"])
        plt.title("ZW tags identified")
        ZWset = Zlinked_freq.union(Zlinked_het)


    elif execute_seq == "100":

        XYset = Xlinked_freq
        ZWset = Zlinked_freq

    elif execute_seq == "010":

        XYset = Xlinked_het
        ZWset = Zlinked_het

    elif execute_seq == "001":

        XYset = Ylinked_tags
        ZWset = Wlinked_tags
    
    
    ## Now just need to output the fasta and its done
    if write_files == True:
    
	    if Parameter_dict['Catalog'].endswith("gz"):
	        catalog = gzip.open(Parameter_dict['Catalog'], 'r').readlines()
	    else:
	        catalog = open(Parameter_dict['Catalog'], 'r').readlines()
	        
	    ## Outdir same as catalog dir:
	    
	    outdir = Parameter_dict['Catalog'].rpartition("/")[0]
	        
	    if len(XYset) > 0:
	        Putative_XYlinked_makers_file = open("%s/%s" %(outdir, "Putative_XYlinked_makers.fa"), 'w')
	
	        ## record which method found which tag
	        for locus in XYset:
	            
	            if locus in Xlinked_het and locus in Xlinked_freq:
	                method = "Xlinked_het_Xlinked_freq"
	            
	            elif locus in Xlinked_freq:
	                method = "Xlinked_freq"
	            
	            elif locus in Xlinked_het:
	                method = "Xlinked_het"
	            
	            elif locus in Ylinked_tags:
	                method = "Ylinked_tags"
	
	            for tag in catalog:
	                if locus.split("_")[0] == tag.split()[2]:
	                    Putative_XYlinked_makers_file.write(">X_linkedLocusID_%s_%s\n" % (locus, method))
	                    Putative_XYlinked_makers_file.write("%s\n" % (tag.split()[8]))
	        Putative_XYlinked_makers_file.close()

	    if len(ZWset) > 0:
	        Putative_ZWlinked_makers_file = open("%s/%s" %(outdir, "Putative_ZWlinked_makers.fa"), 'w')
	        
	        ## record which method found which tag
	        for locus in ZWset:
	
	            if locus in Zlinked_het and locus in Zlinked_freq:
	                method = "Zlinked_het_Zlinked_freq" 
	            
	            elif locus in Zlinked_freq: 
	                method = "Zlinked_freq"

        	    elif locus in Zlinked_het:
	                method = "Zlinked_het"    	
           
	            elif locus in Wlinked_tags:
	                method = "Wlinked_tags"
	                
	            for tag in catalog:
	                if locus.split("_")[0] == tag.split()[2]:
	                    Putative_ZWlinked_makers_file.write(">Z_linked|LocusID_%s_%s\n" % (locus, method))
	                    Putative_ZWlinked_makers_file.write("%s\n" % (tag.split()[8]))
	        Putative_ZWlinked_makers_file.close()
    
    
    
	    ### Merge log files ###
    
	    master_log = open("%s/%s" %(outdir, "Super_Sex_linked_marker.log"), 'w')
	    master_log.write("### Analyses run on at " + time.strftime("%c")+"\n\n")
	    
	    master_log.write("\n\n##### Approach 1. SNP frequencies #####\n")
	    
	    if len(freq_log) > 0:
	        for line in freq_log:
	            master_log.write(line)
	    else:
	        master_log.write("\n\n## NOT USED ##\n")
	        
	        
	    master_log.write("\n\n##### Approach 2. SNP heterozygosities #####\n")
	    
	    if len(het_log) > 0:
	        for line in het_log:
	            master_log.write(line)
	    else:
	        master_log.write("\n\n## NOT USED ##\n")
	    
	    master_log.write("\n\n##### Approach 3. Sex specific tags #####\n")
	    
	    if len(pres_abs_log) > 0:    
	        for line in pres_abs_log:
	            master_log.write(line)
	    else:
	        master_log.write("\n\n## NOT USED ##\n")
	    
	    master_log.write("\n\n##### Final summary #####\n")
	    master_log.write("\n ## After merging tags accross methods ## \n")
	    
	    master_log.write("Final number of XY tags = %s\n" % len(XYset))
	    master_log.write("Final number of ZW tags = %s\n" % len(ZWset))
	    master_log.write("Sex linked tags outputted to fastas 'Putative_XYlinked_makers.fa' and Putative_ZWlinked_makers.fa in the directory %s" % outdir)
	    master_log.close()
       	    print "Sex linked tags outputted to fastas 'Putative_XYlinked_makers.fa' and Putative_ZWlinked_makers.fa"
	    print "in the directory %s" % outdir

    ## output summary
    
    print "\n ## After merging tags accross methods ## \n"
    
    print "Final number of XY tags = %s" % len(XYset)
    print "Final number of ZW tags = %s" % len(ZWset)
                                                       
    
    return XYset, ZWset



def BlastParseExtra(infile, genome_fa, best_hit_criteria, Eval_thresh, get_frags_switch = 1, window_size = 2000, verb = 0):
    """
    
    Usage: BlastParseExtra.py  <blast_xml_output> <genome_fasta> <best_hit_criteria> <Eval_thresh> <get_frags_switch> <window_size> <verb>
    
    <blast_xml_output>  -  absolute path to blast xml output
    <genome_fasta>      -  absolute path to genome fasta file
    <best_hit_criteria> -  factor by which the best e-value must be lower than the second best (recommended: 1e-5)
    <Eval_thresh>       -  e-value threshold for unique alignments
    <window_size>       -  size of the window (in bp) around the mapping coordinates used to extract the subject scaffold segment
    <get_frags_switch>  -  switch (1 = on, 0 = off) for getting fragments around the mapping
    <verb>              -  become verbose
    
    This script first filters the mappings in the <blast_xml_output> for uniq hits with evalues better than
    <Eval_thresh> or for multi hits where the best hit is <best_hit_criteria> orders of magnitude better than 
    the second.
    
    It will then retrieve a segment of the scaffold from <genome_fasta> which is + and - the <window_size> around
    the mapping coordinates for each query. If the ends of the scaffold are not within this window, then the 
    length of the segment will be (length of mapped query sequence + 2 x <window_size>). However if an end of a
    scaffold is within this window, the segment will be trimmed to this length.
    
    """


    import sys
    from Bio.Blast import NCBIXML
    from Bio import SeqIO
    import gzip

    if infile.endswith("gz"):
	handle = gzip.open(infile, 'r')
    else:
	handle = open(infile, 'r')

    blast = NCBIXML.parse(handle)

    good_blast_outs = {}
    multi_counter = 0
    unique_counter = 0
    ## From Alan's script: Returns blast hits only when the best e-value is 5 orders of magnitude better than the second best.

    for record in blast :
        if len(record.alignments)==1:
            if record.alignments[0].hsps[0].expect <= Eval_thresh:
                unique_counter += 1
                good_blast_outs[record.query] = {}
                good_blast_outs[record.query]["Ref_hit_id"] = str(record.alignments[0].hit_def)
                good_blast_outs[record.query]["Evalue"] = float(record.alignments[0].hsps[0].expect)
                good_blast_outs[record.query]["Hit_start_coord"] = int(record.alignments[0].hsps[0].sbjct_start)
                good_blast_outs[record.query]["Hit_end_coord"] = int(record.alignments[0].hsps[0].sbjct_end)
                #print "Uniq\t%s\t%s\t%s\t%s\t%s" % (record.query, good_blast_outs[record.query]["Ref_hit_id"], good_blast_outs[record.query]["Evalue"], good_blast_outs[record.query]["Hit_start_coord"], good_blast_outs[record.query]["Hit_end_coord"])


        elif len(record.alignments)>1:
            if all([record.alignments[0].hsps[0].expect <= Eval_thresh, record.alignments[0].hsps[0].expect < best_hit_criteria * record.alignments[1].hsps[0].expect]):
                multi_counter += 1
                good_blast_outs[record.query] = {}
                good_blast_outs[record.query]["Ref_hit_id"] = str(record.alignments[0].hit_def)
                good_blast_outs[record.query]["Evalue"] = float(record.alignments[0].hsps[0].expect)
                good_blast_outs[record.query]["Hit_start_coord"] = int(record.alignments[0].hsps[0].sbjct_start)
                good_blast_outs[record.query]["Hit_end_coord"] = int(record.alignments[0].hsps[0].sbjct_end)

    if verb == 1:

        print "Number of multi-alingments kept:", multi_counter
        print "Number of unique alingments kept:", unique_counter


    Rtemp_summary_out = open("%s/blast_%s_summary.out" % (infile.rpartition("/")[0], window_size), 'w')
    Rtemp_summary_out.write("query\tsubject\n")

    if get_frags_switch == 1:
        
        Rtemp_summary_out = open("%s/blast_%s_summary.out" % (infile.rpartition("/")[0], window_size), 'w')
        Rtemp_summary_out.write("query\tsubject\n")

        print "Getting subject scaffold segments from %s . . . " % (genome_fa)

        Rtemp_chunks = open("%s/blast_%s_chunks.fa" % (infile.rpartition("/")[0], window_size), 'w')

        segment_counter = 0

        Rtemp = SeqIO.parse(genome_fa, "fasta")

        recorded_scaffolds = []

        for scaffold in Rtemp:

            for query in good_blast_outs:

                if scaffold.id == good_blast_outs[query]["Ref_hit_id"]: ## If the scaffold has a hit

                    Rtemp_summary_out.write("%s\t%s\n" % (query, scaffold.id))

                    if scaffold.id not in recorded_scaffolds:

                        recorded_scaffolds.append(scaffold.id)

                        if good_blast_outs[query]["Hit_start_coord"] - window_size <= 0 and good_blast_outs[query]["Hit_end_coord"] + window_size >= len(scaffold.seq): # if the beginning and if the end of the scaffold is within the rang of the window

                            SeqIO.write(scaffold, Rtemp_chunks, 'fasta') ## just print whole scaffold
                            segment_counter += 1


                        elif good_blast_outs[query]["Hit_start_coord"] - window_size <= 0 and good_blast_outs[query]["Hit_end_coord"] + window_size < len(scaffold.seq): ## or if the begninning is in range of the window but the end isn't

                            SeqIO.write(scaffold[:good_blast_outs[query]["Hit_end_coord"]+ window_size], Rtemp_chunks, 'fasta') ## print from beginning to upper end of window
                            segment_counter += 1


                        elif good_blast_outs[query]["Hit_start_coord"] - window_size > 0 and good_blast_outs[query]["Hit_end_coord"] + window_size >= len(scaffold.seq): ## or the end of the scaffold is in range but the beginning isnt

                            SeqIO.write(scaffold[good_blast_outs[query]["Hit_start_coord"]- window_size:], Rtemp_chunks, 'fasta') ## print from lower window limit to the end of the scaffold
                            segment_counter += 1

                        elif good_blast_outs[query]["Hit_end_coord"] + window_size < len(scaffold.seq) and good_blast_outs[query]["Hit_start_coord"] - window_size > 0: ## or if neither end of the scaffold is in range of the window
                            SeqIO.write(scaffold[good_blast_outs[query]["Hit_start_coord"]- window_size:good_blast_outs[query]["Hit_end_coord"] + window_size], Rtemp_chunks, 'fasta') ## print from lower window limit to the end of the scaffold
                            segment_counter += 1
        
        print "%s sequence scaffold segments are in %s/blast_%s_chunks.fa" % (segment_counter, infile.rpartition("/")[0], window_size)
        
        Rtemp_chunks.close()
    
    Rtemp_summary_out.close()
    
    return good_blast_outs

    

def Sexy_LEPmap_plotter(LEPmap_out, sex, Interesting_markers, length_measure = "D", N_LGs = 10, axis_offset = 10, plot_out_name = "LEP_map_plot.pdf"):
    
    """
    
    LEPmap_out        - Path to outfile of ordered markers in LGs from the "OrderMarkers" module of LEPmap
    sex               - Sex for which to plot map 1 = male, 2 = female, use 3 for both sna 4 for sex-averaged maps. 
    length_measure    - Order linkage groups by number of SNPs in each LG "S" or by the distance "D"
    N_LGs             - Number of linkage groups to plot (will plot up to the Nth best LGs in terms of length)
    axis_offset       - Parameter to control the position of the linkage group labels. 
    
    """
    
    
    from matplotlib import pyplot as plt
    import sys
    
    lmap = open(LEPmap_out, 'r').readlines()
    
    
    ## Get the data into easily accessible dictionaries

    Female_LGs = {}
    Male_LGs = {}

    for line in lmap:
        if line.startswith("#***") or line.startswith("***") :
            LG = line.split()[3]
            Male_LGs[LG] = []
            Female_LGs[LG] = []
        elif not line.startswith("#"):
            marker_number = line.split()[0]
            male_pos = float(line.split()[1])
            female_pos = float(line.split()[2])

            Male_LGs[LG].append((marker_number, male_pos))
            Female_LGs[LG].append((marker_number, female_pos)) 
            
    
    ## find the lengths of the LGs and order them 
    
    LG_length_SNP = []
    LG_length_distance = []
    
    if sex == 1:

        for LG in Male_LGs:
            LG_length_SNP.append((LG, len(Male_LGs[LG])))
            LG_length_distance.append((LG, max(Male_LGs[LG], key=lambda x: x[1])[1]))
    
    elif sex == 2:
        
        for LG in Female_LGs:
            LG_length_SNP.append((LG, len(Female_LGs[LG])))
            LG_length_distance.append((LG, max(Female_LGs[LG], key=lambda x: x[1])[1]))
            
    elif sex == 3:
        
        Male_LG_length_SNP = []
        Male_LG_length_distance = []
        
        Female_LG_length_SNP = []
        Female_LG_length_distance = []
        
        for LG in Male_LGs:
            Male_LG_length_SNP.append((LG, len(Male_LGs[LG])))
            Male_LG_length_distance.append((LG, max(Male_LGs[LG], key=lambda x: x[1])[1]))
        
        for LG in Female_LGs:
            Female_LG_length_SNP.append((LG, len(Female_LGs[LG])))
            Female_LG_length_distance.append((LG, max(Female_LGs[LG], key=lambda x: x[1])[1]))
        

        
        Male_LG_length_SNP_sorted = sorted(Male_LG_length_SNP, key=lambda x: x[1], reverse = True)
        Male_LG_length_distance_sorted = sorted(Male_LG_length_distance, key=lambda x: x[1], reverse = True)

        Female_LG_length_SNP_sorted = sorted(Female_LG_length_SNP, key=lambda x: x[1], reverse = True)
        Femle_LG_length_distance_sorted = sorted(Female_LG_length_distance, key=lambda x: x[1], reverse = True)
        
        Both_LG_length_SNP_sorted = Male_LG_length_SNP_sorted + Femle_LG_length_distance_sorted
        Both_LG_length_distance_sorted = Male_LG_length_distance_sorted + Femle_LG_length_distance_sorted

    LG_length_SNP_sorted = sorted(LG_length_SNP, key=lambda x: x[1], reverse = True)
    LG_length_distance_sorted = sorted(LG_length_distance, key=lambda x: x[1], reverse = True)

    
    ## Get the N longest scaffolds

    plot_LGs_SNP_length = []
    plot_LGs_distance = []

    if sex == 1 or sex == 2:

        for i in LG_length_SNP_sorted[:N_LGs]:
            plot_LGs_SNP_length.append(i[0])

        for i in LG_length_distance_sorted[:N_LGs]:
            plot_LGs_distance.append(i[0])
    
    elif sex == 3: ## No point doing "S" here, there will be no difference. 
        if length_measure == "S":
            sys.exit("Male and female will be the same for S, so not plotting")
        
        Male_plot_LGs_distance = []
        Female_plot_LGs_distance = []
        
        for i in Male_LG_length_SNP_sorted[:N_LGs]:
            Male_plot_LGs_distance.append(i[0])
            
            
        for i in Femle_LG_length_distance_sorted[:N_LGs]:
            Female_plot_LGs_distance.append(i[0])
        
        
    # And now plot 

    plt.figure(figsize=(40,30))
    plt.gca().invert_yaxis()
    plt.axis('off')
    
    N = 10
    LGcounter = 1
    Locus_counter = 0
    Male_Locus_counter = 0
    Female_Locus_counter = 0
    
    if sex == 1:


        if length_measure == "D":

            for LG in plot_LGs_distance:
                x = [LGcounter]*len(Male_LGs[LG])
                y = [i[1] for i in Male_LGs[LG]]
                z = [i for i in Male_LGs[LG]]
                plt.plot(x,y, linewidth = 20, color = "grey", zorder=1)
                plt.scatter(x,y,zorder=2, marker= "_", s = 1000)

                for loc in z:
                    if int(loc[0]) in Interesting_markers:
                        plt.scatter(x[0],loc[1],zorder=3, marker= "_", s = 2000, color = "red")

                plt.text(x[0],min(y)-axis_offset, "LG %s" % LG, horizontalalignment='center')
                Locus_counter += len(y)
                LGcounter += 1
            plt.title("Male linkage map, ordered by LG map length, Nloci = %s" % Locus_counter, fontsize = 20)
            
        elif length_measure == "S":

            for LG in plot_LGs_SNP_length:
                x = [LGcounter]*len(Male_LGs[LG])
                y = [i[1] for i in Male_LGs[LG]]
                z = [i for i in Male_LGs[LG]]
                plt.plot(x,y, linewidth = 20, color = "grey", zorder=1)
                plt.scatter(x,y,zorder=2, marker= "_", s = 1500)

                for loc in z:
                    if int(loc[0]) in Interesting_markers:
                        plt.scatter(x[0],loc[1],zorder=3, marker= "_", s = 2000, color = "red")

                plt.text(x[0],min(y)-axis_offset, "LG %s" % LG, horizontalalignment='center')
                LGcounter += 1
                Locus_counter += len(y)
            plt.title("Male linkage map, ordered by N snps in LG, Nloci = %s" % Locus_counter, fontsize = 20)
                
    elif sex == 2:
        if length_measure == "D":

            for LG in plot_LGs_distance:
                x = [LGcounter]*len(Female_LGs[LG])
                y = [i[1] for i in Female_LGs[LG]]
                z = [i for i in Female_LGs[LG]]
                plt.plot(x,y, linewidth = 20, color = "grey", zorder=1)
                plt.scatter(x,y,zorder=2, marker= "_", s = 1000)

                for loc in z:
                    if int(loc[0]) in Interesting_markers:
                        plt.scatter(x[0],loc[1],zorder=3, marker= "_", s = 2000, color = "red")

                plt.text(x[0],min(y)-axis_offset, "LG %s" % LG, horizontalalignment='center')
                
                LGcounter += 1
                Locus_counter += len(y)
            plt.title("Female linkage map, ordered by LG map length, Nloci = %s" % Locus_counter, fontsize = 20)

        elif length_measure == "S":

            for LG in plot_LGs_SNP_length:
                x = [LGcounter]*len(Female_LGs[LG])
                y = [i[1] for i in Female_LGs[LG]]
                z = [i for i in Female_LGs[LG]]
                plt.plot(x,y, linewidth = 20, color = "grey", zorder=1)
                plt.scatter(x,y,zorder=2, marker= "_", s = 1500)

                for loc in z:
                    if int(loc[0]) in Interesting_markers:
                        plt.scatter(x[0],loc[1],zorder=3, marker= "_", s = 2000, color = "red")

                plt.text(x[0],min(y)-axis_offset, "LG %s" % LG, horizontalalignment='center')
                                
                LGcounter += 1
                Locus_counter += len(y)
            plt.title("Female linkage map, ordered by N snps in LG, Nloci = %s" % Locus_counter, fontsize = 20)

    elif sex == 3:
        if length_measure == "D":

            for LG in Male_plot_LGs_distance:

                xM = [LGcounter]*(len(Male_LGs[LG]))
                yM = [i[1] for i in Male_LGs[LG]]
                zM = [i for i in Male_LGs[LG]]
                

                plt.plot(xM,yM, linewidth = 20, color = "grey", zorder=1)
                plt.scatter(xM,yM,zorder=2, marker= "_", s = 1000)
                
                for loc in zM:
                    if int(loc[0]) in Interesting_markers:
                        plt.scatter(xM[0],loc[1],zorder=3, marker= "_", s = 2000, color = "red")

                plt.text(xM[0],min(yM)-axis_offset, "LG %s" % LG, horizontalalignment='center')
            
                LGcounter += 1
                Male_Locus_counter += len(yM)
                
            LGcounter += 1 ## add a space
            
            for LG in Female_plot_LGs_distance:

                xF = [LGcounter]*(len(Female_LGs[LG]))
                yF = [i[1] for i in Female_LGs[LG]]
                zF = [i for i in Female_LGs[LG]]

                plt.plot(xF,yF, linewidth = 20, color = "grey", zorder=1)
                plt.scatter(xF,yF,zorder=2, marker= "_", s = 1000)

                for loc in zF:
                    if int(loc[0]) in Interesting_markers:
                        plt.scatter(xF[0],loc[1],zorder=3, marker= "_", s = 2000, color = "red")

                plt.text(xF[0],min(yF)-axis_offset, "LG %s" % LG, horizontalalignment='center')
                
                LGcounter += 1
                
                Female_Locus_counter += len(yF)
                
            plt.title("Male and Female linkage maps, ordered by LG map length, Nloci_Male = %s, Nloci_Female = %s" % (Male_Locus_counter, Female_Locus_counter), fontsize = 20)
                
    plt.savefig("%s/%s" % (LEPmap_out.rpartition("/")[0], plot_out_name))




def Null_mapper(sex_linked_markers_path, catalog_path, Genome_fasta, Genome_db, N_its = 10, Species = "", Threads = 1):
    
        
    """
    This function takes a set of sex linked (SL) markers and maps them against the ordered Xenopus genome. It also retrieves <N_its> * random subsets of markers
    from the Stacks RADseq catalog of the same size as the set of SL markers and maps them to the genome. This essentially gives a null distribution for
    the number of hits to expect from randomly selected markers in the genome. Thus, if the number of SL markers which map to a given chromosome is
    significantly higher than this null distrubution, this is good support for these markers indeed being sex linked.    
    
    Note it is not a general function, it is a function which require 1) 10 chromosomes in the reference genome ordering, 2) the 
    right format of scaffold name. 
    
    sex_linked_markers_path - Path to a fasta of sequences of interest
    catalog_path            - Path to the stacks catalog file
    Genome_fasta            - Path to the fasta for the genome
    Genome_db               - Path to the blastn index of the genome
    N_its                   - Number of random mappings to perform
    Species                 - Name of the species
    Threads                 - Number of threads to use in the mapping step
    
    
    """

    
    import gzip
    import random
    import os
    import numpy as np
    from collections import Counter 
    from Bio.Blast.Applications import NcbiblastnCommandline
    import MISC_RAD_tools as MISC
    from matplotlib import pyplot as plt
    
    ## First, get the IDs of the sex linked tags.

    sex_linked_markers = open(sex_linked_markers_path, 'r').readlines()

    SL_tag_IDs = []

    for line in sex_linked_markers:
        if ">" in line:
            SL_tag_ID = line.split("_")[2]
            if SL_tag_ID not in SL_tag_IDs:
                SL_tag_IDs.append(SL_tag_ID)

    print "Got %s sex-linked tag IDs" % (len(set(SL_tag_IDs)))
    
    
    ## Then get the IDs of the tags in the catalog

    catalog = gzip.open(catalog_path, 'r').readlines()

    cat_tag_IDs = []

    for line in catalog[1:]:
        cat_tag_ID = line.split()[2]
        cat_tag_IDs.append(cat_tag_ID)

    print "\nGot %s catalog tag IDs" % (len(cat_tag_IDs))
    
    
    ## Now remove the sex linked markers from the catalog tags
    counter = 0
    for i in SL_tag_IDs:
        if i in cat_tag_IDs:
            cat_tag_IDs.remove(i)
            counter += 1

    print "\nRemoved %s tags from catalog tags" % counter
    
    ## Now get "N" random samples of tags of the same size as the number of sex linked markers.

    random_samples = {}

    N = N_its
    
    print "\nDoing %s Random mapping iterations" % N

    for i in range(N_its):
        random_samples[i] = random.sample(cat_tag_IDs, len(SL_tag_IDs))    

        
    ## Ok, so now I need to get these markers from the catalog and make fasta's for the mapping. 

    wd = "%s/Random_tags/" % sex_linked_markers_path.rpartition("/")[0]

    if not os.path.exists(wd):
        os.mkdir(wd)

    for sample in random_samples:
        out_fasta = open("%s/random_sample_%s.fasta" % (wd,sample), 'w')

        for tag in catalog:
            if tag.split()[2] in random_samples[sample]:
                seq = tag.split()[8]
                out_fasta.write(">%s\n%s\n" % (tag.split()[2], seq))
        out_fasta.close()
        
    ## Now I need to do the mapping and catch the outputs

    ## First get fasta paths

    fasta_paths = []
    N_fastas = 0
    for root, dirs, files in os.walk(wd):
        for fasta in files:
            if fasta.endswith("fasta"):
                if N_fastas <= N:
                    fasta_paths.append("%s/%s" % (root, fasta))
                    N_fastas += 1
    ## Set mapping filtering criteria (can make these options later if needed)

    best_hit_crit = 1e-5
    Eval_threshold = 1e-20
    Window = 2000

    ## Set up dictionary to catch blast results

    blastn_outs_dict = {}  ## fasta paths can be the dictionary keys

    counter = 1
    for fasta in sorted(fasta_paths):

        print "\nBlasting %s (%s of %s)\n" % (fasta, counter, N)

        blast_out_path = "%s_blastn.xml" % fasta.rpartition(".")[0]

        blastn_cline = NcbiblastnCommandline(query=fasta, db=Genome_db, outfmt=5, out=blast_out_path, num_threads = Threads)
        stdout, stderr = blastn_cline()

        blastn_outs_dict[fasta] = MISC.BlastParseExtra(blast_out_path, Genome_fasta , best_hit_crit, Eval_threshold,Window)
        
        counter += 1

    ## And now map the real sex linked tags

    print "\nBlasting real sex linked markers\n"
    
    real_SLM_blast_out_path = "%s_REAL_SL_TAG_blastn_outs.xml" % sex_linked_markers_path.rpartition(".")[0]

    blastn_cline = NcbiblastnCommandline(query=sex_linked_markers_path, db=Genome_db, outfmt=5, out=real_SLM_blast_out_path, num_threads = Threads)
    stdout, stderr = blastn_cline()

    Real_SL_tag_blastouts = MISC.BlastParseExtra(real_SLM_blast_out_path, Genome_fasta , best_hit_crit, Eval_threshold,Window)

    print "\nReal sex linked marker blastn outputs are here: %s" % real_SLM_blast_out_path

    
    ## Now count the number of mappings per chromosome across the randomised marker samples
    
    All_mapping_counts = {}

    All_mapping_counts["Chr01"] = []
    All_mapping_counts["Chr02"] = []
    All_mapping_counts["Chr03"] = []
    All_mapping_counts["Chr04"] = []
    All_mapping_counts["Chr05"] = []
    All_mapping_counts["Chr06"] = []
    All_mapping_counts["Chr07"] = []
    All_mapping_counts["Chr08"] = []
    All_mapping_counts["Chr09"] = []
    All_mapping_counts["Chr10"] = []

    for i in blastn_outs_dict: ## for each mapping run

        ### Get the mapping counts per chromosome

        Xen_chroms = {}

        Xen_chroms["Chr01"] = 0
        Xen_chroms["Chr02"] = 0
        Xen_chroms["Chr03"] = 0
        Xen_chroms["Chr04"] = 0
        Xen_chroms["Chr05"] = 0
        Xen_chroms["Chr06"] = 0
        Xen_chroms["Chr07"] = 0
        Xen_chroms["Chr08"] = 0
        Xen_chroms["Chr09"] = 0
        Xen_chroms["Chr10"] = 0


        
        hits = []
        hit_Xen_scaff = []

        for record in blastn_outs_dict[i]:
            hits.append(blastn_outs_dict[i][record]["Ref_hit_id"])
            hit_Xen_scaff.append(blastn_outs_dict[i][record]["Ref_hit_id"].split("_")[1])

        counts = Counter(hit_Xen_scaff)

        for chrom in Xen_chroms:
            if chrom in counts:
                Xen_chroms[chrom] = counts[chrom]

        for chrom in Xen_chroms:
            All_mapping_counts[chrom].append(Xen_chroms[chrom])

    ## And do the same for the real sex linked markers

    hits = []
    hit_Xen_scaff = []

    for i in Real_SL_tag_blastouts:
        hits.append(Real_SL_tag_blastouts[i]["Ref_hit_id"])
        hit_Xen_scaff.append(Real_SL_tag_blastouts[i]["Ref_hit_id"].split("_")[1])

    print "\nNumber of hits =", len(hits) 

    counts = Counter(hit_Xen_scaff)

    Real_SL_tagXen_chroms = {}

    Real_SL_tagXen_chroms["Chr01"] = 0
    Real_SL_tagXen_chroms["Chr02"] = 0
    Real_SL_tagXen_chroms["Chr03"] = 0
    Real_SL_tagXen_chroms["Chr04"] = 0
    Real_SL_tagXen_chroms["Chr05"] = 0
    Real_SL_tagXen_chroms["Chr06"] = 0
    Real_SL_tagXen_chroms["Chr07"] = 0
    Real_SL_tagXen_chroms["Chr08"] = 0
    Real_SL_tagXen_chroms["Chr09"] = 0
    Real_SL_tagXen_chroms["Chr10"] = 0


    for i in Real_SL_tagXen_chroms:
        if i in counts:
            Real_SL_tagXen_chroms[i] = counts[i]
    
    
    # And now then I just need to plot this as box plots

    fig = plt.figure(figsize=(20,10), frameon = False, edgecolor = 'none')
    ax = plt.subplot(111)
    pos = 1
    col = "white"
    
    
    max_y = max([max(i) for i in All_mapping_counts.values()] + Real_SL_tagXen_chroms.values()) + 3
    
    for i in sorted(All_mapping_counts.keys()):

        if col == "grey":
            col = "white"
        elif col == "white":
            col = "grey"
        
        ax.bar(left = pos-0.5, width = 1, height= max_y, color= col, edgecolor = col, bottom = 0, alpha = 0.2, zorder = 1)
        vio = ax.violinplot(All_mapping_counts[i], positions=[pos], showmeans = True, showextrema = False )
        
        for pc in vio['bodies']:
            pc.set_facecolor('khaki')
            pc.set_edgecolor('black')
            pc.set_zorder(2)    
            pc.set_alpha(1)
        
        ax.scatter(pos, Real_SL_tagXen_chroms[i], s = 60, zorder = 3, c = "royalblue")
        
        if Real_SL_tagXen_chroms[i] > np.percentile(All_mapping_counts[i], 99):
            ax.plot(pos-0.05, max_y - 1, "*r", zorder = 3)
            ax.plot(pos+0.05, max_y - 1, "*r", zorder = 3)
        elif Real_SL_tagXen_chroms[i] > np.percentile(All_mapping_counts[i], 95):
            ax.plot(pos, max_y - 1, "*r", zorder = 3)
        
        pos += 1

    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Only show ticks on the left spines
    ax.yaxis.set_ticks_position('left')

    ax.set_xlim(0, 11)
    ax.set_ylim(-1, max_y)
    ax.set_xticks(range(1,11))
    ax.set_xticklabels(sorted(All_mapping_counts.keys()))
    ax.set_xlabel("Xenopus chromosome", labelpad=20)
    ax.set_ylabel("Number of mapped markers",labelpad=20)
    ax.get_xaxis().set_tick_params(which='both', direction='out', pad = 15, bottom = 'off', top = 'off')
    plt.title(Species)

    plt.savefig("%s/Per_chromosome_mappings.svg" %  sex_linked_markers_path.rpartition("/")[0])
    
    plt.show()
    
    print "\nAll done, figure saved as svg here: %s/Per_chromosome_mappings_test.svg" %  sex_linked_markers_path.rpartition("/")[0]
    
    return Real_SL_tagXen_chroms, All_mapping_counts


### Function for plotting the missing data in a VCF

### Function for plotting the missing data in a VCF

def Summary_plotter(input_vcf_path, switches = 111, sort = True):
    
    """
    Summary plotter plots some standard metrics for a dataset. Plots missing data per sample, coverage per sample and heterozygosity Vs. coverage accross samples and loci. 
    
    input_vcf_path = full path to vcf
    switches = vector of four 0 or 1 switches for the missing data, coverage and heterozygosity plots in that order. Example: "1101".
    sort = True or False - whether to sort the samples by the plotted metric (True), or keep order in vcf (False)
    """
    
    from matplotlib import pyplot as plt
    from matplotlib import gridspec
    import numpy as np
    import operator
    
    
    
    input_vcf = open(input_vcf_path, 'r').readlines()
    
    samples_ordered = [] ## get samples in the order they are in the vcf
    
    per_sample_missing_data = {} ## record missing data per sample
    per_sample_coverage = {} ## record coverage per sample
    per_sample_heterozygosity = {} ## recordN heterozygous loci per sample
    
    per_locus_coverage = {} ## record the average coverage per locus across samples
    per_locus_heterozygosity = {}  ## record the proportion of heterozygous samples for each locus
    
    
    format_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    
    N_loci = 0
    
    for line in input_vcf:
        
        if line.startswith("#CHROM"):
            for field in line.split():
                if field not in format_fields:
                    if field not in per_sample_missing_data:
                        per_sample_missing_data[field] = 0
                        per_sample_coverage[field] = []
                        per_sample_heterozygosity[field] = 0
                        
                        samples_ordered.append(field)
                        

        elif not line.startswith("#"):
            
            N_loci += 1
            
            record_ID = "%s_%s" % (line.split()[2], line.split()[1])
            
            NS = int(line.split()[7].split(";")[0].split("=")[1])
            
            record = line.split()[9:]
            
            N_samples_het = 0
            
            record_coverage_list = []
            
            sample_index = 0
            
            for sample_field in record:
                
                sample_name = samples_ordered[sample_index]
                                
                DP = sample_field.split(":")[1]
                
                per_sample_coverage[sample_name].append(int(DP))
                record_coverage_list.append(int(DP))
                
                
                GT = sample_field.split(":")[0]
                
                if GT == "1/0" or GT == "0/1":
                    
                    per_sample_heterozygosity[sample_name] += 1
                    
                    N_samples_het += 1                    
                
                elif GT == "./.":
                    per_sample_missing_data[sample_name] += 1
                    
                sample_index += 1
        
            per_locus_heterozygosity[record_ID] = N_samples_het/NS
            per_locus_coverage[record_ID] = np.mean(record_coverage_list)

    
    ### Get the means of the samples
    
    per_sample_coverage_means = {}
        
    for sample in per_sample_coverage:
        per_sample_coverage_means[sample] = np.mean(per_sample_coverage[sample])
    
    

    ####### Make plots ###############
    
    fig = plt.figure(figsize = (30,70))  

    
    if switches[0] == "1": ## make the missing data plot
        
        
        if sort == True:

            ## order the missing data

            missing_sorted_sample_names = [i[0] for i in sorted(per_sample_missing_data.items(), key=operator.itemgetter(1), reverse=True)]
    
        else:
            missing_sorted_sample_names = samples_ordered
        
         
        ## Plot figure
        
        #fig = plt.figure(figsize = (30,10))
        
        ax1 = plt.subplot2grid((4,2),(0,0),colspan=2) 
        
        #ax = fig.add_subplot(111)
        
        bar_pos = 1
        
        x_tick_pos = []
        x_tick_labs = []
        
        for sample in missing_sorted_sample_names:
            x_tick_labs.append(sample)
            x_tick_pos.append(bar_pos)
            ax1.bar(bar_pos, per_sample_missing_data[sample], color = "royalblue", align = "center")
            bar_pos += 1
        
        ax1.hlines(np.mean(per_sample_missing_data.values()), 0, len(per_sample_missing_data))
        
        plt.xticks(x_tick_pos, x_tick_labs, rotation = 90)
        plt.title("Missing data per sample")
        #plt.show()
        
    
    
    if switches[1] == "1": ## make the heterozygosity plot
                

        ## Plot figure
        
        #fig = plt.figure(figsize = (30,10))
        
        #ax = fig.add_subplot(111)
        
        ax2 = plt.subplot2grid((4,2),(1,0),colspan=2)
        
        bar_pos = 1
        
        x_tick_pos = []
        x_tick_labs = []
        
        prop_hets = {}
        
        ## First pass to calculate the heterozyousity as a proportion across all loci for each sample
        
        for sample in per_sample_heterozygosity:
            
            if per_sample_missing_data[sample] == 0:
                N_called = N_loci
            else:
                N_called = N_loci - per_sample_missing_data[sample]
            
            if N_called == 0:
                N_called = 1
            
            het_prop = per_sample_heterozygosity[sample]/N_called ## calculate proportion heterozygous
            
            prop_hets[sample] = het_prop
            
            
        if sort == True:

            ## order the missing data

            heterozygosity_sorted_sample_names = [i[0] for i in sorted(prop_hets.items(), key=operator.itemgetter(1), reverse=True)]  
        else:
            heterozygosity_sorted_sample_names = samples_ordered
            

        
        ## Second pass to plot
        
        for sample in heterozygosity_sorted_sample_names:    
            
            x_tick_labs.append(sample)
            x_tick_pos.append(bar_pos)
            
            ax2.bar(bar_pos, prop_hets[sample] , color = "orange", alpha = 0.7, align = "center")

            bar_pos += 1
        
        
        
        ax2.hlines(np.mean(prop_hets.values()), 0, len(per_sample_missing_data))
        
        plt.xticks(x_tick_pos, x_tick_labs, rotation = 90)
        plt.title("N heterozygous loci per sample")
        #plt.show()
    
    
    if switches[2] == "1": ## make the coverage plot
        
        
        overall_mean = np.mean(per_sample_coverage_means.values())
        
        
        if sort == True:

            ## order the missing data

            coverage_sorted_sample_names = [i[0] for i in sorted(per_sample_coverage_means.items(), key=operator.itemgetter(1), reverse=True)]  
        else:
            coverage_sorted_sample_names = samples_ordered
             
        
        
        #fig = plt.figure(figsize = (30,10))
        
        #ax = fig.add_subplot(111)
        
        ax3 = plt.subplot2grid((4,2),(2,0),colspan=2)
        
        
        vio_pos = 1
        
        x_tick_pos = []
        x_tick_labs = []
        
        for sample in coverage_sorted_sample_names:
            
            x_tick_pos.append(vio_pos)
            x_tick_labs.append(sample)
            
            if sum(per_sample_coverage[sample]) == 0:
                per_sample_coverage[sample].append(1)
                        
            
            vio = ax3.violinplot(per_sample_coverage[sample], positions = [vio_pos], showmeans=True,showextrema = False )
                        
            for pc in vio['bodies']:
                pc.set_facecolor('darkgreen')
                pc.set_edgecolor('black')
                pc.set_alpha(0.5)
                
            vio["cmeans"].set_color("black")
            
            vio_pos += 1
            
        plt.xticks(x_tick_pos, x_tick_labs, rotation = 90)
        plt.title("Locus coverage per sample")
    
            
    if switches[3] == "1": ## make the coverage vs Heterozygosity plots (per sample and per locus). Might take a little time for big datasets
        
        ## no sorting required here
        
        #fig = plt.figure(figsize = (30,15))
        
        #ax = fig.add_subplot(121)
        
        ax4 = plt.subplot2grid((4,2),(3,0))

        
        for sample in per_sample_coverage:
            
            if per_sample_missing_data[sample] == 0:
                N_called = N_loci
            else:
                N_called = N_loci - per_sample_missing_data[sample]
            
            if N_called == 0:
                N_called = 1
            
            ax4.scatter(per_sample_heterozygosity[sample]/N_called, per_sample_coverage_means[sample], color = "black")
            ax4.text((per_sample_heterozygosity[sample]/N_called)+0.002, per_sample_coverage_means[sample], sample[:10])
            
        plt.ylabel("Coverage (reads)")
        plt.xlabel("Prop. heterozygous loci")
        plt.title("Coverage Vs Heterozygosity")

        
        
        #ax1 = fig.add_subplot(122)
        
        ax5 = plt.subplot2grid((4,2),(3,1))
        
        for record in per_locus_coverage:
            
            ax5.scatter(per_locus_heterozygosity[record], per_locus_coverage[record], color = "black")
            
        plt.title("Per locus coverage vs heterozygosity")
        plt.xlabel("Prop. heterozygous samples")
        plt.ylabel("Mean coverage across samples")
        
    
    outfile = "%s_VCF_summary_stats.pdf" % input_vcf_path.rpartition(".")[0]
    
    plt.savefig(outfile)
    
    plt.show()
        
    print "#### DONE ####\n"
    
    print "Number of Loci = %s" %N_loci
    
    print "Average coverage = %s reads" % overall_mean
    
    
        
        
