from __future__ import division

def SL_snp_freq(myvcfpath, popmappath, catalog_tags_file, X_or_Z_freq_threshold = 0.4, sample_presence_cutoff = 0.75, coverage_threshold = 3, maf_threshold = 0.05, homgam_REF_freq = 0.95, verbose = True, plot = True):    
   
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
    
    if verbose == True:
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
    if verbose == True:
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

            if all([percent_samples_present >= sample_presence_cutoff, femREF_count > 0, fem_genotypes > 0, malREF_count > 0, male_genotypes > 0]):
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

    if verbose == True:
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

    if plot == True:
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
    
    if verbose == True:
        print "\n***DONE!***\n"


    return set(Putative_Xlinked_makers), set(Putative_Zlinked_makers), freq_ratios_log



def SL_snp_het(myvcfpath, popmappath, catalog_tags_file, homogamtic_homozygosity_threshold = 0.9, heterogametic_heterozygosity_threshold = 0.5, sample_presence_cutoff = 0.75, coverage_threshold = 3, maf_threshold = 0.05, verbose = True):

    
   
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
    
    
    if verbose == True:
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
    if verbose == True:
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

    if verbose == True:
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

    if verbose == True:
        print "\n ### DONE! ### \n"
   
 
    return set(Putative_Xlinked_makers), set(Putative_Zlinked_makers), het_approach_log




def SL_tag_finder(catalog_tags_file, popmappath, sex_presence_thresh = 0.5, verbose = True):
    
    import gzip
    import math

    ## Note popmappath here should point to a file that contains usual two popmap columns plus a
    ## a column for the ID of each sample . . . 
    
    if verbose == True:
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

    log.append("\nSUMMARY:\nNumber of males: %s\n" % (female_count))
    log.append("Number of males: %s\n" % (male_count))
    log.append("Number of Putative Y linked tags: %s\n" % (len(putative_Ylinked_tags)))
    log.append("Number of Putative W linked tags: %s\n" % (len(putative_Wlinked_tags)))
    
    if verbose == True:    
	## Print summary to STDOUT

        print "\nSUMMARY:\nNumber of males: %s" % (female_count)
        print "Number of males: %s" % (male_count)
        print "Number of Putative Y linked tags: %s" % (len(putative_Ylinked_tags))
        print "Number of Putative W linked tags: %s" % (len(putative_Wlinked_tags))
        print "\n ### DONE! ###\n"
    
    
    ## return list of sex linked tags
    
    return set(putative_Ylinked_tags), set(putative_Wlinked_tags), log




def Super_SLM_finder(Parameter_dict, execute_seq = "111", write_files = True, verbose = True, plot = True):
    
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
## Define a program to help parallelise the analyses
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
						 Parameter_dict['homogametic_REF_allele_freq'],\
						 verbose = verbose,\
						 plot = plot)
        
    ## Heterozygosity approach ## ---------------------------------------------------------
    
    if execute_seq[1] == 1 or execute_seq[1] == "1":
        Xlinked_het, Zlinked_het, het_log = SL_snp_het(Parameter_dict['VCF'],\
                   Parameter_dict['Pop_map'],\
                   Parameter_dict['Catalog'],\
                   Parameter_dict['homogamtic_homozygosity_threshold'],\
                   Parameter_dict['heterogamtic_heterozygosity_threshold'],\
                   Parameter_dict['sample_presence_cutoff2'],\
                   Parameter_dict['coverage_threshold2'],\
                   Parameter_dict['maf_threshold2'],\
                   verbose = verbose)
    
    
    ## Presence absence approach ## --------------------------------------------------------
    
    if execute_seq[2] == 1 or execute_seq[2] == "1":
        Ylinked_tags, Wlinked_tags, pres_abs_log = SL_tag_finder(Parameter_dict['Catalog'],\
                      Parameter_dict['Pop_map'],\
                      Parameter_dict['sex_presence_threshold'],\
                      verbose = verbose)

    
    ## Merge the tag lists and draw Venn diagrams for loci identified by each method. 
    
    if execute_seq == "111":

        XYset = Xlinked_freq.union(Xlinked_het,Ylinked_tags)
        ZWset = Zlinked_freq.union(Zlinked_het,Wlinked_tags)

    if plot == True:
            ## For XY tags
            fig = plt.figure(figsize= (15,15))
            fig.add_subplot(1,2,1)
            a = venn3([Xlinked_freq,Xlinked_het, Ylinked_tags], ["Xlinked_freq","Xlinked_het", "Ylinked_tags"])
            plt.title("XY tags identified")

            ## For ZW tags
            fig.add_subplot(1,2,2)
            b = venn3([Zlinked_freq, Zlinked_het, Wlinked_tags], ["Zlinked_freq","Zlinked_het", "Wlinked_tags"])
            plt.title("ZW tags identified")


    elif execute_seq == "101":
        
        XYset = Xlinked_freq.union(Ylinked_tags)
        ZWset = Zlinked_freq.union(Wlinked_tags)

        if plot == True:
            ## For XY tags
            fig = plt.figure(figsize= (15,15))
            fig.add_subplot(1,2,1)
            venn2([Xlinked_freq, Ylinked_tags], ["Xlinked_freq", "Ylinked_tags"])
            plt.title("XY tags identified")

            ## For ZW tags
            fig.add_subplot(1,2,2)
            venn3([Zlinked_freq, Wlinked_tags], ["Zlinked_freq", "Wlinked_tags"])
            plt.title("ZW tags identified")

    elif execute_seq == "011":
        XYset = Xlinked_het.union(Ylinked_tags)
        ZWset = Zlinked_het.union(Wlinked_tags)

        if plot == True:
            ## For XY tags
            fig = plt.figure(figsize= (15,15))
            fig.add_subplot(1,2,1)
            venn2([Xlinked_het, Ylinked_tags], ["Xlinked_freq", "Ylinked_tags"])
            plt.title("XY tags identified")

            ## For ZW tags
            fig.add_subplot(1,2,2)
            venn2([Zlinked_het, Wlinked_tags], ["Zlinked_freq", "Wlinked_tags"])
            plt.title("ZW tags identified")

    elif execute_seq == "110":
        XYset = Xlinked_freq.union(Xlinked_het)
        ZWset = Zlinked_freq.union(Zlinked_het)
        
        if plot == True:
            ## For XY tags
            fig = plt.figure(figsize= (15,15))
            fig.add_subplot(1,2,1)
            venn2([Xlinked_freq,Xlinked_het], ["Xlinked_freq", "Ylinked_tags"])
            plt.title("XY tags identified")

            ## For ZW tags
            fig.add_subplot(1,2,2)
            venn2([Zlinked_freq,Zlinked_het], ["Zlinked_freq", "Wlinked_tags"])
            plt.title("ZW tags identified")

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
	    
	    outdir = Parameter_dict['VCF'].rpartition("/")[0]
	        
	    if len(XYset) > 0:
	        Putative_XYlinked_makers_file = open("%s/%s_%s_%s_%s_%s_%s_%s_%s_%s.fasta" %(outdir, "XY",Parameter_dict['Pop_map'].rpartition("/")[2].rpartition(".")[0] ,Parameter_dict['sample_presence_cutoff1'],Parameter_dict['coverage_threshold1'],Parameter_dict['homogametic_REF_allele_freq'],Parameter_dict['X_or_Z_freq_threshold'],Parameter_dict['homogamtic_homozygosity_threshold'],Parameter_dict['heterogamtic_heterozygosity_threshold'],Parameter_dict['sex_presence_threshold']), 'w')
	
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
	        Putative_ZWlinked_makers_file = open("%s/%s_%s_%s_%s_%s_%s_%s_%s_%s.fasta" %(outdir, "ZW",Parameter_dict['Pop_map'].rpartition("/")[2].rpartition(".")[0] ,Parameter_dict['sample_presence_cutoff1'],Parameter_dict['coverage_threshold1'],Parameter_dict['homogametic_REF_allele_freq'],Parameter_dict['X_or_Z_freq_threshold'],Parameter_dict['homogamtic_homozygosity_threshold'],Parameter_dict['heterogamtic_heterozygosity_threshold'],Parameter_dict['sex_presence_threshold']), 'w')
	        
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
	                if locus == tag.split()[2]:
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

    if verbose == True:
        print "Sex linked tags outputted to fastas 'Putative_XYlinked_makers.fa' and Putative_ZWlinked_makers.fa"
        print "in the directory %s" % outdir

    ## output summary
    if verbose == True:
        print "\n ## After merging tags accross methods ## \n"
        print "Final number of XY tags = %s" % len(XYset)
        print "Final number of ZW tags = %s" % len(ZWset)
                                                       
    ## Make a more detailed dictionary for outputting

    detailed_out_dict = {}
    detailed_out_dict["XY"] = {}
    detailed_out_dict["XY"]["freq"] = Xlinked_freq
    detailed_out_dict["XY"]["het"] = Xlinked_het
    detailed_out_dict["XY"]["Ytags"] = Ylinked_tags


    detailed_out_dict["ZW"] = {}
    detailed_out_dict["ZW"]["freq"] = Zlinked_freq
    detailed_out_dict["ZW"]["het"] = Zlinked_het
    detailed_out_dict["ZW"]["Wtags"] = Wlinked_tags


    
    return XYset, ZWset, detailed_out_dict


