from __future__ import division
import sys
import vcf
from collections import Counter
import pprint as pp
from pydoc import help


# In[1]:


def MST_input_maker(vcf_path, sex_info_file, offspring_presence_threshold, mendelian_threshold):
    
    """
    Usage: MST_input_maker  <vcf_path>  <sex_info_file>  <offspring_presence_threshold>  <mendelian_threshold>
    
    ### Informative locus filters:
    -Loci are used only if the number of offspring present is equal to or above the percentage specified by the user
     e.g. a threshold of 0.75 will only keep loci which are present in 75% of offspring or more.
    -Loci are used only if they are heterozygous in one sex and homozygous in another.
    
    ###Null allele filters:

    Null allele in Parents:
    -If an offspring is found to have an allele that isn't present in the parent, this locus is discarded

    Null allele in offspring:
    -If only 1 offspring is homozygous for the parental minor allele then that sample is assumed to possess a null 
     allele and the data for that individual is coded as missing. **Note, the sample presence/absence filter works 
     after this step, so some loci may be pushed over the "missing data" limit by this null allele filter. 
    -If two or more samples are found to be homozygous for the minor parental allele, it is more likely that this 
     is due to the homozygous parent having a null allele. In this case, the locus is again discarded.
    -Finally, loci with non-mendelian segregation proportions indicative of allele dropout or similar are also filtered 
     using the <mendelian_threshold> set by the user. It is recommended that this threshold be set to 0.75. This means
     that, if a site is heterozygous in a parent and more than 75% of offspring are found to be either heterozygous or 
     homozygous for the minor allele then this locus is discarded. 
    
    ### sex_info_file format example

	Sample_1	Male
	Sample_2	Female
	Sample_3	Offspring
	Sample_4	Offspring
	Sample_5	Offspring
	...

     The 2nd column labels can be:
	M, m, Male, male, MALE, Father, father (for the male)
	F, f, Female, female, FEMALE, Mother, mother (for the female)
	O, o, Offspring, offspring, OFFSPRING (for the offspring)

     """
    

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
    
    myvcf = vcf.Reader(altered_vcf)
        
    MST_path = "%s/" % (vcf_path.rpartition("/")[0])
    
    loci_used_for_male_map = 0
    loci_used_for_female_map = 0
    
    Loc_w_null_alleles = 0
    
    offspring = []
    
    femMSTfile = open("%sMST_female_map_input.txt" % (MST_path), 'w')
    malMSTfile = open("%sMST_male_map_input.txt" % (MST_path), 'w')
    log_file = open("%sMST_input_maker.log" % (MST_path), 'w')


    sex_file = open(sex_info_file, 'r').readlines()

    for line in sex_file: ## Note, only samples in the sex_info file are used
        
        sample_name = line.split()[0]
        sex = line.split()[1]                    
        
        if sex == "F" or sex == "f" or sex == "Female" or sex == "female" or sex == "FEMALE" or sex == "Mother" or sex == "mother" or sex == "MOTHER":
            Mother = sample_name
            log_file.write("Mother = %s" % (sample_name))
        
        elif sex == "M" or sex == "m" or sex == "Male" or sex == "male" or sex == "MALE" or sex == "Father" or sex == "father" or sex == "FATHER" :
            Father = sample_name
            log_file.write("\nFather = %s\n" % (sample_name))        
    
        elif sex == "O" or sex == "o" or sex == "Offspring" or sex == "offspring" or sex == "OFFSPRING":
            offspring.append(sample_name)        
        
    offspring = sorted(offspring) ## sort the list of offspring - used to retain correct order throughout script
    
    
    femmap_file_lines = []
    malmap_file_lines = []
    
    ## Counters for log_file summary information:
    locus_counter = 0 
    loci_w_too_many_missing = 0
    loci_w_data_missing_in_par = 0
    loci_w_both_par_hom = 0
    loci_w_both_par_het = 0
    loci_w_too_many_par_alleles = 0
    loc_w_excess_hom = 0    
    loci_w_non_par_allele_in_offspr = 0
    loci_w_multi_offspr_hom_for_par_minor_allele = 0


    for record in myvcf:
        locus_counter += 1
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
    
        if record.ID == None:
            if record.CHROM == None:
                Loc_Id = "%s_%s" % (locus_counter, record.POS) ## If ID and CHROM not in vcf, then just use locus order.
                ID_line = "\nUsed locus order and position as locus IDs"
            else:
                Loc_Id = "%s_%s" % (record.CHROM, record.POS) ## if ID not in vcf, use CHROM as the major ID. 
                ID_line = "\nUsed scaffold (CHROM) and position as locus IDs"
        else:        
            Loc_Id = "%s_%s" % (record.ID, record.POS) ## If record ID and POS there, use both as the ID
            ID_line = "\nNOTE: Used marker ID (ID) and position as locus IDs\n"

    
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
            loci_w_data_missing_in_par += 1
        elif len(Par_allele_counts) > 2: ## If there are more than two alleles in the parents, discard locus
            log_file.write("\nparents contain more than 2 alleles, locus not used")
            loci_w_too_many_par_alleles += 1
        elif pat_het == True and mat_het == True: ## If both parents are heterozygous, discard locus
            log_file.write("\nBoth parents heterozygous, locus not used")
            loci_w_both_par_het += 1
        elif pat_het == False and mat_het == False: ## If both parents are homozygous, discard locus
            log_file.write("\nBoth parents homozygous, locus not used")
            loci_w_both_par_hom += 1
        elif len(Par_allele_counts) == 2 and not 2 in Par_allele_counts.values():
            for allele in Par_allele_counts:
                if Par_allele_counts[allele] == 3:
                    Par_major_allele = allele
                elif Par_allele_counts[allele] == 1:
                    Par_minor_allele = allele
        
        ## Get offspring information -----------------------------------------------------------------------
               
            for sample in record.samples:            
                name = sample.sample
                genotype = sample['GT']            
        
                if name in offspring:  ## Only samples included in the pop_map file are used
                    N_offspring += 1
                        
                    ### Filter sample genotypes for signs of allele droupout 
                                    
                    if genotype == None or genotype == "./.":                    
                            
                        offspring_gt_types.append(sample.gt_type) # record hom or het
                        offspring_gt_bases_dict[name] = (sample.gt_type, sample.gt_bases)
                                            
                    elif genotype == Mat_genotype or genotype == Pat_genotype or genotype == Mat_genotype[::-1] or genotype == Pat_genotype[::-1]:
                            
                        offspring_gt_types.append(sample.gt_type) # record hom or het
                        offspring_gt_bases_dict[name] = (sample.gt_type, sample.gt_bases)
                                
                    elif sample.gt_bases == "%s/%s" % (Par_minor_allele,Par_minor_allele):
                        log_file.write("\n%s is homozygous for minor parental allele (%s), likely allele dropout in offspring (if rare) or parent (if common)" % (name, sample.gt_bases))
                        offspr_w_allele_dropout += 1
                            
                        offspring_gt_types.append(None) # Count as missing data
                        offspring_gt_bases_dict[name] = (None, None)
                        N_offsp_hom_for_minor_par_allele += 1
                            
                    else:
                        off_ref = sample.gt_bases.split("/")[0]
                        off_alt = sample.gt_bases.split("/")[1]
                                
                        if off_ref not in Parental_alleles or off_alt not in Parental_alleles:
                            offspring_gt_types.append(sample.gt_type) # record hom or het
                            offspring_gt_bases_dict[name] = (sample.gt_type, sample.gt_bases)
                                
                            log_file.write("\n%s contains non-parental allele in genotype (%s), likely allele dropout in a Parent" % (name, sample.gt_bases))
                            Non_par_allele = True
                
            counted = Counter(offspring_gt_types) ## Count the numbers of each genotype at this locus
            
            ## Caluculate the percentage of offspring missing at the locus
            
            if None in counted:
                perc_offspring_missing = counted[None]/N_offspring
            else:
                perc_offspring_missing = 0
            
            ## calculate the percentages of homozygotes/heterozygotes for mendelian frequency filters below
            
            if (counted[0]+counted[2]+counted[1]) > 0:
		
            
                perc_hom_0 = counted[0]/(counted[0]+counted[2]+counted[1])
                perc_hom_2 = counted[2]/(counted[0]+counted[2]+counted[1])
                perc_het = counted[1]/(counted[0]+counted[2]+counted[1])
                perc_off_w_all_dropout = offspr_w_allele_dropout/N_offspring
                      
            
                ## Filter loci with non-mendelian genotype proportions and write data to repsective map files
            
                if Non_par_allele == True:
                    log_file.write("\nNon parental allele found, locus not used")
                    loci_w_non_par_allele_in_offspr += 1
                elif N_offsp_hom_for_minor_par_allele > 1: ## if there is more than one offspring homozygous for the minor parental allele, then it is likely to be allele dropout in the parent, locus is discarded
                    log_file.write("\nMore than one offspring homozygous for minor allele, locus not used")
                    loci_w_multi_offspr_hom_for_par_minor_allele += 1
                elif perc_hom_0 > mendelian_threshold or perc_hom_2 > mendelian_threshold:
                    log_file.write("\nHomozygosity excess (>%s)" % (mendelian_threshold))
                    log_file.write("\nLocus discarded, Perc_hom_0= %s, Perc_hom_2= %s" % (perc_hom_0, perc_hom_2))
                    loc_w_excess_hom += 1
                elif perc_offspring_missing > 1-offspring_presence_threshold: ## put this filter last to allow other filters to work first and give accurate error report
                    log_file.write("\nToo many missing genotypes (%s), locus not used" % (perc_offspring_missing))
                    loci_w_too_many_missing += 1
    
                # Get data for Female map--------------------------------------------------------------------------------------------------------               
                    
                elif mat_het == True and pat_het == False:
                    log_file.write("\n\n###Used in Female map###\n")
                    loci_used_for_female_map += 1
                    femalemap_offspring_haploid_genotypes.append("\n%s" % (Loc_Id)) ## make line for the record (not used unless criteria below are met)
                    femalemap_offspring_haploid_genotypes_compli.append("\ncompli_%s" % (Loc_Id)) ## make line for the complimentary record (not usd unless criteria below are met)
                
                    for i in offspring: ## retain order throughout
                        if i in sample_names: ## make sure all samples in the pop_codes file are in the vcf
                            if offspring_gt_bases_dict[i][1] == Mat_genotype_bases:
                                off_MST = "B"
                                off_compli = "a"
                                femalemap_offspring_haploid_genotypes.append("B")
                                femalemap_offspring_haploid_genotypes_compli.append("a")
                            elif offspring_gt_bases_dict[i][1] == Pat_genotype_bases:
                                off_MST = "A"
                                off_compli = "b"
                                femalemap_offspring_haploid_genotypes.append("A")
                                femalemap_offspring_haploid_genotypes_compli.append("b")
                            elif offspring_gt_bases_dict[i][1] == None:
                                off_MST = "U"
                                off_compli = "U"
                                femalemap_offspring_haploid_genotypes.append("U")
                                femalemap_offspring_haploid_genotypes_compli.append("U")
                            else:
                                log_file.write("\nNo genotype match, Mat_gt = %s, Pat_gt = %s, off_gt = %s" %(Mat_genotype_bases, Pat_genotype_bases, offspring_gt_bases_dict[i][1]))
                                femalemap_offspring_haploid_genotypes.append("U")
                                femalemap_offspring_haploid_genotypes_compli.append("U")
                                
                    femmap_file_lines.append("\t".join(femalemap_offspring_haploid_genotypes))
                    femmap_file_lines.append("\t".join(femalemap_offspring_haploid_genotypes_compli))
                    
                    log_file.write("\n\nPercentage homozygous offspring for reference allele: %s" % (perc_hom_0))
                    log_file.write("\nPercentage homozygous offspring for alternative allele: %s" % (perc_hom_2))
                    log_file.write("\nPercentage heterozygous offspring for reference allele: %s" % (perc_het))
                    log_file.write("\nPerc offspring with probable allele dropout: %s" % (offspr_w_allele_dropout/N_offspring))
		
                            
                # Do the same to get data forMale map----------------------------------------------------------------------------------------------                
            
                elif pat_het == True and mat_het == False: 
                    log_file.write("\n\n###Used in Male map###\n")
                    loci_used_for_male_map += 1
                    malemap_offspring_haploid_genotypes.append("\n%s" % (Loc_Id)) ## make line for the record (not used unless criteria below are met)
                    malemap_offspring_haploid_genotypes_compli.append("\ncompli_%s" % (Loc_Id)) ## make line for the complimentary record (not used unless criteria below are met)            
                        
                    for i in offspring:
                        if i in sample_names:
                            if offspring_gt_bases_dict[i][1] == Mat_genotype_bases:
                                off_MST = "A"
                                off_compli = "b"
                                malemap_offspring_haploid_genotypes.append("A")
                                malemap_offspring_haploid_genotypes_compli.append("b")
                            elif offspring_gt_bases_dict[i][1] == Pat_genotype_bases:
                                off_MST = "B"
                                off_compli = "a"
                                malemap_offspring_haploid_genotypes.append("B")
                                malemap_offspring_haploid_genotypes_compli.append("a")
                            elif offspring_gt_bases_dict[i][1] == None:
                                off_MST = "U"
                                off_compli = "U"
                                malemap_offspring_haploid_genotypes.append("U")
                                malemap_offspring_haploid_genotypes_compli.append("U")
                            else:
                                log_file.write("\nNo genotype match, Mat_gt = %s, Pat_gt = %s, off_gt = %s" %(Mat_genotype_bases, Pat_genotype_bases, offspring_gt_bases_dict[i][1]))
                                malemap_offspring_haploid_genotypes.append("\tU")
                                malemap_offspring_haploid_genotypes_compli.append("\tU")
                    
                    malmap_file_lines.append("\t".join(malemap_offspring_haploid_genotypes))
                    malmap_file_lines.append("\t".join(malemap_offspring_haploid_genotypes_compli))
                    
                
                    log_file.write("\n\nPercentage homozygous offspring for reference allele: %s" % (perc_hom_0))
                    log_file.write("\nPercentage homozygous offspring for alternative allele: %s" % (perc_hom_2))
                    log_file.write("\nPercentage heterozygous offspring for reference allele: %s" % (perc_het))
                    log_file.write("\nPerc offspring with probable allele dropout= %s" % (offspr_w_allele_dropout/N_offspring))
            else:       
                log_file.write("\n\nNo valid genotypes in the offspring after filters")

    log_file.write("\n\nSUMMARY:\nN Loci in VCF: %s" % (locus_counter))
    log_file.write("\nN loci used in female map: %s" % (loci_used_for_female_map))
    log_file.write("\nN loci used in male map: %s" % (loci_used_for_male_map))
    
    
    ## Make file headers
    N_offspring = 0
    Sample_headers = []        
    # Sample_headers.append("locus_name")
    for i in offspring:
        if i in sample_names:
            if len(Sample_headers) == 0:
                N_offspring += 1
                Sample_headers.append("%s" % (i))
            elif len(Sample_headers) > 0:
                N_offspring += 1
                Sample_headers.append("\t%s" % (i))
    Sample_header_line = ''.join(Sample_headers)
    
    
    #fem_header = "    population_type DH\n    population_name Female_map\n    distance_function kosambi\n    cut_off_p_value 0.000005\n    no_map_dist 30\n    no_map_size 1\n    missing_threshold 0.25\n    estimation_before_clustering no\n    detect_bad_data yes\n    objective_function ML\n    number_of_loci %s\n    number_of_individual %s\n\n" % (loci_used_for_female_map*2, N_offspring)
    
   # male_header = "    population_type DH\n    population_name Female_map\n    distance_function kosambi\n    cut_off_p_value 0.000005\n    no_map_dist 30\n    no_map_size 1\n    missing_threshold 0.25\n    estimation_before_clustering no\n    detect_bad_data yes\n    objective_function ML\n    number_of_loci %s\n    number_of_individual %s\n\n" % (loci_used_for_male_map*2, N_offspring)
    
    #femMSTfile.write(fem_header)
    femMSTfile.write(Sample_header_line)
    
    for line in femmap_file_lines:
        femMSTfile.write(line)
    
    #malMSTfile.write(male_header)
    malMSTfile.write(Sample_header_line)
    
    for line in malmap_file_lines:
        malMSTfile.write(line)
        
        
    femMSTfile.close()
    malMSTfile.close()
    log_file.close()
    
    print "Loci in VCF: %s" % (locus_counter)               
    print "N loci used in female map: %s" % (loci_used_for_female_map)
    print "N loci used in male map: %s" % (loci_used_for_male_map)
    print "N loci with too much missing data: %s (%1.2f percent)" % (loci_w_too_many_missing, (loci_w_too_many_missing/locus_counter)*100)
    print "N loci with missing data in one or both parents: %s (%1.2f percent)" % (loci_w_data_missing_in_par, (loci_w_data_missing_in_par/locus_counter)*100)
    print "N loci with two homozygous parents: %s (%1.2f percent)" % (loci_w_both_par_hom, (loci_w_both_par_hom/locus_counter)*100)
    print "N loci with two heterozygous parents: %s (%1.2f percent)" % (loci_w_both_par_het, (loci_w_both_par_het/locus_counter)*100)
    print "N loci with more than two different alleles in parents: %s (%1.2f percent)" % (loci_w_too_many_par_alleles, (loci_w_too_many_par_alleles/locus_counter)*100) ## Check this!!
    print "N loci with homozygosity excess in offspring: %s (%1.2f percent)" % (loc_w_excess_hom, (loc_w_excess_hom/locus_counter)*100)
    print "N loci with allele in offspring not found in parents: %s (%1.2f percent)" % (loci_w_non_par_allele_in_offspr, (loci_w_non_par_allele_in_offspr/locus_counter)*100)
    print "N loci where more than one offspring is homozygous for minor parental allele: %s (%1.2f percent)" % (loci_w_multi_offspr_hom_for_par_minor_allele, (loci_w_multi_offspr_hom_for_par_minor_allele/locus_counter)*100)
    print ID_line

# In[2]:

## Cline args

if len(sys.argv) == 1:
    # sys.exit(help(VCF_MSTmap_converter))
    print MST_input_maker.__doc__
    
elif len(sys.argv) < 5: ## If not enough args are supplied print error message
    sys.exit("\n##Error, not enough arguments, run script with no arguments to see help message\n")

elif len(sys.argv) == 5:
    vcf_file = sys.argv[1]
    sex_info = sys.argv[2]
    offsp_pres = float(sys.argv[3])
    mend_thresh = float(sys.argv[4])

    MST_input_maker(vcf_file, sex_info, offsp_pres, mend_thresh)


# In[ ]:



