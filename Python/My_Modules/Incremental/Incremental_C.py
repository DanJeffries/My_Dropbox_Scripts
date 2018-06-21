import sys
import os
import subprocess
import matplotlib
import re
import Incremental_U as Inc
import pprint
import pylab
import gzip
import csv
import random
import numpy as np
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt
from matplotlib.patches import Circle, Ellipse
from itertools import chain
from collections import Iterable

## Plotting parameters ##

plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams['figure.figsize'] = (20.0, 20.0)


def natural_key(string_): 
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]



def MakeAndRun_Cstacks_cline(Setup_dict):

    Cstacks_test_dir = Setup_dict["Ustacks_outs_dir"] ## change if you want them to go somewhere else

    ### Making the Cstacks command line and getting IDs
    samples_in_directory = {}
    sample_cline = []
    for sample in os.listdir(Setup_dict["Ustacks_outs_dir"]):
        if "tags" in sample:
            
            sample_prefix = sample.split(".tags")[0]
            
            if sample.endswith("gz"):
                f = gzip.open("%s/%s" % (Setup_dict["Ustacks_outs_dir"],sample), 'r')
            else:
                f = open("%s/%s" % (Setup_dict["Ustacks_outs_dir"],sample), 'r')
            
            for line in f.readlines()[:10]:
                if "consensus" in line:
                    sample_ID = line.split()[1]
                    samples_in_directory[sample_prefix] = sample_ID
            
            if len(Setup_dict["sample_ids"]) > 0:
                if int(sample_ID) in Setup_dict["sample_ids"]:
                    sample_cline.append("-s %s/%s" % (Setup_dict["Ustacks_outs_dir"],sample_prefix))
                    
    samples = " ".join(sample_cline)

    for n_value in Setup_dict["n_vals"]:
        outdir = "%s/n_%s" % (Cstacks_test_dir, n_value)

        if not os.path.exists(outdir):
            try:
                os.makedirs(outdir)
            except:
                os.mkdir(outdir)

        print "\nn value = %s" % n_value
        Cstacks_cline = "cstacks -b %s -n %s %s %s -o %s" % (Setup_dict["batch_ID"], n_value, samples, Setup_dict["threads"], outdir)

        print Cstacks_cline

        if Setup_dict["run_cline_switch"] == 1:
            print "Running Command line"
            log = open("%s/Cstacks.log" % outdir, 'w')
            subprocess.call(Cstacks_cline.split(), stdout=log, stderr=log )
            
    return samples_in_directory

## -------------------------------------------------------------------------------------------------------

def Get_sample_IDs(directory, verbose = True):
    """ 
    Gets sample IDs from tags.tsv files in a directory
    
    """             
    sample_info = {}
    for root, dirs, samples in os.walk(directory):
        for sample in samples:
            if "tags.tsv" in sample and "catalog" not in sample:

                sample_prefix = sample.split(".")[0]

                if not sample_prefix in sample_info:
                    if sample.endswith("gz"):
                        sample_file = gzip.open("%s/%s" % (root, sample))
                    else:
                        sample_file = open("%s/%s" % (root, sample))
                    for line in sample_file.readlines()[:10]:
                        if "consensus" in line:
                            sample_ID = line.split()[1]

                            sample_info[sample_prefix] = sample_ID
                            
    if verbose == True:
        for sample in sorted(sample_info.keys()):
            print "%s, ID = %s" % (sample, sample_info[sample])
        
    return sample_info


#--------------------------------------------------------------------
alignment = {'horizontalalignment':'center', 'verticalalignment':'center'}

#--------------------------------------------------------------------
def get_labels(data, fill="number"):
    """
    to get a dict of labels for groups in data
    input
      data: data to get label for
      fill = ["number"|"logic"|"both"], fill with number, logic label, or both
    return
      labels: a dict of labels for different sets
    example:
    In [12]: get_labels([range(10), range(5,15), range(3,8)], fill="both")
    Out[12]:
    {'001': '001: 0',
     '010': '010: 5',
     '011': '011: 0',
     '100': '100: 3',
     '101': '101: 2',
     '110': '110: 2',
     '111': '111: 3'}
    """

    N = len(data)

    sets_data = [set(data[i]) for i in range(N)]  # sets for separate groups
    s_all = set(chain(*data))                             # union of all sets

    # bin(3) --> '0b11', so bin(3).split('0b')[-1] will remove "0b"
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i] for i in range(N) if  key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if  key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value

    if fill == "number":
        labels = {k: len(set_collections[k]) for k in set_collections}
    elif fill == "logic":
        labels = {k: k for k in set_collections}
    elif fill == "both":
        labels = {k: ("%s: %d" % (k, len(set_collections[k]))) for k in set_collections}
    else:  # invalid value
        raise Exception("invalid value for fill")

    return labels



def venn4(data=None, names=None, fill="number", show_names=True, show_plot=True, axes = None, **kwds):

    if (data is None) or len(data) != 4:
        raise Exception("length of data should be 4!")
    if (names is None) or (len(names) != 4):
        names = ("set 1", "set 2", "set 3", "set 4")

    labels = get_labels(data, fill=fill)

    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else: # default figure size
        figsize = (10, 10)

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 4:
        colors = kwds['colors']
    else:
        colors = ['r', 'g', 'b', 'c']

    # draw ellipse, the coordinates are hard coded in the rest of the function
    
    patches = []
    width, height = 170, 110  # width and height of the ellipses
    patches.append(Ellipse((170, 170), width, height, -45, color=colors[0], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -45, color=colors[1], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -135, color=colors[2], alpha=0.5))
    patches.append(Ellipse((230, 170), width, height, -135, color=colors[3], alpha=0.5))
    for e in patches:
        axes.add_patch(e)
    axes.set_xlim(80, 320); axes.set_ylim(80, 320)
    axes.set_xticks([]); axes.set_yticks([]);
    axes.set_aspect("equal")

    ### draw text
    # 1
    pylab.text(120, 200, labels['1000'], fontsize=10, **alignment)
    pylab.text(280, 200, labels['0100'], fontsize=10, **alignment)
    pylab.text(155, 250, labels['0010'], fontsize=10, **alignment)
    pylab.text(245, 250, labels['0001'], fontsize=10, **alignment)
    # 2
    pylab.text(200, 115, labels['1100'],fontsize=10,  **alignment)
    pylab.text(140, 225, labels['1010'],fontsize=10,  **alignment)
    pylab.text(145, 155, labels['1001'],fontsize=10,  **alignment)
    pylab.text(255, 155, labels['0110'],fontsize=10,  **alignment)
    pylab.text(260, 225, labels['0101'],fontsize=10,  **alignment)
    pylab.text(200, 240, labels['0011'],fontsize=10,  **alignment)
    # 3
    pylab.text(235, 205, labels['0111'],fontsize=10,  **alignment)
    pylab.text(165, 205, labels['1011'],fontsize=10,  **alignment)
    pylab.text(225, 135, labels['1101'],fontsize=10,  **alignment)
    pylab.text(175, 135, labels['1110'],fontsize=10,  **alignment)
    # 4
    pylab.text(200, 175, labels['1111'],fontsize=10,  **alignment)
    # names of different groups
    if show_names:
        pylab.text(110, 110, names[0], fontsize=16, **alignment)
        pylab.text(290, 110, names[1], fontsize=16, **alignment)
        pylab.text(130, 275, names[2], fontsize=16, **alignment)
        pylab.text(270, 275, names[3], fontsize=16, **alignment)

        
    #leg_names = [names[0], names[2], names[3], names[1]]    ## the legend function messes up the names order so gave it manually here
    #leg = ax.legend(leg_names, loc='best', fancybox=True)
    #leg.get_frame().set_alpha(0.5)


    return axes
    
    
    #if not outdir == False and not n_val == False:
     #   pylab.savefig("%s/Venn4_n=%s" % (outdir, n_val))
    #
    #if show_plot:
     #   pylab.show()


def Tag_share(directory, sample_ID_list = [], n_vals = []):
    
    print "Number of samples specified: %s\n" % len(sample_ID_list)
    
    ## Adjust global fonts to make figures a bit neater

    font = {'family' : 'sans','weight' : 'bold','size'   : 4}
    matplotlib.rc('font', **font)

    ## Define lists
    
    Sample1 = []
    Sample2 = []
    Sample3 = []
    plot_y = []
    shared = []
        
    ## Get the catalog files
    
    catalog_files = []
    for root, dirs, files in os.walk(directory):
        for fil in files:
            if 'batch' in fil and 'tags' in fil:
                cat_nval = root.rpartition("/")[2].split("_")[1]
                if int(cat_nval) in n_vals:
                    catalog_files.append("%s/%s" % (root, fil))

    catalog_files = sorted(catalog_files, key = natural_key)    
    
    ## Get locus sharing data for each catalog.
    catalogs_processed = 0
    
    fig = plt.figure(figsize= (20,len(catalog_files)*5))
    #fig.subplots_adjust(hspace=0.5)
    
    plot_number = 1
    total_shared = {}
    
    for catalog in catalog_files:
        
        cat_id = catalog.rpartition("/")[0].rpartition("/")[2]
        
        Sample_set_1 = []
        Sample_set_2 = []
        Sample_set_3 = []
        Sample_set_4 = []
        
        f = gzip.open(catalog, 'rb') 
        csv_read = csv.reader(f,delimiter ="\t")
        seqNumber = 1 ## Start seqNumber
        
        Loc_dict = {}
        sample_ids = []
        
        for line in csv_read:
            
            if "consensus" in line:
                loc_ID = line[2]
                Loc_dict[loc_ID] = []
                    
                locfield = line[8] ## Isolate and split the relevant column
                    
                seqNumber += 1

                ## Assign the seqNumber to a list depending on which individuals it is present in

                for sample in locfield.split(","):
                    sample_ID = sample.split("_")[0]
                        
                    if sample_ID not in sample_ids:
                        sample_ids.append(sample_ID) ## list of all sample IDs in catalog
                        
                    Loc_dict[loc_ID].append(sample_ID)
                    #print loc_ID, sample, sample_ID

                Loc_dict[loc_ID] = set(Loc_dict[loc_ID])
        
        
        ### Get Venn data --------------------------------------------------------------------------------------
        
        if len(sample_ID_list) == 0:  ## if no list specified
            
            if len(sample_ids) < 4:
                N_samples_to_get = len(sample_ids)
            else:
                N_samples_to_get = 4
            

            ## Get 4 sample sets from the samples present. Or take all samples if there are fewer than 4.
            if catalogs_processed == 0: ## this means that it will only get random IDs once, to keep the samples used the same in all plots.
            
                rand_smpls = [ sample_ids[i] for i in sorted(random.sample(xrange(len(sample_ids)), N_samples_to_get))]
                
                sample_ID_list = rand_smpls
                
            catalogs_processed += 1      
            
            Venn_sets = {}
            
            for sample in rand_smpls:
                Venn_sets[sample] = []
                for locus in Loc_dict:
                    if sample in Loc_dict[locus]:
                        Venn_sets[sample].append(locus)
            
                
        elif len(sample_ID_list) == 1:
            sys.exit("Only specified one sample ID, need at least two")
            
        elif len(sample_ID_list) > 4:
            sys.exit("Specified too many sample IDs, need at most 4")
        
        elif len(sample_ID_list) > 1: ## If list was specified

            Venn_sets = {}
            
            for sample in sample_ID_list:
                sample = str(sample)
                
                if sample not in sample_ids: ## check that IDs are in there. 
                    sys.exit("Specified sample ID (%s) not in catalog" % sample)
                else:
                    Venn_sets[sample] = []
                    for locus in Loc_dict:
                        if sample in Loc_dict[locus]:
                            Venn_sets[sample].append(locus)
                        
                        
        
        ## plot venns

        if len(Venn_sets) == 2:
            
            fig.add_subplot(len(catalog_files)/2+1,2,plot_number)
            lab1, vals1 = Venn_sets.items()[0]
            lab2, vals2 = Venn_sets.items()[1]
            
            shared_all = set(vals1) & set(vals2)
            total_shared[cat_id] = len(shared_all)
            
            v2 = venn2([set(vals1), set(vals2)], set_labels= (str(lab1), str(lab2)), normalize_to=3.0)
            
            for text in v2.set_labels:
                text.set_fontsize(15)
            for text in v2.subset_labels:
                text.set_fontsize(10)
            
            plt.title("%s, Total tags = %s" % (catalog.rpartition("/")[0].rpartition("/")[2],seqNumber), fontsize = 15)
            plot_number += 1
            
        elif len(Venn_sets) == 3:
            
            fig.add_subplot(len(catalog_files)/2+1,2,plot_number)
            lab1, vals1 = Venn_sets.items()[0]
            lab2, vals2 = Venn_sets.items()[1]
            lab3, vals3 = Venn_sets.items()[2]
            
            shared_all = set(vals1) & set(vals2) & set(vals3)
            total_shared[cat_id] = len(shared_all)
            
            v3 = venn3([set(vals1), set(vals2), set(vals3)], set_labels= ("ID=%s" % str(lab1), "ID=%s" % str(lab2), "ID=%s" % str(lab3)))
            for text in v3.set_labels:
                text.set_fontsize(15)
            for text in v3.subset_labels:
		if not text == None:
		    text.set_fontsize(10)
            
            plt.title("%s, Total tags = %s" % (catalog.rpartition("/")[0].rpartition("/")[2],seqNumber), fontsize = 15)
            plot_number += 1
            
        elif len(Venn_sets) == 4:
            
            lab1, vals1 = Venn_sets.items()[0]
            lab2, vals2 = Venn_sets.items()[1]
            lab3, vals3 = Venn_sets.items()[2]
            lab4, vals4 = Venn_sets.items()[3]
            
            shared_all = set(vals1) & set(vals2) & set(vals3) & set(vals4)
            total_shared[cat_id] = len(shared_all)
            
            ax = pylab.subplot(len(catalog_files)/2+1,2,plot_number)
            v4 = venn4([set(vals1), set(vals2), set(vals3), set(vals4)], names= ["ID=%s" % str(lab1), "ID=%s" % str(lab2), "ID=%s" % str(lab3), "ID=%s" % str(lab4)], axes = ax) 
            ax.set_title("%s, Total tags = %s" % (catalog.rpartition("/")[0].rpartition("/")[2],seqNumber), fontsize = 15)
            plot_number += 1
            
            print "N tags in catalog %s = %s," % (catalog.rpartition("/")[0].rpartition("/")[2], seqNumber), "Number of loci shared by all = %s" % len(shared_all)

    
    
    ## for the shared by all plot
    counter = range(1,100,1)
    y = []
    xlabs = []
    
    for key, val in total_shared.items():
        for i in counter:
            if key == "n_%s" % i:
                y.append(val)
                xlabs.append(i)
    
    fig.add_subplot(len(catalog_files)/2+1,2,plot_number)
    plt.bar(xlabs, y , color = "blue", alpha = 0.4, align = "center", width = 0.6)
    plt.ylim(np.mean(y)-(np.std(y)*2), np.mean(y)+(np.std(y)*2))
    plt.xticks(xlabs, [str(i) for i in xlabs])
    plt.xlabel("n value", fontsize = 10)
    plt.title("Loci shared between all samples in venn", fontsize = 12)
    plt.savefig("%s/Tag_sharing_IDs_%s.pdf" % (directory, "_".join([str(i) for i in sample_ID_list])))
    plt.show()
    plt.close()


    
def Incremental_C(Analysis_parameters):
    
    
    """
    Structure of the Analyses_parameter dictionary is:
   
    Analysis_parameters = {}
    Analysis_parameters["batch_ID"] = 1
    Analysis_parameters["n_vals"] = [1,2,3,4]
    Analysis_parameters["threads"] = 7
    Analysis_parameters["Ustacks_outs_dir"] = "/home/djeffrie/Data/RADseq/Rdal_test/Cstacks_tests_N3/"
    Analysis_parameters["sample_ids"] = []   ## add specific sample ID's of samples you want to look at, if this is not all samples in the directory given above.
    Analysis_parameters["run_cline_switch"] = 0     ## to switch on (1) to run stacks cline or off (0) to run pipeline on existing outputs
    Analysis_parameters["Get_IDs"] = True ## To get IDs. Will not continue with pipeline.
    """
    
    if Analysis_parameters["Get_IDs"] == True:
        print "\n ## Sample IDs:\n"
    
        ids = Get_sample_IDs(Analysis_parameters["Ustacks_outs_dir"])
    else:
        print "\n### ----- Making and running command lines ----- ###"

        samples = MakeAndRun_Cstacks_cline(Analysis_parameters)

        print "\n\n#### ---- Plotting the tag sharing patterns between individuals ---- ####\n"

        Tag_share(Analysis_parameters["Ustacks_outs_dir"], Analysis_parameters["sample_ids"], Analysis_parameters["n_vals"])
