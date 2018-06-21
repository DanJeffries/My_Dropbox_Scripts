
# coding: utf-8

# In[139]:

float("1e-20")


# In[176]:

import sys
from Bio.Blast import NCBIXML
from Bio import SeqIO
  

def BlastParseExtra(infile, genome_fa, best_hit_criteria, Eval_thresh, window_size):
    """
    
    Usage: BlastParseExtra.py  <blast_xml_output> <genome_fasta> <best_hit_criteria> <Eval_thresh> <window_size> 
    
    <blast_xml_output>  -  absolute path to blast xml output
    <genome_fasta>      -  absolute path to genome fasta file
    <best_hit_criteria> -  factor by which the best e-value must be lower than the second best (recommended: 1e-5)
    <Eval_thresh>       -  e-value threshold for unique alignments
    <window_size>       -  size of the window (in bp) around the mapping coordinates used to extract the subject scaffold segment
    
    This script first filters the mappings in the <blast_xml_output> for uniq hits with evalues better than
    <Eval_thresh> or for multi hits where the best hit is <best_hit_criteria> orders of magnitude better than 
    the second.
    
    It will then retrieve a segment of the scaffold from <genome_fasta> which is + and - the <window_size> around
    the mapping coordinates for each query. If the ends of the scaffold are not within this window, then the 
    length of the segment will be (length of mapped query sequence + 2 x <window_size>). However if an end of a
    scaffold is within this window, the segment will be trimmed to this length.
    
    """


    blast = NCBIXML.parse(open(infile,"r"))

    good_blast_outs = {}
    multi_counter = 0
    unique_counter = 0
    ## From Alan's script: Returns blast hits only when the best e-value is 5 orders of magnitude better than the second best.

    for record in blast :
        if len(record.alignments)>1:
            if record.alignments[0].hsps[0].expect < best_hit_criteria * record.alignments[1].hsps[0].expect:
                multi_counter += 1
                good_blast_outs[record.query] = {}
                good_blast_outs[record.query]["Ref_hit_id"] = str(record.alignments[0].hit_def)
                good_blast_outs[record.query]["Evalue"] = float(record.alignments[0].hsps[0].expect)
                good_blast_outs[record.query]["Hit_start_coord"] = int(record.alignments[0].hsps[0].sbjct_start)
                good_blast_outs[record.query]["Hit_end_coord"] = int(record.alignments[0].hsps[0].sbjct_end)
                print "Multi\t%s\t%s\t%s\t%s\t%s" % (record.query, good_blast_outs[record.query]["Ref_hit_id"], good_blast_outs[record.query]["Evalue"], good_blast_outs[record.query]["Hit_start_coord"], good_blast_outs[record.query]["Hit_end_coord"])
                
        elif len(record.alignments)==1:
            if float(record.alignments[0].hsps[0].expect) < Eval_thresh:
                unique_counter += 1
                good_blast_outs[record.query] = {}
                good_blast_outs[record.query]["Ref_hit_id"] = str(record.alignments[0].hit_def)
                good_blast_outs[record.query]["Evalue"] = float(record.alignments[0].hsps[0].expect)
                good_blast_outs[record.query]["Hit_start_coord"] = int(record.alignments[0].hsps[0].sbjct_start)
                good_blast_outs[record.query]["Hit_end_coord"] = int(record.alignments[0].hsps[0].sbjct_end)
                print "Uniq\t%s\t%s\t%s\t%s\t%s" % (record.query, good_blast_outs[record.query]["Ref_hit_id"], good_blast_outs[record.query]["Evalue"], good_blast_outs[record.query]["Hit_start_coord"], good_blast_outs[record.query]["Hit_end_coord"])
    print "Number of multi-alingments kept:", multi_counter
    print "Number of unique alingments kept:", unique_counter
    
    
    Rtemp_chunks = open("%s/blast_%s_chunks.fa" % (infile.rpartition("/")[0], window_size), 'w')
    
    print "Getting subject scaffold segments from %s . . . " % (genome_fa)
    
    segment_counter = 0
    
    Rtemp = SeqIO.parse(genome_fa, "fasta")
    
    for scaffold in Rtemp:
        
        for query in good_blast_outs:
                
            if scaffold.id == good_blast_outs[query]["Ref_hit_id"]: ## If the scaffold has a hit
                                
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
                
    Rtemp_chunks.close()
    
    print "%s sequence scaffold segments are in %s/blast_%s_chunks.fa" % (segment_counter, infile.rpartition("/")[0], window_size)


# In[ ]:

## Cline args

if len(sys.argv) == 1:
    
    print BlastParseExtra.__doc__
    
elif len(sys.argv) < 6: ## If not enough args are supplied print error message
    sys.exit("\n##Error, not enough arguments, run script with no arguments to see help message\n")

elif len(sys.argv) == 6:
    in_path = sys.argv[1]
    gen_fa = sys.argv[2]
    multi_criteria = float(sys.argv[3])
    Eval = float(sys.argv[4])
    window = int(sys.argv[5])

    BlastParseExtra(in_path, gen_fa, multi_criteria, Eval, window) ## RUN!


# ### args for use without cline
# in_path = "/home/djeffrie/Data/Ribe_LM/Male_LG3_Rtemp_blasthits.xml"
# gen_fa = "/home/djeffrie/Data/Genomes/Rtemp/Rtemp_new.fa"
# multi_criteria = 1e-5
# Eval = 1e-20
# window = 1000
# 
# BlastParseExtra(in_path, gen_fa, multi_criteria, Eval, window)
