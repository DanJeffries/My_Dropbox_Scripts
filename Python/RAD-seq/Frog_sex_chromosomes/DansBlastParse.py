
# coding: utf-8
import sys
from Bio.Blast import NCBIXML
from Bio import SeqIO


def DansBlastParse(blast_xml_output, best_hit_criteria, Eval_thresh):
    """
    
    Usage: DansBlastParse.py  <blast_xml_output>  <best_hit_criteria>  <Evalue_threshold>
    
    <blast_xml_output>  -  Path should be absolute
    <best_hit_criteria> -  Orders of magnitude higher that best hit has to be in multi alignment (scientific format)
    <Evalue_threshold>  -  (1e-20 is default)
    
    
    This script first filters the mappings in the <blast_xml_output>, for only those that are uniq and map with an 
    e-value higher than <Evalue_threshold> or mappings where there are multiple hits but where the best hit is at 
    least 5 orders of magnitude higher than the second.
    
    
    """

    blast = NCBIXML.parse(open(blast_xml_output,"r"))

    good_blast_counts = 0
    uniq_counts = 0
    
    ## From Alan's script: Returns blast hits only when the best e-value is 5 orders of magnitude better than the second best.
    for record in blast :
        #print len(record.alignments)
        if len(record.alignments)==1:  ## Extra else statement - this wasn't in Alan's original script.
            if float(record.alignments[0].hsps[0].expect) < Eval_thresh:
                subject = str(record.alignments[0].hit_def)
                Evalue = float(record.alignments[0].hsps[0].expect)
                Hit_start_coord = int(record.alignments[0].hsps[0].sbjct_start)
                Hit_end_coord = int(record.alignments[0].hsps[0].sbjct_end)
                print "Unique\t%s\t%s\t%s\t%s\t%s" % (record.query, subject, Evalue, Hit_start_coord, Hit_end_coord)
                good_blast_counts += 1
                uniq_counts += 1
        
        elif len(record.alignments)>1:
            if all([record.alignments[0].hsps[0].expect < best_hit_criteria * record.alignments[1].hsps[0].expect, record.alignments[0].hsps[0].expect <= Eval_thresh]):
                
                subject = str(record.alignments[0].hit_def)
                Evalue = float(record.alignments[0].hsps[0].expect)
                Hit_start_coord = int(record.alignments[0].hsps[0].sbjct_start)
                Hit_end_coord = int(record.alignments[0].hsps[0].sbjct_end)
                print "Multi\t%s\t%s\t%s\t%s\t%s" % (record.query, subject, Evalue, Hit_start_coord, Hit_end_coord)
                good_blast_counts += 1

    print "Number of good blast hits:", good_blast_counts
    print "Number of unique mappings:", uniq_counts



## Cline args

if len(sys.argv) == 1:
    
    print DansBlastParse.__doc__
    
elif len(sys.argv) < 4: ## If not enough args are supplied print error message
    sys.exit("\n##Error, not enough arguments, run script with no arguments to see help message\n")

elif len(sys.argv) == 4:
    in_path = sys.argv[1]
    muti_hit_criteria = float(sys.argv[2])
    Eval = float(sys.argv[3])


    DansBlastParse(in_path, muti_hit_criteria, Eval) ## RUN!


