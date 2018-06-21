from Bio import SeqIO
import sys
import gzip

### Progress bar!

def printProgress (iteration, total, prefix = '', suffix = '', decimals = 1, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    formatStr       = "{0:." + str(decimals) + "f}"
    percents        = formatStr.format(100 * (iteration / float(total)))
    filledLength    = int(round(barLength * iteration / float(total)))
    bar             = '|' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()



print "\nUsage: Filtering_reads_from_fasta.py <readpool.fast(a/q/gz)> <read_IDs.txt> [N_sequences_in_file (optional)]"
print "\n***Note: all paths must be absolute"
## CLINE options ====================================================

if len(sys.argv) < 3:
    sys.exit("\n#### ERROR: Not enough arguments####\n")

elif len(sys.argv) >= 3:
    fasta = sys.argv[1]
    Ids = open(sys.argv[2], 'r').readlines()

if len(sys.argv) == 4:
    N_seqs = sys.argv[3]
    p_bar = True
else:
    p_bar = False

## Open reads file ===================================================
    
if fasta.endswith("a"):
    fasta_it = SeqIO.parse(fasta, "fasta")
    file_format = "fasta"

elif fasta.endswith("q"):
    fasta_it = SeqIO.parse(fasta, "fastq")
    file_format = "fastq"

elif fasta.endswith("a.gz"):
    fasta_handl = gzip.open(fasta, 'r')
    fasta_it = SeqIO.parse(fasta_handl, "fasta")
    file_format = "fasta"

elif fasta.endswith("q.gz"):
    fasta_handl = gzip.open(fasta, 'r')
    fasta_it = SeqIO.parse(fasta_handl, "fastq")
    file_format = "fastq"

else:
    sys.exit("\nERROR: File format no recognised, expecting .fa, .fq, .fasta, .fastq, or any of the gzipped equivalients (.gz)\n")
    
## Open output file ==================================================

outfile_path = "%s/Mt_filtered_reads.%s" % (fasta.rpartition("/")[0],file_format)

filtered_fasta = open(outfile_path, 'w')

## Write all reads not mapped to output file =========================

print "\nWriting un-filtered reads to file\n"
if p_bar == True:
    i = 0
    l = N_seqs
    printProgress(i, l, prefix = 'Progress:', suffix = 'Complete', barLength = 50)

for record in fasta_it:
    if record.id not in Ids:
        SeqIO.write(record, filtered_fasta, "fastq")
    if p_bar == True:
	i += 1
	printProgress(i, l, prefix = 'Progress:', suffix = 'Complete', barLength = 50)    
filtered_fasta.close()

print "\n\n ### Filtered file is here: %s\n" % outfile_path
