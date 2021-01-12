from Bio import SeqIO
import sys

genome_path = sys.argv[1]

genome = SeqIO.parse(open(genome_path, 'r'), 'fasta')
out_prefix = genome_path.rpartition(".")[0]


for seq in genome:

    for Chr_index in range(1,21,1):
        
        if "scaffold%s," % Chr_index in seq.description:
            chr_file = open("%s_%s.fasta" % (out_prefix, Chr_index), 'w')   
            SeqIO.write(seq, chr_file, 'fasta')
            chr_file.close()
            print "Chr %s written to %s" % (Chr_index,"%s_%s.fasta" % (out_prefix, Chr_index))

print "Done"
