def Repeat_extractor(infile):
    
    """
    
    Extracts the repeats from PB reads masked by DAZZLER

    USAGE: 
         Rep_Extractor <full/input/path>

    OUTPUTS:
         FASTA file called "Extracted_Repeats.fa.gz", containing the sequence of the repeats named 
         by read and location within the read (in input directory).

    """


    import gzip
    
    if infile.endswith("gz"):
    	TEs_in = gzip.open(infile, 'r')
    else:
        TEs_in = open(infile, 'r')
    
    outpath = "%s/%s" % (infile.rpartition("/")[0], "Extracted_Repeats.fa.gz")
    
    TEs_out = gzip.open(outpath, 'w')

    N_lines = 0

    print "\nWorking . . . \n"

    for line in TEs_in:

        N_lines += 1

        if line.startswith(">prolog"):

            if N_lines > 1:

                if repeats == True:

                    for repeat in Repeat_coords:

                        TEs_out.write("%s_%s_%s\n" % (read_ID, Repeat_coords[repeat]["strt"], Repeat_coords[repeat]["end"]))
                        TEs_out.write("%s\n" % whole_seq[int(Repeat_coords[repeat]["strt"]):int(Repeat_coords[repeat]["end"])])

                        #print "%s_%s_%s" % (read_ID, Repeat_coords[repeat]["strt"], Repeat_coords[repeat]["end"])
                        #print whole_seq[int(Repeat_coords[repeat]["strt"]):int(Repeat_coords[repeat]["end"])]
                    #print Repeat_coords
                    #print whole_seq
                    #print N_lines

            repeats = False

            read_ID = line.strip()

            seq = []


        elif line.startswith("> rep"):


            repeats = True

            Repeat_coords = {}

            N_reps = 0

            for rep in line.split():

                if rep.startswith("["):
                    N_reps += 1

                    Repeat_coords[N_reps] = {}
                    Repeat_coords[N_reps]["strt"] = rep.split("[")[1].split(",")[0]
                    Repeat_coords[N_reps]["end"] = rep.split("]")[0].split(",")[1]


        elif all([repeats == True, not line.startswith(">")]):

            New_read = False

            seq.append(line.strip())
            whole_seq = "".join(seq)


    TEs_out.close()

    print "Done . . . Repeats are here: %s\n" % outpath


### CLINE
import sys

if len(sys.argv) < 2:
    print Repeat_extractor.__doc__
    sys.exit("ERROR: Need an input file (and nothing else)\n")

input_file = sys.argv[1]

Repeat_extractor(input_file)


