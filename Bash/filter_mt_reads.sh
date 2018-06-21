#!/bin/bash

query_fasta=$1
genome_fasta=$2
filtered_fasta_name=$3


printf "\nUsage:\n\nfilter_mt_reads.sh <query fasta> <genome fasta> <filtered fastq prefix>\n"



printf "\n 1. Indexing Genome #########\n"
bwa index $2  ## Index genome


printf "\n\n 2. Mapping reads #########\n"
bwa aln -t 8 $2 $1 >  BWA_outs.sai ## Outputs in .sai format


printf "\n\n 3. Convertin output to .sam\n"
bwa samse $2  BWA_outs.sai $1   > BWA_outs.sam  ## Convert to sam

printf "\n\n 4. Getting unmapped reads\n"
samtools view -f 4 BWA_outs.sam -S |  cut -f1 > mt_mapped_reads.txt  ## Get read Ids for reads that didn't map


printf "\n\n 5. Outputting filtered fasta file\n"

gunzip $1 -c > unzipped.temp

for i in $(cat mt_mapped_reads.txt); do grep -A 3 $i unzipped.temp ; done > $3.fq ## Get all reads which didn't map.

rm unzipped.temp

printf "\n\n########## Done ############\n"


