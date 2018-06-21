#!/bin/bash


pop=/home/generic/Christoph/pairwise/sub_2/populationmap_min20000_pstacks_min4_ind_per_pop_FROM_MBUNA
stacks_dir=/media/generic/DATA/Cichlids/stacks_data/christoph-temp
threads=10
min_depth=5

mkdir stacks_data
cd stacks_data
for sample in $(cut -f 1 $pop); do echo $stacks_dir/$sample |perl -ne 'chomp;$full_file="$_".".snps.tsv.gz";if (-e $full_file){`ln -sv $full_file .`;}else{print "WARNING: $_.snps.tsv.gz does not exist\n";}'; done
for sample in $(cut -f 1 $pop); do echo $stacks_dir/$sample |perl -ne 'chomp;$full_file="$_".".matches.tsv.gz";if (-e $full_file){`ln -sv $full_file .`;}else{print "WARNING: $_.matches.tsv.gz does not exist\n";}'; done
for sample in $(cut -f 1 $pop); do echo $stacks_dir/$sample |perl -ne 'chomp;$full_file="$_".".tags.tsv.gz";if (-e $full_file){`ln -sv $full_file .`;}else{print "WARNING: $_.tags.tsv.gz does not exist\n";}'; done

ln -sv $stacks_dir/batch_1.catalog.* .

cd ..

cp $pop intermediate
p=2

for population1 in $(cut -f 2 $pop| sort |uniq)
do
	current=$population1
#	echo -e "current: $current"
#	echo -e "populationmap is $pop"
	cat intermediate | grep $current -v > intermediate2	#create new populationmap excluding the current population
	cp intermediate2 intermediate	#create backup of the second populationmap
	for population2 in $(cat intermediate | cut -f 2 | grep "$current" -v| sort |uniq)	
	do
		echo -e "\n==============\n$population1-$population2"
		mkdir $population1-$population2
		cd $population1-$population2
		grep "$population1" $pop > populationmap
		grep "$population2" $pop >> populationmap
		
		for r in $(seq 0.8 0.1 0.8)
		do
			mkdir r_$r-p_$p
			cd r_$r-p_$p
			date
#			cmd="valgrind --tool=massif --pages-as-heap=yes --massif-out-file=massif.out populations -b 1 -P ../../stacks_data -M ../populationmap -t $threads -r $r -p $p -m $min_depth -f p_value --fstats -k --fasta --genomic --vcf --structure --plink --phylip --write_single_snp --bootstrap_fst"
			cmd="populations -b 1 -P ../../stacks_data -M ../populationmap -t $threads -r $r -p $p -m $min_depth -f p_value --fstats -k --fasta --genomic --vcf --structure --plink --phylip --write_single_snp --bootstrap_fst"
#			cmd="/home/chrishah/src/STACKS/stacks-1.24/bin/bin/populations -b 1 -P ../../stacks_data -M ../populationmap -t $threads -r $r -p $p -m $min_depth -f p_value --fstats -k --fasta --genomic --vcf --structure --plink --phylip --bootstrap"
			echo -e "\n$cmd\n"
			$cmd
			for file in $(ls -1 ../../stacks_data/batch_1.* | grep "catalog" -v); do mv -v $file .; done
			gzip *.*
			cd ..
		done
		cd ..
	done
done
rm -v intermediate*
date
