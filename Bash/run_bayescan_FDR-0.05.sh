#!/bin/bash

man_file=/media/chrishah/STORAGE/RAD/popgen/Fst-outlier/Diplotaxodon/Bayescan/manfile
threads=7
pr_odds=10
FDR="0.05"
Full_pop_map="/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/pop_codes.txt"
stacks_dir="/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/Pure_cru_M_2/"
current_dir=$(pwd)

uniq_pops=$(cut -f1 $Full_pop_map |sort | uniq)

spid=/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/BAYESCAN/Pure_cru_M2_r07_p17_Ho_filtered/STRUCTURE_to_GESTE_BAYE_SCAN.spid

for $pop1 in $uniq_pops
do
	for $pop2 in $uniq_pops
	do
		if [$pop1 != $pop2]
		then
			test_ID=$pop1-$pop2 
			mkdir $test_ID
			python Pairwise_pop_map_maker.py ./master_pop_codes.txt ./$test_ID $pop1 $pop2
			## Makes pop_codes for each pairwise comp and puts in its folder
			echo "Running populations for $pop1 and $pop2"

			## Now Run stacks for the pop_maps that I have just made!
			populations -b1 -P $stacks_dir -r0.7 -p2 -M $test_ID/pop_codes.txt --vcf --structure --write_single_snp
			mv $stacks_dir/batch_1* $test_ID/
			mv $test_ID/batch_1.catalog* $stacks_dir
			echo "Populations has finished for $pop1 and $pop2"

			structure=/$test_ID/batch_1.structure.tsv
			
			cd $test_ID
	        	if [ ! -e "$test_ID.structure" ]
        		then
                		echo -e "fetching structure data from $structure\n"
                		cat $structure | grep "#" -v > $test_ID.structure ## use zcat if gzipped stacks outputs
        		fi
        		if [ ! -e "$test_ID.bayescan" ]
        		then
				echo -e "producing bayescan file from structure file\n"
                		java -Xmx1024m -Xms512m -jar ~/RAD_programs/PGDSpider_2.0.4.0/PGDSpider_2.0.4.0/PGDSpider2-cli.jar -inputfile $test_ID.structure -outputfile $test_ID.bayescan -spid $spid
        		fi
			if [ ! -d "prior_odds_$pr_odds" ]
		        then
                		echo -e "\nrunning Bayescan (prior odds for neutral model: $pr_odds)\n"
	                	mkdir prior_odds_$pr_odds
        	        	cd prior_odds_$pr_odds
               	 		cmd="BayeScan2.1 -snp -od . -pr_odds $pr_odds -threads $threads ../$test_ID.bayescan"
                		echo -e "$cmd\n"
                		$cmd
                		cd ..
			fi
			current=$pwd
			if [ ! -e "prior_odds_$pr_odds/$test_ID-$pr_odds-FDR-$FDR-outlier_bayescan.list" ]
        		then
                		echo -e "\nattempting outlier identification with FDR: $FDR\n"
		                cd prior_odds_$pr_odds
		                echo -e "
source(\"~/RAD_programs/BayeScan2.1/R functions//plot_R.r\")
setwd(\"$current/prior_odds_$pr_odds\")

result<-plot_bayescan('$test_ID.baye_fst.txt',FDR=$FDR)
write(result\$outliers, file = \"$test_ID-$pr_odds-FDR-$FDR-outlier_bayescan.list\", ncolumns = 1)" > R-FDR-$FDR.script

				Rscript R-FDR-$FDR.script
		                cd ..
                		declare -i size
		                size=$(wc prior_odds_$pr_odds/$test_ID-$pr_odds-FDR-$FDR-outlier_bayescan.list | perl -ne 'chomp; @a=split(" "); print "$a[2]\n";')
				if [ "$size" -gt 1 ]
		                then
                		        echo -e "\nBayescan identified candidate outlier loci\n"
                        		for locus in $(cat prior_odds_$pr_odds/$test_ID-$pr_odds-FDR-$FDR-outlier_bayescan.list)
                        		do
                                		declare -i locus
                                		locus=$((locus+2))
                                		head -n 1 $test_ID.structure |cut -f $locus | tr '_' '\t' | cut -f 1 >> prior_odds_$pr_odds/temp
                        		done
                        		cat prior_odds_$pr_odds/temp | sort -n |uniq > $test_ID-$pr_odds-FDR-$FDR.outlier_stacks_ID.list
                		else
                        		echo -e "\nno outliers found with neutral model prior odds $pr_odds\n"
                		fi
        		fi
        		cd ..
		fi

	done
done




